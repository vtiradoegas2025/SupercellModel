/**
 * @file slab.cpp
 * @brief Implementation for the boundary_layer module.
 *
 * Provides executable logic for the boundary_layer runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/boundary_layer subsystem.
 */

#include "slab.hpp"
#include <iostream>
#include <algorithm>

/**
 * @brief Initializes the slab boundary layer scheme.
 */

SlabScheme::SlabScheme()
    : h_(500.0), theta_m_(300.0), qv_m_(0.01),
      entrainment_coeff_(0.2), min_h_(100.0), max_h_(2000.0) 
      {}

/**
 * @brief Returns the required fields for the slab boundary layer scheme.
 */
int SlabScheme::required_fields() const 
{
    return static_cast<int>(BoundaryLayerRequirements::BASIC) |
           static_cast<int>(BoundaryLayerRequirements::MOISTURE);
}

/**
 * @brief Initializes the slab boundary layer scheme.
 */
void SlabScheme::initialize(const BoundaryLayerConfig& cfg) 
{
    min_h_ = 100.0;
    max_h_ = cfg.pbl_max_height;
    entrainment_coeff_ = 0.2;

    std::cout << "Initialized Slab Boundary Layer:" << std::endl;
    std::cout << "  Max PBL height: " << max_h_ << " m" << std::endl;
    std::cout << "  Entrainment coefficient: " << entrainment_coeff_ << std::endl;
}

/**
 * @brief Computes the column for the slab boundary layer scheme.
 */
void SlabScheme::compute_column(
    const BoundaryLayerConfig& cfg,
    const SurfaceConfig& sfc,
    const BoundaryLayerColumnStateView& col,
    BoundaryLayerTendencies& tend,
    BoundaryLayerDiagnostics* diag_opt
) 
{
    const size_t nz = col.theta->size();

    tend.dudt_pbl.assign(nz, 0.0);
    tend.dvdt_pbl.assign(nz, 0.0);
    tend.dthetadt_pbl.assign(nz, 0.0);
    tend.dqvdt_pbl.assign(nz, 0.0);
    tend.dtkedt_pbl.assign(nz, 0.0);

    surface_fluxes::BulkFluxes fluxes;

    if (cfg.apply_surface_fluxes) 
    {
        fluxes = surface_fluxes::compute_surface_fluxes(
            cfg,
            sfc, col.u_sfc, col.v_sfc, col.theta_sfc, col.qv_sfc,
            (*col.theta)[0], (*col.qv)[0], (*col.p)[0], (*col.rho)[0], col.z_sfc
        );
    }

    update_slab_state(col, sfc, fluxes, cfg.dt_pbl);

    apply_slab_tendencies(col, tend, cfg);

    if (diag_opt) 
    {
        diag_opt->pbl_height = h_;
        diag_opt->ustar = fluxes.ustar;
        diag_opt->tau_u = fluxes.tau_u;
        diag_opt->tau_v = fluxes.tau_v;
        diag_opt->H = fluxes.H;
        diag_opt->E = fluxes.E;
        diag_opt->Cd = fluxes.Cd;
        diag_opt->Ch = fluxes.Ch;
        diag_opt->Ce = fluxes.Ce;

        diag_opt->K_m.assign(nz + 1, 0.0);
        diag_opt->K_h.assign(nz + 1, 0.0);
        diag_opt->K_e.assign(nz + 1, 0.0);
    }
}

/**
 * @brief Updates the slab state.
 */
void SlabScheme::update_slab_state(
    const BoundaryLayerColumnStateView& col,
    const SurfaceConfig& sfc,
    const surface_fluxes::BulkFluxes& fluxes,
    double dt
) 
{

    const double cp = boundary_layer_constants::cp;
    const double g = boundary_layer_constants::g;

    double surface_theta_tend = fluxes.H / (cp * (*col.rho)[0]) / h_;

    double theta_above = theta_m_;
    double qv_above = qv_m_;

    if (h_ < (*col.z_int)[col.theta->size()]) 
    {
        for (size_t k = 0; k < col.theta->size(); ++k) 
        {
            if ((*col.z_int)[k] > h_) 
            {
                theta_above = (*col.theta)[k];
                qv_above = (*col.qv)[k];
                break;
            }
        }
    }

    double delta_theta = theta_above - theta_m_;
    double entrainment_flux = entrainment_coeff_ * fluxes.ustar * delta_theta;
    double entrainment_theta_tend = entrainment_flux / h_;

    theta_m_ += (surface_theta_tend + entrainment_theta_tend) * dt;
    qv_m_ += (fluxes.E / (*col.rho)[0] / h_) * dt;

    h_ = diagnose_pbl_height(col, BoundaryLayerConfig{});

    h_ = std::max(min_h_, std::min(max_h_, h_));
    theta_m_ = std::max(250.0, theta_m_);
    qv_m_ = std::max(0.0, qv_m_);
}

/**
 * @brief Diagnoses the PBL height.
 */
double SlabScheme::diagnose_pbl_height(
    const BoundaryLayerColumnStateView& col,
    const BoundaryLayerConfig& cfg
) 
{
    const double Ri_crit = 0.25;

    for (size_t k = 1; k < col.theta->size(); ++k) 
    {
        double z_level = 0.5 * ((*col.z_int)[k] + (*col.z_int)[k-1]);
        double Ri_b = surface_fluxes::bulk_richardson_number(
            theta_m_, (*col.theta)[k-1],
            col.u_sfc, col.v_sfc, z_level
        );

        if (Ri_b > Ri_crit) 
        {
            return z_level;
        }
    }

    return max_h_;
}

/**
 * @brief Applies slab-mixed-layer tendencies to column prognostic fields.
 */
void SlabScheme::apply_slab_tendencies(
    const BoundaryLayerColumnStateView& col,
    BoundaryLayerTendencies& tend,
    const BoundaryLayerConfig& cfg
) 

{
    const size_t nz = col.theta->size();

    for (size_t k = 0; k < nz; ++k) 
    {
        double z_level = 0.5 * ((*col.z_int)[k] + (*col.z_int)[k + 1]);

        if (z_level <= h_) 
        {
            const double tau_relax = 1.0 / 300.0;

            tend.dthetadt_pbl[k] = tau_relax * (theta_m_ - (*col.theta)[k]);
            tend.dqvdt_pbl[k] = tau_relax * (qv_m_ - (*col.qv)[k]);

            tend.dudt_pbl[k] = -0.0001 * (*col.u)[k];
            tend.dvdt_pbl[k] = -0.0001 * (*col.v)[k];
        }
    }
}
