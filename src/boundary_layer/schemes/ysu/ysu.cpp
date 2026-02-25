/**
 * @file ysu.cpp
 * @brief Implementation for the boundary_layer module.
 *
 * Provides executable logic for the boundary_layer runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/boundary_layer subsystem.
 */

#include "ysu.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

YSUScheme::YSUScheme()
    : pblfac_(1.0), cn_(0.75), ck_(0.1), ce_(0.5), c0_(0.15) {
}

/**
 * @brief Returns the required fields for the YSU boundary layer scheme.
 */
int YSUScheme::required_fields() const 
{
    return static_cast<int>(BoundaryLayerRequirements::BASIC) |
           static_cast<int>(BoundaryLayerRequirements::MOISTURE);
}


/**
 * @brief Initializes the YSU boundary layer scheme.
 */
void YSUScheme::initialize(const BoundaryLayerConfig& cfg) 
{
    pblfac_ = 1.0;
    cn_ = 0.75;
    ck_ = cfg.enable_nonlocal_transport ? 0.1 : 0.0;
    ce_ = 0.5;
    c0_ = 0.15;

    std::cout << "Initialized YSU Boundary Layer:" << std::endl;
    std::cout << "  Nonlocal transport: " << (ck_ > 0.0 ? "enabled" : "disabled") << std::endl;
}

/**
 * @brief Computes the column for the YSU boundary layer scheme.
 */
void YSUScheme::compute_column(
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

    double h = diagnose_pbl_height(col, cfg);

    K_m_.resize(nz + 1);
    K_h_.resize(nz + 1);
    compute_k_profile(col, h, fluxes.ustar, K_m_, K_h_);

    apply_diffusion_tendencies(col, K_m_, K_h_, tend);

    if (cfg.enable_nonlocal_transport) 
    {
        std::vector<double> nonlocal_theta(nz, 0.0);
        std::vector<double> nonlocal_qv(nz, 0.0);
        compute_nonlocal_transport(col, fluxes, h, nonlocal_theta, nonlocal_qv);

        for (size_t k = 0; k < nz; ++k) 
        {
            tend.dthetadt_pbl[k] += nonlocal_theta[k];
            tend.dqvdt_pbl[k] += nonlocal_qv[k];
        }
    }

    if (diag_opt) 
    {
        diag_opt->pbl_height = h;
        diag_opt->ustar = fluxes.ustar;
        diag_opt->tau_u = fluxes.tau_u;
        diag_opt->tau_v = fluxes.tau_v;
        diag_opt->H = fluxes.H;
        diag_opt->E = fluxes.E;
        diag_opt->Cd = fluxes.Cd;
        diag_opt->Ch = fluxes.Ch;
        diag_opt->Ce = fluxes.Ce;
        diag_opt->K_m = K_m_;
        diag_opt->K_h = K_h_;
        diag_opt->K_e.assign(nz + 1, 0.0);
    }
}

/**
 * @brief Diagnoses the PBL height.
 */
double YSUScheme::diagnose_pbl_height(
    const BoundaryLayerColumnStateView& col,
    const BoundaryLayerConfig& cfg
) 

{
    const double Ri_crit = 0.25;
    const double min_h = 100.0;
    const double max_h = cfg.pbl_max_height;

    for (size_t k = 1; k < col.theta->size(); ++k) 
    {
        double z_level = 0.5 * ((*col.z_int)[k] + (*col.z_int)[k-1]);
        double Ri_b = surface_fluxes::bulk_richardson_number(
            (*col.theta)[0], (*col.theta)[k],
            (*col.u)[0], (*col.v)[0], z_level
        );

        if (Ri_b > Ri_crit) 
        {
            return std::max(min_h, std::min(max_h, z_level));
        }
    }

    return max_h;
}

/**
 * @brief Computes the K-profile diffusivities.
 */
void YSUScheme::compute_k_profile(
    const BoundaryLayerColumnStateView& col,
    double h, double ustar,
    std::vector<double>& K_m,
    std::vector<double>& K_h
) 

{
    const size_t nz = col.theta->size();
    const double kappa = boundary_layer_constants::kappa;

    K_m[0] = c0_ * ustar * boundary_layer_constants::kappa * col.z_sfc;
    K_h[0] = K_m[0];

    for (size_t k = 1; k <= nz; ++k) 
    {
        double z_level = (*col.z_int)[k-1];

        if (z_level <= h) 
        {
            double phi_m = z_level / h;
            double phi_h = phi_m;

            double Ri_z = 0.0;

            if (Ri_z > 0.0) 
            {
                phi_m = 1.0 + 5.0 * Ri_z;
                phi_h = phi_m;
            }

            K_m[k] = kappa * ustar * z_level / phi_m;
            K_h[k] = kappa * ustar * z_level / phi_h;
        } 
        else 
        {
            K_m[k] = 1.0;
            K_h[k] = 1.0;
        }

        K_m[k] = std::min(K_m[k], 100.0);
        K_h[k] = std::min(K_h[k], 100.0);
    }
}


/**
 * @brief Computes the nonlocal transport.
 */
void YSUScheme::compute_nonlocal_transport(
    const BoundaryLayerColumnStateView& col,
    const surface_fluxes::BulkFluxes& fluxes,
    double h,
    std::vector<double>& nonlocal_theta,
    std::vector<double>& nonlocal_qv
) 

{
    const size_t nz = col.theta->size();
    if (nz == 0) return;
    const double h_safe = std::max(h, 1.0);



    for (size_t k = 0; k < nz; ++k) 
    {
        double z_level = 0.5 * ((*col.z_int)[k] + (*col.z_int)[k + 1]);

        if (z_level <= h) 
        {
            const double rho_k = std::max((*col.rho)[k], 1.0e-6);
            double gamma_theta = ck_ * fluxes.H / (boundary_layer_constants::cp * rho_k);
            double gamma_qv = ck_ * fluxes.E / rho_k;

            double profile_factor = 1.0 - z_level / h_safe;

            nonlocal_theta[k] = gamma_theta * profile_factor / h_safe;
            nonlocal_qv[k] = gamma_qv * profile_factor / h_safe;
        }
    }
}

/**
 * @brief Applies the diffusion tendencies.
 */
void YSUScheme::apply_diffusion_tendencies(
    const BoundaryLayerColumnStateView& col,
    const std::vector<double>& K_m,
    const std::vector<double>& K_h,
    BoundaryLayerTendencies& tend
) {
    const size_t nz = col.theta->size();
    if (nz == 0 || col.z_int->size() < nz + 1) return;

    auto compute_tendency = [&](const std::vector<double>& profile,
                                const std::vector<double>& K,
                                size_t k,
                                double rho_k) -> double
    {
        if (nz < 2) return 0.0;

        const double rho_safe = std::max(rho_k, 1.0e-6);
        const double dz_up = std::max((*col.z_int)[k + 1] - (*col.z_int)[k], 1.0e-6);

        if (k == 0)
        {
            const double flux_up = -K[1] * (profile[1] - profile[0]) / dz_up;
            const double flux_dn = 0.0;
            return (flux_dn - flux_up) / (rho_safe * dz_up);
        }

        const double dz_dn = std::max((*col.z_int)[k] - (*col.z_int)[k - 1], 1.0e-6);
        if (k == nz - 1)
        {
            const double flux_up = 0.0;
            const double flux_dn = -K[k] * (profile[k] - profile[k - 1]) / dz_dn;
            return (flux_dn - flux_up) / (rho_safe * dz_dn);
        }

        const double flux_up = -K[k + 1] * (profile[k + 1] - profile[k]) / dz_up;
        const double flux_dn = -K[k] * (profile[k] - profile[k - 1]) / dz_dn;
        const double dz_cell = 0.5 * (dz_up + dz_dn);
        return (flux_dn - flux_up) / (rho_safe * std::max(dz_cell, 1.0e-6));
    };


    for (size_t k = 0; k < nz; ++k) 
    {
        const double rho_k = (*col.rho)[k];

        tend.dudt_pbl[k] = compute_tendency(*col.u, K_m, k, rho_k);
        tend.dvdt_pbl[k] = compute_tendency(*col.v, K_m, k, rho_k);
        tend.dthetadt_pbl[k] = compute_tendency(*col.theta, K_h, k, rho_k);
        tend.dqvdt_pbl[k] = compute_tendency(*col.qv, K_h, k, rho_k);

        if (!std::isfinite(tend.dudt_pbl[k])) tend.dudt_pbl[k] = 0.0;
        if (!std::isfinite(tend.dvdt_pbl[k])) tend.dvdt_pbl[k] = 0.0;
        if (!std::isfinite(tend.dthetadt_pbl[k])) tend.dthetadt_pbl[k] = 0.0;
        if (!std::isfinite(tend.dqvdt_pbl[k])) tend.dqvdt_pbl[k] = 0.0;
    }
}
