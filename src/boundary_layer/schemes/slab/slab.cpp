#include "slab.hpp"
#include <iostream>
#include <algorithm>

/*This constructor initializes the slab boundary layer scheme. 
Takes in the initial PBL height, initial potential temperature, 
initial water vapor mixing ratio, entrainment coefficient, 
minimum PBL height, and maximum PBL height and initializes 
the slab boundary layer scheme.*/

SlabScheme::SlabScheme()
    : h_(500.0), theta_m_(300.0), qv_m_(0.01),
      entrainment_coeff_(0.2), min_h_(100.0), max_h_(2000.0) 
      {}

/*This function returns the required fields for the slab boundary layer scheme.*/
int SlabScheme::required_fields() const 
{
    return static_cast<int>(BoundaryLayerRequirements::BASIC) |
           static_cast<int>(BoundaryLayerRequirements::MOISTURE);
}

/*This function initializes the slab boundary layer scheme. Takes in 
the boundary layer configurationand initializes the slab boundary 
layer scheme.*/
void SlabScheme::initialize(const BoundaryLayerConfig& cfg) 
{
    min_h_ = 100.0;
    max_h_ = cfg.pbl_max_height;
    entrainment_coeff_ = 0.2;  // tunable parameter

    std::cout << "Initialized Slab Boundary Layer:" << std::endl;
    std::cout << "  Max PBL height: " << max_h_ << " m" << std::endl;
    std::cout << "  Entrainment coefficient: " << entrainment_coeff_ << std::endl;
}

/*This function computes the column for the slab boundary layer scheme. 
Takes in the boundary layer configuration, surface configuration, column state, 
and column tendencies and computes the boundary layer physics.
Passes out the boundary layer tendencies and diagnostics to the calling function 
for use in the simulation.*/
void SlabScheme::compute_column(
    const BoundaryLayerConfig& cfg,
    const SurfaceConfig& sfc,
    const BoundaryLayerColumnStateView& col,
    BoundaryLayerTendencies& tend,
    BoundaryLayerDiagnostics* diag_opt
) 
{
    const size_t nz = col.theta->size();

    // Resize tendency arrays
    tend.dudt_pbl.assign(nz, 0.0);
    tend.dvdt_pbl.assign(nz, 0.0);
    tend.dthetadt_pbl.assign(nz, 0.0);
    tend.dqvdt_pbl.assign(nz, 0.0);
    tend.dtkedt_pbl.assign(nz, 0.0);  // not used in slab

    // Compute surface fluxes if enabled
    surface_fluxes::BulkFluxes fluxes;

    // If the surface fluxes are to be applied, compute the surface fluxes.
    if (cfg.apply_surface_fluxes) 
    {
        fluxes = surface_fluxes::compute_bulk_fluxes(
            sfc, col.u_sfc, col.v_sfc, col.theta_sfc, col.qv_sfc,
            (*col.theta)[0], (*col.qv)[0], (*col.p)[0], (*col.rho)[0], col.z_sfc
        );
    }

    // Update slab state
    update_slab_state(col, sfc, fluxes, cfg.dt_pbl);

    // Apply tendencies to grid
    apply_slab_tendencies(col, tend, cfg);

    // If the diagnostics are to be computed, compute the diagnostics.
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

        // Simple eddy diffusivity profile
        diag_opt->K_m.assign(nz + 1, 0.0);
        diag_opt->K_h.assign(nz + 1, 0.0);
        diag_opt->K_e.assign(nz + 1, 0.0);
    }
}

/*This function updates the slab state. 
Takes in the column state, the surface configuration, the surface fluxes, and the time step
and updates the slab state.*/
void SlabScheme::update_slab_state(
    const BoundaryLayerColumnStateView& col,
    const SurfaceConfig& sfc,
    const surface_fluxes::BulkFluxes& fluxes,
    double dt
) 
{
    // Simple slab model evolution
    // Heat budget: d(theta_m)/dt = (surface_flux + entrainment_flux) / h

    const double cp = boundary_layer_constants::cp;
    const double g = boundary_layer_constants::g;

    // Surface heat flux tendency [K/s]
    double surface_theta_tend = fluxes.H / (cp * (*col.rho)[0]) / h_;

    // Simple entrainment parameterization
    // Find inversion strength above PBL
    double theta_above = theta_m_;
    double qv_above = qv_m_;

    // If the PBL height is less than the height of the column, find the level above the PBL.
    if (h_ < (*col.z_int)[col.theta->size()]) 
    {
        // Iterate over the vertical levels.
        for (size_t k = 0; k < col.theta->size(); ++k) 
        {
            // If the level is above the PBL, set the potential temperature and water vapor mixing ratio to the level above the PBL.
            if ((*col.z_int)[k] > h_) 
            {
                theta_above = (*col.theta)[k];
                qv_above = (*col.qv)[k];
                break;
            }
        }
    }

    // Entrainment flux (simplified for now). Feel free to change/build upon this.
    double delta_theta = theta_above - theta_m_;
    double entrainment_flux = entrainment_coeff_ * fluxes.ustar * delta_theta;
    double entrainment_theta_tend = entrainment_flux / h_;

    // Update mixed-layer properties
    theta_m_ += (surface_theta_tend + entrainment_theta_tend) * dt;
    qv_m_ += (fluxes.E / (*col.rho)[0] / h_) * dt;  // moisture tendency

    // Update PBL height (simple diagnostic)
    h_ = diagnose_pbl_height(col, BoundaryLayerConfig{});

    // Ensure bounds
    h_ = std::max(min_h_, std::min(max_h_, h_));
    theta_m_ = std::max(250.0, theta_m_);
    qv_m_ = std::max(0.0, qv_m_);
}

/*This function diagnoses the PBL height. 
Takes in the column state and boundary layer configuration and diagnoses the PBL height.
Passes out the PBL height to the calling function for use in the simulation.*/
double SlabScheme::diagnose_pbl_height(
    const BoundaryLayerColumnStateView& col,
    const BoundaryLayerConfig& cfg
) 
{
    // Simple PBL height diagnostic based on bulk Ri
    const double Ri_crit = 0.25;  // critical Richardson number

    // Iterate over the vertical levels 
    for (size_t k = 1; k < col.theta->size(); ++k) 
    {
        double z_level = 0.5 * ((*col.z_int)[k] + (*col.z_int)[k-1]);
        double Ri_b = surface_fluxes::bulk_richardson_number(
            theta_m_, (*col.theta)[k-1],
            col.u_sfc, col.v_sfc, z_level
        );

        // If the bulk Richardson number is greater than the critical Richardson number, return the level height.
        if (Ri_b > Ri_crit) 
        {
            return z_level;
        }
    }

    return max_h_;  // if no inversion found
}

void SlabScheme::apply_slab_tendencies(
    const BoundaryLayerColumnStateView& col,
    BoundaryLayerTendencies& tend,
    const BoundaryLayerConfig& cfg
) 

{
    const size_t nz = col.theta->size();

    // Iterate over the vertical levels.
    for (size_t k = 0; k < nz; ++k) 
    {
        // Get the level height
        double z_level = 0.5 * ((*col.z_int)[k] + (*col.z_int)[k-1]);

        // If the level height is less than or equal to the PBL height, apply the slab tendencies.
        if (z_level <= h_) 
        {
            // Within PBL: relax toward slab values
            const double tau_relax = 1.0 / 300.0;  // relaxation timescale [1/s]

            tend.dthetadt_pbl[k] = tau_relax * (theta_m_ - (*col.theta)[k]);
            tend.dqvdt_pbl[k] = tau_relax * (qv_m_ - (*col.qv)[k]);

            // Simple momentum mixing (reduced diffusion)
            tend.dudt_pbl[k] = -0.0001 * (*col.u)[k];  // weak damping
            tend.dvdt_pbl[k] = -0.0001 * (*col.v)[k];
        }
        // Above PBL: no tendencies (stable stratification)
    }
}
