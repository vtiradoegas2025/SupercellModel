#include "ysu.hpp"
#include <iostream>
#include <algorithm>

YSUScheme::YSUScheme()
    : pblfac_(1.0), cn_(0.75), ck_(0.1), ce_(0.5), c0_(0.15) {
}

/*This function returns the required fields for the YSU boundary layer scheme.*/
int YSUScheme::required_fields() const 
{
    return static_cast<int>(BoundaryLayerRequirements::BASIC) |
           static_cast<int>(BoundaryLayerRequirements::MOISTURE);
}


/*This function initializes the YSU boundary layer scheme. Takes in the boundary layer configuration
and initializes the YSU boundary layer scheme.*/
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

/*This function computes the column for the YSU boundary layer scheme. 
Takes in the boundary layer configuration, surface configuration, column state, 
and column tendencies and computes the boundary layer physics.
Passes out the boundary layer tendencies and diagnostics to the calling function 
for use in the simulation.*/
void YSUScheme::compute_column(
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
    tend.dtkedt_pbl.assign(nz, 0.0);  // not used in YSU

    // Compute surface fluxes
    surface_fluxes::BulkFluxes fluxes;

    // If the surface fluxes are to be applied, compute the surface fluxes.
    if (cfg.apply_surface_fluxes) 
    {
        fluxes = surface_fluxes::compute_bulk_fluxes(
            sfc, col.u_sfc, col.v_sfc, col.theta_sfc, col.qv_sfc,
            (*col.theta)[0], (*col.qv)[0], (*col.p)[0], (*col.rho)[0], col.z_sfc
        );
    }

    // Diagnose PBL height
    double h = diagnose_pbl_height(col, cfg);

    // Compute K-profile diffusivities
    K_m_.resize(nz + 1);
    K_h_.resize(nz + 1);
    compute_k_profile(col, h, fluxes.ustar, K_m_, K_h_);

    // Apply diffusion tendencies
    apply_diffusion_tendencies(col, K_m_, K_h_, tend);

    // If the nonlocal transport is enabled, compute the nonlocal transport.
    if (cfg.enable_nonlocal_transport) 
    {
        std::vector<double> nonlocal_theta(nz, 0.0);
        std::vector<double> nonlocal_qv(nz, 0.0);
        compute_nonlocal_transport(col, fluxes, h, nonlocal_theta, nonlocal_qv);

        // Iterate over the vertical levels.
        for (size_t k = 0; k < nz; ++k) 
        {
            tend.dthetadt_pbl[k] += nonlocal_theta[k];
            tend.dqvdt_pbl[k] += nonlocal_qv[k];
        }
    }

    // If the diagnostics are to be computed, compute the diagnostics.
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
        diag_opt->K_e.assign(nz + 1, 0.0);  // not used in YSU
    }
}

/*This function diagnoses the PBL height. 
Takes in the column state and boundary layer configuration and diagnoses the PBL height.
Passes out the PBL height to the calling function for use in the simulation.*/
double YSUScheme::diagnose_pbl_height(
    const BoundaryLayerColumnStateView& col,
    const BoundaryLayerConfig& cfg
) 

{
    // YSU PBL height diagnosis based on critical Richardson number
    const double Ri_crit = 0.25;
    const double min_h = 100.0;
    const double max_h = cfg.pbl_max_height;

    // Iterate over the vertical levels.
    for (size_t k = 1; k < col.theta->size(); ++k) 
    {
        // Get the level height
        double z_level = 0.5 * ((*col.z_int)[k] + (*col.z_int)[k-1]);
        double Ri_b = surface_fluxes::bulk_richardson_number(
            (*col.theta)[0], (*col.theta)[k],
            (*col.u)[0], (*col.v)[0], z_level
        );

        // If the bulk Richardson number is greater than the critical Richardson number, return the level height.
        if (Ri_b > Ri_crit) 
        {
            return std::max(min_h, std::min(max_h, z_level));
        }
    }

    return max_h;
}

/*This function computes the K-profile diffusivities. 
Takes in the column state, the PBL height, the surface friction velocity, and the K-profile diffusivities
and computes the K-profile diffusivities.*/
void YSUScheme::compute_k_profile(
    const BoundaryLayerColumnStateView& col,
    double h, double ustar,
    std::vector<double>& K_m,
    std::vector<double>& K_h
) 

{
    // Initialize the K-profile diffusivities structure
    const size_t nz = col.theta->size();
    const double kappa = boundary_layer_constants::kappa;

    // Surface values
    K_m[0] = c0_ * ustar * boundary_layer_constants::kappa * col.z_sfc;
    K_h[0] = K_m[0];

    // Iterate over the vertical levels.
    for (size_t k = 1; k <= nz; ++k) 
    {
        // Get the level height
        double z_level = (*col.z_int)[k-1];  // interface height

        // If the level height is less than or equal to the PBL height, compute the K-profile diffusivities.
        if (z_level <= h) 
        {
            // Within PBL: K increases with height
            double phi_m = z_level / h;
            double phi_h = phi_m;

            // Stability correction (simplified)
            double Ri_z = 0.0;  // could compute local Ri

            // If the bulk Richardson number is greater than 0, apply the stability correction.
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
            // Above PBL: small constant diffusivity
            K_m[k] = 1.0;  // m²/s
            K_h[k] = 1.0;
        }

        // Apply maximum limits
        K_m[k] = std::min(K_m[k], 100.0);
        K_h[k] = std::min(K_h[k], 100.0);
    }
}


/*This function computes the nonlocal transport. 
Takes in the column state, the surface fluxes, the PBL height, 
the nonlocal theta, and the nonlocal qv
and computes the nonlocal transport.*/
void YSUScheme::compute_nonlocal_transport(
    const BoundaryLayerColumnStateView& col,
    const surface_fluxes::BulkFluxes& fluxes,
    double h,
    std::vector<double>& nonlocal_theta,
    std::vector<double>& nonlocal_qv
) 

{
    // Initialize the nonlocal transport structure
    const size_t nz = col.theta->size();

    // Simplified nonlocal transport (counter-gradient flux)
    // In YSU, this represents the nonlocal heat transport term


    for (size_t k = 0; k < nz; ++k) 
    {
        double z_level = 0.5 * ((*col.z_int)[k] + (*col.z_int)[k-1]);

        // If the level height is less than or equal to the PBL height, compute the nonlocal transport.
        if (z_level <= h) 
        {
            // Within PBL: nonlocal flux contribution
            double gamma_theta = ck_ * fluxes.H / (boundary_layer_constants::cp * (*col.rho)[k]);
            double gamma_qv = ck_ * fluxes.E / (*col.rho)[k];

            // Counter-gradient term (simplified profile)
            double profile_factor = 1.0 - z_level / h;  // decreases with height

            nonlocal_theta[k] = gamma_theta * profile_factor / h;
            nonlocal_qv[k] = gamma_qv * profile_factor / h;
        }
    }
}

/*This function applies the diffusion tendencies. 
Takes in the column state, the K-profile diffusivities, and the boundary layer tendencies
and applies the diffusion tendencies.*/
void YSUScheme::apply_diffusion_tendencies(
    const BoundaryLayerColumnStateView& col,
    const std::vector<double>& K_m,
    const std::vector<double>& K_h,
    BoundaryLayerTendencies& tend
) {
    const size_t nz = col.theta->size();

    // Compute diffusion tendencies using finite differences
    // d(φ)/dt = d/dz(K dφ/dz)

    // Iterate over the vertical levels.
    for (size_t k = 0; k < nz; ++k) 
    {
        // Vertical gradients (centered differences)
        double dz_k = (*col.z_int)[k+1] - (*col.z_int)[k];

        // Momentum diffusion
        double du_dz = (k > 0 && k < nz-1) ?
            ((*col.u)[k+1] - (*col.u)[k-1]) / (2 * dz_k) :
            ((*col.u)[k+1] - (*col.u)[k]) / dz_k;

        double dv_dz = (k > 0 && k < nz-1) ?
            ((*col.v)[k+1] - (*col.v)[k-1]) / (2 * dz_k) :
            ((*col.v)[k+1] - (*col.v)[k]) / dz_k;

        double K_m_avg = 0.5 * (K_m[k] + K_m[k+1]);
        tend.dudt_pbl[k] = (K_m[k+1] * du_dz - K_m[k] * du_dz) / (dz_k * (*col.rho)[k]);
        tend.dvdt_pbl[k] = (K_m[k+1] * dv_dz - K_m[k] * dv_dz) / (dz_k * (*col.rho)[k]);

        // Scalar diffusion (theta and qv)
        double dtheta_dz = (k > 0 && k < nz-1) ?
            ((*col.theta)[k+1] - (*col.theta)[k-1]) / (2 * dz_k) :
            ((*col.theta)[k+1] - (*col.theta)[k]) / dz_k;

        double dqv_dz = (k > 0 && k < nz-1) ?
            ((*col.qv)[k+1] - (*col.qv)[k-1]) / (2 * dz_k) :
            ((*col.qv)[k+1] - (*col.qv)[k]) / dz_k;

        double K_h_avg = 0.5 * (K_h[k] + K_h[k+1]);
        tend.dthetadt_pbl[k] = (K_h[k+1] * dtheta_dz - K_h[k] * dtheta_dz) / (dz_k * (*col.rho)[k]);
        tend.dqvdt_pbl[k] = (K_h[k+1] * dqv_dz - K_h[k] * dqv_dz) / (dz_k * (*col.rho)[k]);
    }
}
