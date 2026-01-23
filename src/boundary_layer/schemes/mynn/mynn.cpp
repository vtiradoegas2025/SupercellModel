#include "mynn.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>

/*This file contains the implementation of the MYNN boundary layer scheme.
This module is used to compute the boundary layer physics of the simulation.
Takes in the boundary layer configuration, surface configuration, column state, and column tendencies
and computes the boundary layer physics.
Passes out the boundary layer tendencies and diagnostics to the calling function for use in the simulation.*/

/*This constructor initializes the MYNN boundary layer scheme.*/
MYNNScheme::MYNNScheme()
    : a1_(1.18), a2_(0.665), b1_(24.0), b2_(15.0),
      c1_(0.137), c2_(0.75), c3_(0.352), c4_(0.0), c5_(0.2),
      ce1_(1.0), ce2_(1.33), lmax_(500.0) {
}

/*This function returns the required fields for the MYNN boundary layer scheme.*/
int MYNNScheme::required_fields() const 
{
    return static_cast<int>(BoundaryLayerRequirements::BASIC) |
           static_cast<int>(BoundaryLayerRequirements::MOISTURE) |
           static_cast<int>(BoundaryLayerRequirements::TKE);
}

/*This function initializes the MYNN boundary layer scheme Takes in the boundary layer configuration
and initializes the MYNN boundary layer scheme.*/
void MYNNScheme::initialize(const BoundaryLayerConfig& cfg) 
{
    lmax_ = cfg.pbl_max_height * 0.5;  // mixing length limit

    std::cout << "Initialized MYNN Boundary Layer:" << std::endl;
    std::cout << "  Prognostic TKE: enabled" << std::endl;
    std::cout << "  Max mixing length: " << lmax_ << " m" << std::endl;
}

/*This function computes the column for the MYNN boundary layer scheme. 
Takes in the boundary layer configuration, surface configuration, column state, 
and column tendencies and computes the boundary layer physics.
Passes out the boundary layer tendencies and diagnostics to the calling function 
for use in the simulation.*/
void MYNNScheme::compute_column(
    const BoundaryLayerConfig& cfg,
    const SurfaceConfig& sfc,
    const BoundaryLayerColumnStateView& col,
    BoundaryLayerTendencies& tend,
    BoundaryLayerDiagnostics* diag_opt
) 
{
    const size_t nz = col.theta->size();

    // If the TKE array is not the same size as the number of vertical levels, 
    // initialize it to a background TKE.
    if (tke_.size() != nz) 
    {
        tke_.assign(nz, 0.1);  // background TKE [m²/s²]
    }

    // Resize tendency arrays
    tend.dudt_pbl.assign(nz, 0.0);
    tend.dvdt_pbl.assign(nz, 0.0);
    tend.dthetadt_pbl.assign(nz, 0.0);
    tend.dqvdt_pbl.assign(nz, 0.0);
    tend.dtkedt_pbl.assign(nz, 0.0);

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

    // Compute mixing length
    std::vector<double> l_mix(nz);
    compute_mixing_length(col, h, l_mix);

    // Compute eddy diffusivities
    std::vector<double> K_m(nz + 1, 0.0);
    std::vector<double> K_h(nz + 1, 0.0);
    std::vector<double> K_e(nz + 1, 0.0);
    compute_eddy_diffusivities(col, l_mix, K_m, K_h, K_e);

    // Update TKE prognostic equation
    update_tke(col, fluxes, K_m, K_h, K_e, cfg.dt_pbl, tend.dtkedt_pbl);

    // Apply diffusion tendencies
    apply_diffusion_tendencies(col, K_m, K_h, tend);

    // Iterate over the vertical levels and update the TKE field.
    for (size_t k = 0; k < nz; ++k) 
    {
        tke_[k] += tend.dtkedt_pbl[k] * cfg.dt_pbl;
        tke_[k] = std::max(0.001, tke_[k]);  // minimum TKE
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
        diag_opt->K_m = K_m;
        diag_opt->K_h = K_h;
        diag_opt->K_e = K_e;
    }
}

/*This function diagnoses the PBL height. 
Takes in the column state and boundary layer configuration 
and diagnoses the PBL height.*/
double MYNNScheme::diagnose_pbl_height(
    const BoundaryLayerColumnStateView& col,
    const BoundaryLayerConfig& cfg
) 
{
    // MYNN PBL height based on TKE minimum or Ri criterion
    const double min_h = 100.0;
    const double max_h = cfg.pbl_max_height;

    // Find level where TKE drops below threshold
    const double tke_min = 0.01;  // minimum TKE threshold [m²/s²]

    // Iterate over the vertical levels.
    for (size_t k = 0; k < col.theta->size(); ++k)
    {
        // If the TKE drops below the threshold, return the level where the TKE dropped below .
        if (tke_[k] < tke_min) 
        {
            double z_level = 0.5 * ((*col.z_int)[k] + (*col.z_int)[k-1]);
            return std::max(min_h, std::min(max_h, z_level));
        }
    }

    return max_h;
}

/*This function computes the mixing length. 
Takes in the column state, the PBL height, and the mixing length 
and computes the mixing length.*/
void MYNNScheme::compute_mixing_length(
    const BoundaryLayerColumnStateView& col,
    double h,
    std::vector<double>& l_mix
) 
{
    // Initialize the mixing length structure
    const size_t nz = col.theta->size();
    const double kappa = boundary_layer_constants::kappa;

    // Iterate over the vertical levels.
    for (size_t k = 0; k < nz; ++k) 
    {
        double z_level = 0.5 * ((*col.z_int)[k] + (*col.z_int)[k-1]);

        // Blackadar mixing length formulation (simplified)
        double l_blackadar = kappa * z_level / (1.0 + kappa * z_level / lmax_);

        // Limit by distance to boundaries
        double l_wall = kappa * z_level;
        double l_top = kappa * (h - z_level);

        l_mix[k] = std::min({l_blackadar, l_wall, l_top, lmax_});
        l_mix[k] = std::max(l_mix[k], 1.0);  // minimum mixing length
    }
}

/*This function computes the eddy diffusivities. 
Takes in the column state, the mixing length, and the eddy diffusivities 
and computes the eddy diffusivities.*/
void MYNNScheme::compute_eddy_diffusivities(
    const BoundaryLayerColumnStateView& col,
    const std::vector<double>& l_mix,
    std::vector<double>& K_m,
    std::vector<double>& K_h,
    std::vector<double>& K_e
) {
    const size_t nz = col.theta->size();

    // Iterate over the vertical levels.
    for (size_t k = 1; k <= nz; ++k) 
    {
        // Get the cell index
        size_t k_cell = k - 1;  // cell index

        // Get the TKE at the cell center
        double e_k = tke_[k_cell];  // TKE at cell center

        // MYNN stability functions (simplified)
        double q2 = 2.0 * e_k;  // twice TKE
        double gh = 0.0;  // stability parameter (could compute from Ri)

        // Stability functions
        double sm = a1_ * (1.0 - c1_ * gh) / ((1.0 - gh) * (1.0 + a2_ * gh));
        double sh = a1_ * (1.0 - c2_ * gh) / ((1.0 - gh) * (1.0 + a2_ * gh));

        // Limit stability functions
        sm = std::max(sm, 0.01);
        sh = std::max(sh, 0.01);

        // Eddy diffusivities
        double l_k = l_mix[k_cell];
        K_m[k] = l_k * sqrt(q2) * sm;
        K_h[k] = l_k * sqrt(q2) * sh;
        K_e[k] = l_k * sqrt(q2) * sh;  // same as heat for TKE

        // Apply limits
        K_m[k] = std::min(K_m[k], 1000.0);
        K_h[k] = std::min(K_h[k], 1000.0);
        K_e[k] = std::min(K_e[k], 1000.0);
    }

    // Surface values
    K_m[0] = 0.0;  // no diffusivity at surface
    K_h[0] = 0.0;
    K_e[0] = 0.0;
}

/*This function updates the TKE field. 
Takes in the column state, the surface fluxes, the eddy diffusivities, the time step, and the TKE tendency
and updates the TKE field.*/
void MYNNScheme::update_tke(
    const BoundaryLayerColumnStateView& col,
    const surface_fluxes::BulkFluxes& fluxes,
    const std::vector<double>& K_m,
    const std::vector<double>& K_h,
    const std::vector<double>& K_e,
    double dt,
    std::vector<double>& tke_tend
) 
{
    // Initialize the TKE tendency structure
    const size_t nz = col.theta->size();

    // TKE prognostic equation: d(e)/dt = P_s + P_b - ε + ∇·(K_e ∇e)
    // where P_s = shear production, P_b = buoyancy production, ε = dissipation

    for (size_t k = 0; k < nz; ++k) 
    {
        // Shear production: K_m * (du/dz)²
        double dz_k = (*col.z_int)[k+1] - (*col.z_int)[k];
        double du_dz = ((*col.u)[k] - (k > 0 ? (*col.u)[k-1] : 0.0)) / dz_k;
        double dv_dz = ((*col.v)[k] - (k > 0 ? (*col.v)[k-1] : 0.0)) / dz_k;
        double S2 = du_dz * du_dz + dv_dz * dv_dz;  // shear squared

        double P_s = K_m[k+1] * S2;

        // Buoyancy production (simplified)
        double dtheta_dz = ((*col.theta)[k] - (k > 0 ? (*col.theta)[k-1] : (*col.theta)[k])) / dz_k;
        double P_b = -K_h[k+1] * dtheta_dz * boundary_layer_constants::g / (*col.theta)[k];

        // Dissipation: ε = (q²)^{3/2} / (B1 * l)
        double q2 = 2.0 * tke_[k];
        double q3 = q2 * sqrt(q2);
        double l_mix = 100.0;  // simplified mixing length
        double eps = q3 / (b1_ * l_mix);

        // TKE diffusion (simplified central difference)
        double de_dz = (k < nz-1) ? (tke_[k+1] - tke_[k-1]) / (2 * dz_k) : 0.0;
        double dK_dz = K_e[k+1] - K_e[k];
        double diffusion = (K_e[k+1] * de_dz + dK_dz * de_dz) / (*col.rho)[k];

        // Total tendency
        tke_tend[k] = P_s + P_b - eps + diffusion / (*col.rho)[k];

        // Add surface production at lowest level
        if (k == 0) {
            double ustar3 = fluxes.ustar * fluxes.ustar * fluxes.ustar;
            tke_tend[k] += ustar3 / (boundary_layer_constants::kappa * col.z_sfc);
        }
    }
}

/*This function applies the diffusion tendencies. 
Takes in the column state, the eddy diffusivities, and the boundary layer tendencies
and applies the diffusion tendencies.*/
void MYNNScheme::apply_diffusion_tendencies(
    const BoundaryLayerColumnStateView& col,
    const std::vector<double>& K_m,
    const std::vector<double>& K_h,
    BoundaryLayerTendencies& tend
) 
{
    // Initialize the boundary layer tendencies structure
    const size_t nz = col.theta->size();

    // Iterate over the vertical levels.
    for (size_t k = 0; k < nz; ++k) 
    {
        // Get the cell index
        size_t k_cell = k - 1;  // cell index

        double dz_k = (*col.z_int)[k+1] - (*col.z_int)[k];

        // Momentum diffusion
        double K_m_avg = 0.5 * (K_m[k] + K_m[k+1]);
        double du_dz = (k < nz-1) ?
            ((*col.u)[k+1] - (*col.u)[k]) / dz_k : 0.0;
        double dv_dz = (k < nz-1) ?
            ((*col.v)[k+1] - (*col.v)[k]) / dz_k : 0.0;

        tend.dudt_pbl[k] = K_m_avg * du_dz / dz_k / (*col.rho)[k];
        tend.dvdt_pbl[k] = K_m_avg * dv_dz / dz_k / (*col.rho)[k];

        // Scalar diffusion
        double K_h_avg = 0.5 * (K_h[k] + K_h[k+1]);
        double dtheta_dz = (k < nz-1) ?
            ((*col.theta)[k+1] - (*col.theta)[k]) / dz_k : 0.0;
        double dqv_dz = (k < nz-1) ?
            ((*col.qv)[k+1] - (*col.qv)[k]) / dz_k : 0.0;

        tend.dthetadt_pbl[k] = K_h_avg * dtheta_dz / dz_k / (*col.rho)[k];
        tend.dqvdt_pbl[k] = K_h_avg * dqv_dz / dz_k / (*col.rho)[k];
    }
}
