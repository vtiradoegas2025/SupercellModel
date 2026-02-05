#include "tke.hpp"
#include <iostream>
#include <algorithm>


/*This file contains the implementation of the TKE turbulence scheme.
It manages the initialization of the TKE turbulence scheme and the computation of the TKE turbulence scheme.*/


/*This function initializes the TKE turbulence scheme.
Takes in the required fields and returns the required fields.*/
TKEScheme::TKEScheme()
    : Ce1_(1.0), Ce2_(1.33), c_k_(0.1), c_eps_(0.19),
      Pr_t_(0.7), Sc_t_(0.7), c_l_(0.15), l_max_(500.0) {
}

/*This function gets the required fields.
Takes in the required fields and returns the required fields.*/
int TKEScheme::required_fields() const 
{
    return static_cast<int>(TurbulenceRequirements::BASIC) |
           static_cast<int>(TurbulenceRequirements::TKE);
}

/*This function initializes the TKE turbulence scheme.
Takes in the configuration and initializes the TKE turbulence scheme.*/ 
void TKEScheme::initialize(const TurbulenceConfig& cfg) 
{
    Ce1_ = 1.0;
    Ce2_ = 1.33;
    c_k_ = 0.1;
    c_eps_ = 0.19;
    Pr_t_ = cfg.Pr_t;
    Sc_t_ = cfg.Sc_t;
    c_l_ = 0.15;
    l_max_ = 500.0;

    std::cout << "Initialized TKE Turbulence (1.5-order closure):" << std::endl;
    std::cout << "  Prognostic TKE: enabled" << std::endl;
    std::cout << "  Ce1 = " << Ce1_ << ", Ce2 = " << Ce2_ << std::endl;
    std::cout << "  Max mixing length = " << l_max_ << " m" << std::endl;
}

/*This function computes the TKE turbulence scheme.
Takes in the configuration, grid, state, tendencies, and diagnostics and computes the TKE turbulence scheme.*/
void TKEScheme::compute(
    const TurbulenceConfig& cfg,
    const GridMetrics& grid,
    const TurbulenceStateView& state,
    TurbulenceTendencies& tend,
    TurbulenceDiagnostics* diag_opt
) 
{
    // If the TKE field is not set, initialize the TKE field.
    if (tke_.empty()) 
    {
        tke_.resize(state.NR, state.NTH, state.NZ, 0.1f); // background TKE
    }

    // Initialize tendency arrays
    tend.dudt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);
    tend.dvdt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);
    tend.dwdt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);
    tend.dthetadt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);
    tend.dqvdt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);
    tend.dtkedt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);

    // Compute eddy coefficients from current TKE
    compute_eddy_coefficients_from_tke(cfg, grid, state);

    // Apply momentum diffusion
    eddy_viscosity::compute_momentum_diffusion_tendencies(
        state, grid, nu_t_, tend.dudt_sgs, tend.dvdt_sgs, tend.dwdt_sgs
    );

    // Apply scalar diffusion (temperature)
    eddy_viscosity::compute_scalar_diffusion_tendencies(
        state, grid, K_theta_, *state.theta, tend.dthetadt_sgs
    );

    // If the moisture is not set, return.
    if (state.qv) 
    {
        eddy_viscosity::compute_scalar_diffusion_tendencies(
            state, grid, K_q_, *state.qv, tend.dqvdt_sgs
        );
    }

    // Update prognostic TKE
    update_tke_prognostic(cfg, grid, state, tend.dtkedt_sgs);

    // If the diagnostics are requested, fill the diagnostics.  
    if (diag_opt) 
    {
        diag_opt->nu_t = nu_t_;
        diag_opt->K_theta = K_theta_;
        diag_opt->K_q = K_q_;
        diag_opt->K_tke = K_tke_;

        // Compute TKE budget terms for diagnostics
        Field3D shear_prod, buoy_prod;
        shear_prod.resize(state.NR, state.NTH, state.NZ, 0.0f);
        buoy_prod.resize(state.NR, state.NTH, state.NZ, 0.0f);
        diag_opt->dissipation.resize(state.NR, state.NTH, state.NZ, 0.0f);

        compute_tke_production(state, grid, nu_t_, shear_prod, buoy_prod);

        diag_opt->shear_prod = shear_prod;
        diag_opt->buoy_prod = buoy_prod;

        // Iterate over the rows, columns, and levels and compute the dissipation.
        for (int i = 0; i < state.NR; ++i) 
        {

            for (int j = 0; j < state.NTH; ++j)
            {

                for (int k = 0; k < state.NZ; ++k) 
                {
                    double e = static_cast<double>(tke_[i][j][k]);
                    double l = compute_mixing_length(cfg, grid, e,
                        eddy_viscosity::compute_brunt_vaisala_frequency(state, i, j, k), i, j, k);
                    double eps = c_eps_ * std::pow(e, 1.5) / l;
                    diag_opt->dissipation[i][j][k] = static_cast<float>(eps);
                }
            }
        }
    }
}


/*This function computes the mixing length.
Takes in the configuration, grid, TKE, Brunt-Väisälä frequency, and the row, column, and level and computes the mixing length.*/
double TKEScheme::compute_mixing_length(
    const TurbulenceConfig& cfg,
    const GridMetrics& grid,
    double e,
    double N,
    int i, int j, int k
) 
{
    // Filter width (mixing length limit)
    double Delta = eddy_viscosity::compute_filter_width(cfg, grid, i, j, k);

    // Stability-limited mixing length
    double l_stability = (N > 1e-6) ? c_l_ * std::sqrt(e) / N : Delta;

    // Take minimum of grid-based and stability limits
    double l = std::min(Delta, l_stability);

    // Apply maximum limit
    return std::min(l, l_max_);
}

/*This function computes the eddy coefficients from the TKE.
Takes in the configuration, grid, and state and computes the eddy coefficients from the TKE.*/
void TKEScheme::compute_eddy_coefficients_from_tke(
    const TurbulenceConfig& cfg,
    const GridMetrics& grid,
    const TurbulenceStateView& state
) {
    // Initialize coefficient arrays
    nu_t_.resize(state.NR, state.NTH, state.NZ, 0.0f);
    K_theta_.resize(state.NR, state.NTH, state.NZ, 0.0f);
    K_q_.resize(state.NR, state.NTH, state.NZ, 0.0f);
    K_tke_.resize(state.NR, state.NTH, state.NZ, 0.0f);

    // Iterate over the rows, columns, and levels and compute the eddy coefficients from the TKE.
    for (int i = 0; i < state.NR; ++i) 
    {
        // Iterate over the columns and compute the eddy coefficients from the TKE.
        for (int j = 0; j < state.NTH; ++j) 
        {
            // Iterate over the levels and compute the eddy coefficients from the TKE.
            for (int k = 0; k < state.NZ; ++k) 
            {
                // Compute the TKE.
                double e = static_cast<double>(tke_[i][j][k]);
                double N = eddy_viscosity::compute_brunt_vaisala_frequency(state, i, j, k);

                // Mixing length
                double l = compute_mixing_length(cfg, grid, e, N, i, j, k);

                // Eddy viscosity: nu_t = c_k * l * sqrt(e)
                double nu_t = c_k_ * l * std::sqrt(std::max(e, 1e-6));

                // Apply limits
                nu_t = std::min(nu_t, static_cast<double>(cfg.nu_t_max));

                // Store eddy viscosity
                nu_t_[i][j][k] = static_cast<float>(nu_t);

                // Compute eddy diffusivities
                double K_theta, K_q, K_tke;
                eddy_viscosity::compute_eddy_diffusivities(nu_t, Pr_t_, Sc_t_, K_theta, K_q, K_tke);

                // Apply limits
                K_theta = std::min(K_theta, static_cast<double>(cfg.K_max));
                K_q = std::min(K_q, static_cast<double>(cfg.K_max));
                K_tke = std::min(K_tke, static_cast<double>(cfg.K_max));

                // Store diffusivities
                K_theta_[i][j][k] = static_cast<float>(K_theta);
                K_q_[i][j][k] = static_cast<float>(K_q);
                K_tke_[i][j][k] = static_cast<float>(K_tke);
            }
        }
    }
}

/*This function updates the TKE prognostic equation.
Takes in the configuration, grid, state, and the TKE tendency and updates the TKE prognostic equation.*/
void TKEScheme::update_tke_prognostic(
    const TurbulenceConfig& cfg,
    const GridMetrics& grid,
    const TurbulenceStateView& state,
    Field3D& dtke_dt
) 
{
    // Initialize TKE tendency
    dtke_dt.resize(state.NR, state.NTH, state.NZ, 0.0f);

    // Compute TKE production terms
    Field3D shear_prod, buoy_prod;
    shear_prod.resize(state.NR, state.NTH, state.NZ, 0.0f);
    buoy_prod.resize(state.NR, state.NTH, state.NZ, 0.0f);

    compute_tke_production(state, grid, nu_t_, shear_prod, buoy_prod);

    // Iterate over the rows, columns, and levels and update the TKE prognostic equation.
    for (int i = 0; i < state.NR; ++i) 
    {
        // Iterate over the columns and update the TKE prognostic equation.
        for (int j = 0; j < state.NTH; ++j) 
        {
            // Iterate over the levels and update the TKE prognostic equation.
            for (int k = 0; k < state.NZ; ++k) 
            {
                double e = static_cast<double>(tke_[i][j][k]);
                double N = eddy_viscosity::compute_brunt_vaisala_frequency(state, i, j, k);
                double l = compute_mixing_length(cfg, grid, e, N, i, j, k);

                // Shear production
                double P_s = static_cast<double>(shear_prod[i][j][k]);

                // Buoyancy production
                double P_b = static_cast<double>(buoy_prod[i][j][k]);

                // Dissipation: ε = Cε * (e^{3/2}) / l
                double eps = Ce2_ * std::pow(std::max(e, 1e-6), 1.5) / l;

                // TKE diffusion (simplified - would need proper flux divergence)
                double de_dz = 0.0;  // simplified
                double diffusion = static_cast<double>(K_tke_[i][j][k]) * de_dz / 100.0; // assume 100m scale

                // Total tendency: d(e)/dt = P_s + P_b - ε + diffusion
                double dtke = P_s + P_b - eps + diffusion;

                dtke_dt[i][j][k] = static_cast<float>(dtke);
            }
        }
    }
}

/*This function computes the TKE production.
Takes in the state, grid, the eddy viscosity, and the shear and buoyancy production and computes the TKE production.*/
void TKEScheme::compute_tke_production(
    const TurbulenceStateView& state,
    const GridMetrics& grid,
    const Field3D& nu_t,
    Field3D& shear_prod,
    Field3D& buoy_prod
) 
{
    // Iterate over the rows, columns, and levels and compute the TKE production.
    for (int i = 0; i < state.NR; ++i) 
    {
        // Iterate over the columns and compute the TKE production.
        for (int j = 0; j < state.NTH; ++j) 
        {
            // Iterate over the levels and compute the TKE production.
            for (int k = 0; k < state.NZ; ++k) 
            {
                // Shear production: P_s = 2 * nu_t * |S|²
                eddy_viscosity::StrainRate S = eddy_viscosity::compute_strain_rate_3d(
                    state, grid, i, j, k);
                double nu_local = static_cast<double>(nu_t[i][j][k]);
                double P_s = 2.0 * nu_local * S.magnitude * S.magnitude;

                shear_prod[i][j][k] = static_cast<float>(P_s);

                // Buoyancy production: P_b = -nu_t/Pr_t * (g/θ) * dθ/dz

                // If the level is not the last level, compute the buoyancy production.
                if (k < state.NZ - 1) 
                {
                    double theta_k = static_cast<double>((*state.theta)[i][j][k]);
                    double theta_kp1 = static_cast<double>((*state.theta)[i][j][k+1]);
                    double dtheta_dz = (theta_kp1 - theta_k) / 100.0;  // assume 100m layer

                    double P_b = -nu_local / 0.7 * (9.81 / theta_k) * dtheta_dz;
                    buoy_prod[i][j][k] = static_cast<float>(P_b);
                }
            }
        }
    }
}
