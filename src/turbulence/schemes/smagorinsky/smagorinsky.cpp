#include "smagorinsky.hpp"
#include <iostream>
#include <algorithm>

/*This file contains the implementation of the Smagorinsky turbulence scheme.
It manages the initialization of the Smagorinsky turbulence scheme and the computation of the Smagorinsky turbulence scheme.*/

/*This function initializes the Smagorinsky turbulence scheme.
Takes in the required fields and returns the required fields.*/
SmagorinskyScheme::SmagorinskyScheme()
    : Cs_(turbulence_constants::C_s_default),
      Pr_t_(turbulence_constants::Pr_t_default),
      Sc_t_(turbulence_constants::Sc_t_default),
      nu_t_max_(1000.0),
      K_max_(1000.0) {
}

/*This function gets the required fields.
Takes in the required fields and returns the required fields.*/
int SmagorinskyScheme::required_fields() const 
{
    return static_cast<int>(TurbulenceRequirements::BASIC);
}

/*This function initializes the Smagorinsky turbulence scheme.
Takes in the configuration and initializes the Smagorinsky turbulence scheme.*/
void SmagorinskyScheme::initialize(const TurbulenceConfig& cfg) 
{
    Cs_ = cfg.Cs;
    Pr_t_ = cfg.Pr_t;
    Sc_t_ = cfg.Sc_t;
    nu_t_max_ = cfg.nu_t_max;
    K_max_ = cfg.K_max;

    std::cout << "Initialized Smagorinsky Turbulence:" << std::endl;
    std::cout << "  Cs = " << Cs_ << std::endl;
    std::cout << "  Pr_t = " << Pr_t_ << ", Sc_t = " << Sc_t_ << std::endl;
    std::cout << "  Max nu_t = " << nu_t_max_ << " m²/s" << std::endl;
}

/*This function computes the Smagorinsky turbulence scheme.
Takes in the configuration, grid, state, tendencies, and diagnostics and computes the Smagorinsky turbulence scheme.*/
void SmagorinskyScheme::compute(
    const TurbulenceConfig& cfg,
    const GridMetrics& grid,
    const TurbulenceStateView& state,
    TurbulenceTendencies& tend,
    TurbulenceDiagnostics* diag_opt
) 
{
    // Initialize tendency arrays
    tend.dudt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);
    tend.dvdt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);
    tend.dwdt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);
    tend.dthetadt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);
    tend.dqvdt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);
    tend.dtkedt_sgs.resize(state.NR, state.NTH, state.NZ, 0.0f);

    // Compute eddy coefficients
    compute_eddy_coefficients(cfg, grid, state);

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

    // If the diagnostics are requested, fill the diagnostics.
    if (diag_opt) 
    {
        diag_opt->nu_t = nu_t_;
        diag_opt->K_theta = K_theta_;
        diag_opt->K_q = K_q_;
        diag_opt->K_tke.resize(state.NR, state.NTH, state.NZ, 0.0f);
        diag_opt->Cs_eff.resize(state.NR, state.NTH, state.NZ, static_cast<float>(Cs_));
    }
}

/*This function computes the eddy coefficients.
Takes in the configuration, grid, and state and computes the eddy coefficients.*/
void SmagorinskyScheme::compute_eddy_coefficients(
    const TurbulenceConfig& cfg,
    const GridMetrics& grid,
    const TurbulenceStateView& state
) 
{
    // Initialize coefficient arrays
    nu_t_.resize(state.NR, state.NTH, state.NZ, 0.0f);
    K_theta_.resize(state.NR, state.NTH, state.NZ, 0.0f);
    K_q_.resize(state.NR, state.NTH, state.NZ, 0.0f);

    // Iterate over the rows, columns, and levels and compute the eddy coefficients.
    for (int i = 0; i < state.NR; ++i)
    {
        // Iterate over the columns and compute the eddy coefficients.
        for (int j = 0; j < state.NTH; ++j) 
        {
            // Iterate over the levels and compute the eddy coefficients.
            for (int k = 0; k < state.NZ; ++k) 
            {
                // Compute the strain rate.
                // Compute strain rate
                eddy_viscosity::StrainRate S = eddy_viscosity::compute_strain_rate_3d(state, grid, i, j, k);

                // Filter width
                double Delta = eddy_viscosity::compute_filter_width(cfg, grid, i, j, k);

                // Stability correction
                double stability_factor = apply_stability_correction(cfg, state, i, j, k);

                // Eddy viscosity
                double nu_t = eddy_viscosity::compute_smagorinsky_viscosity(
                    Cs_, Delta, S.magnitude, stability_factor
                );

                // Apply limits
                nu_t = std::min(nu_t, nu_t_max_);

                // Store eddy viscosity
                nu_t_[i][j][k] = static_cast<float>(nu_t);

                // Compute eddy diffusivities
                double K_theta, K_q, K_tke;
                eddy_viscosity::compute_eddy_diffusivities(nu_t, Pr_t_, Sc_t_, K_theta, K_q, K_tke);

                // Apply limits
                K_theta = std::min(K_theta, K_max_);
                K_q = std::min(K_q, K_max_);

                // Store diffusivities
                K_theta_[i][j][k] = static_cast<float>(K_theta);
                K_q_[i][j][k] = static_cast<float>(K_q);
            }
        }
    }
}

/*This function applies the stability correction.
Takes in the configuration, state, and the row, column, and level and applies the stability correction.*/
double SmagorinskyScheme::apply_stability_correction(
    const TurbulenceConfig& cfg,
    const TurbulenceStateView& state,
    int i, int j, int k
) 
{
    // If the stability correction is "ri", apply the stability correction.
    if (cfg.stability_correction == "ri") 
    {
        // Compute local Richardson number (simplified)
        double N = eddy_viscosity::compute_brunt_vaisala_frequency(state, i, j, k);
        double S_mag = eddy_viscosity::compute_strain_rate_3d(state, state.NR > 0 ? GridMetrics{} : GridMetrics{}, i, j, k).magnitude;

        double Ri = (S_mag > 1e-10) ? (N * N) / (S_mag * S_mag) : 0.0;

        return eddy_viscosity::stability_correction_ri(Ri);
    }
    else 
    {
        // No stability correction
        return 1.0;
    }
}
