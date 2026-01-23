#pragma once
#include "../../base/eddy_viscosity.hpp"
#include "turbulence_base.hpp"

/*This file contains the declaration of the TKE turbulence scheme.
It manages the initialization of the TKE turbulence scheme and the computation of the TKE turbulence scheme.*/
class TKEScheme : public TurbulenceSchemeBase 
{
private:
    // TKE scheme parameters (Deardorff-style)
    double Ce1_ = 1.0;        // dissipation coefficient
    double Ce2_ = 1.33;       // dissipation coefficient
    double c_k_ = 0.1;        // eddy viscosity coefficient
    double c_eps_ = 0.19;     // dissipation coefficient
    double Pr_t_ = 0.7;       // turbulent Prandtl number
    double Sc_t_ = 0.7;       // turbulent Schmidt number

    // Mixing length parameters
    double c_l_ = 0.15;       // mixing length coefficient
    double l_max_ = 500.0;    // maximum mixing length [m]

    // Prognostic TKE field
    std::vector<std::vector<std::vector<float>>> tke_;  // TKE [m²/s²]

    // Computed fields
    std::vector<std::vector<std::vector<float>>> nu_t_;     // eddy viscosity
    std::vector<std::vector<std::vector<float>>> K_theta_;  // temperature diffusivity
    std::vector<std::vector<std::vector<float>>> K_q_;      // moisture diffusivity
    std::vector<std::vector<std::vector<float>>> K_tke_;    // TKE diffusivity

public:
    TKEScheme();

    std::string name() const override { return "tke"; }
    int required_fields() const override;

    /*This function initializes the TKE turbulence scheme.
    Takes in the configuration and initializes the TKE turbulence scheme.*/
    void initialize(const TurbulenceConfig& cfg) override;

    /*This function computes the TKE turbulence scheme.
    Takes in the configuration, grid, state, tendencies, and diagnostics and computes the TKE turbulence scheme.*/
    void compute(
        const TurbulenceConfig& cfg,
        const GridMetrics& grid,
        const TurbulenceStateView& state,
        TurbulenceTendencies& tend,
        TurbulenceDiagnostics* diag_opt = nullptr) override;

private:
    
    /*This function computes the mixing length.
    Takes in the configuration, grid, TKE, Brunt-Väisälä frequency, and the row, column, and level and computes the mixing length.*/
    double compute_mixing_length(
        const TurbulenceConfig& cfg,
        const GridMetrics& grid,
        double e,              // TKE
        double N,              // Brunt-Väisälä frequency
        int i, int j, int k
    );

    /*This function computes the eddy coefficients from the TKE.
    Takes in the configuration, grid, and state and computes the eddy coefficients from the TKE.*/
    void compute_eddy_coefficients_from_tke(
        const TurbulenceConfig& cfg,
        const GridMetrics& grid,
        const TurbulenceStateView& state
    );

    /*This function updates the TKE prognostic equation.
    Takes in the configuration, grid, state, and the TKE tendency and updates the TKE prognostic equation.*/
    void update_tke_prognostic(
        const TurbulenceConfig& cfg,
        const GridMetrics& grid,
        const TurbulenceStateView& state,
        std::vector<std::vector<std::vector<float>>>& dtke_dt
    );

    /*This function computes the TKE production.
    Takes in the state, grid, the eddy viscosity, and the shear and buoyancy production and computes the TKE production.*/
    void compute_tke_production(
        const TurbulenceStateView& state,
        const GridMetrics& grid,
        const std::vector<std::vector<std::vector<float>>>& nu_t,
        std::vector<std::vector<std::vector<float>>>& shear_prod,
        std::vector<std::vector<std::vector<float>>>& buoy_prod
    );
};
