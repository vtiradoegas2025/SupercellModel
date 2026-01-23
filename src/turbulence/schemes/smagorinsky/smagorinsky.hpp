#pragma once
#include "../../base/eddy_viscosity.hpp"
#include "turbulence_base.hpp"

/*This file contains the declaration of the Smagorinsky turbulence scheme.
It manages the initialization of the Smagorinsky turbulence scheme and the computation of the Smagorinsky turbulence scheme.*/

class SmagorinskyScheme : public TurbulenceSchemeBase {
private:
    // Smagorinsky scheme parameters
    double Cs_;              // Smagorinsky coefficient
    double Pr_t_;            // turbulent Prandtl number
    double Sc_t_;            // turbulent Schmidt number
    double nu_t_max_;        // maximum eddy viscosity
    double K_max_;           // maximum eddy diffusivity

    // Computed fields
    std::vector<std::vector<std::vector<float>>> nu_t_;     // eddy viscosity
    std::vector<std::vector<std::vector<float>>> K_theta_;  // temperature diffusivity
    std::vector<std::vector<std::vector<float>>> K_q_;      // moisture diffusivity

public:
    /*This function initializes the Smagorinsky turbulence scheme.
    Takes in the required fields and returns the required fields.*/
    SmagorinskyScheme();

    /*This function gets the name of the Smagorinsky turbulence scheme.
    Takes in the name and returns the name of the Smagorinsky turbulence scheme.*/
    std::string name() const override { return "smagorinsky"; }

    /*This function gets the required fields.
    Takes in the required fields and returns the required fields.*/
    int required_fields() const override;

    /*This function initializes the Smagorinsky turbulence scheme.
    Takes in the configuration and initializes the Smagorinsky turbulence scheme.*/
    void initialize(const TurbulenceConfig& cfg) override;

    void compute(
        const TurbulenceConfig& cfg,
        const GridMetrics& grid,
        const TurbulenceStateView& state,
        TurbulenceTendencies& tend,
        TurbulenceDiagnostics* diag_opt = nullptr) override;

private:

    /*This function computes the eddy coefficients.
    Takes in the configuration, grid, and state and computes the eddy coefficients.*/
   
    void compute_eddy_coefficients(
        const TurbulenceConfig& cfg,
        const GridMetrics& grid,
        const TurbulenceStateView& state
    );

    /*This function applies the stability correction.
    Takes in the configuration, state, and the row, column, and level and applies the stability correction.*/
    double apply_stability_correction(
        const TurbulenceConfig& cfg,
        const TurbulenceStateView& state,
        int i, int j, int k
    );
};
