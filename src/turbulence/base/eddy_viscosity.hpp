#pragma once
#include <vector>
#include "../../../include/turbulence_base.hpp"
#include "../../../include/field3d.hpp"

// Eddy viscosity and turbulence utilities for turbulence schemes
namespace eddy_viscosity 
{

/*This struct contains the strain rate tensor.
Takes in the strain rate components and the magnitude and returns the strain rate tensor.*/
struct StrainRate 
{
    double S11, S12, S13;  // strain rate components
    double S21, S22, S23;  // (symmetric tensor)
    double S31, S32, S33;
    double magnitude;      // |S| = sqrt(2*Sij*Sij)
};

/*This function computes the 3D strain rate tensor at a grid point.
Takes in the state, grid, and the row, column, and level and computes the 3D strain rate tensor.*/
StrainRate compute_strain_rate_3d(
    const TurbulenceStateView& state,
    const GridMetrics& grid,
    int i, int j, int k  // grid indices
);

/*This function computes the filter width (Smagorinsky mixing length).
Takes in the configuration, grid, and the row, column, and level and computes the filter width.*/
double compute_filter_width(
    const TurbulenceConfig& cfg,
    const GridMetrics& grid,
    int i, int j, int k  // grid indices (for potential future variations)
);

/*This function computes the Smagorinsky viscosity.
Takes in the Smagorinsky coefficient, filter width, strain magnitude, and stability factor and computes the Smagorinsky viscosity.*/
double compute_smagorinsky_viscosity(
    double Cs,           // Smagorinsky coefficient
    double Delta,        // filter width [m]
    double strain_mag,   // |S| [1/s]
    double stability_factor = 1.0  // stability correction factor
);

/*This function computes the eddy diffusivities from the eddy viscosity.
Takes in the eddy viscosity, turbulent Prandtl number, turbulent Schmidt number, and the eddy diffusivities and computes the eddy diffusivities.*/
void compute_eddy_diffusivities(
    double nu_t,         // eddy viscosity [m²/s]
    double Pr_t,         // turbulent Prandtl number
    double Sc_t,         // turbulent Schmidt number
    double& K_theta,     // output: temperature diffusivity [m²/s]
    double& K_q,         // output: moisture diffusivity [m²/s]
    double& K_tke        // output: TKE diffusivity [m²/s]
);

/*This function computes the scalar diffusion tendency.
Takes in the state, grid, the eddy diffusivities, and the scalar and computes the scalar diffusion tendency.*/
double compute_scalar_diffusion_tendency(
    const TurbulenceStateView& state,
    const GridMetrics& grid,
    const Field3D& K_field,  // diffusivity field
    const Field3D& phi,      // scalar field
    int i, int j, int k,  // grid indices
    int var_index        // 0=u, 1=v, 2=w, 3=theta, 4=qv, 5=tke
);

/*This function computes the momentum diffusion tendencies.
Takes in the state, grid, the eddy viscosity, and the momentum and computes the momentum diffusion tendencies.*/
void compute_momentum_diffusion_tendencies(
    const TurbulenceStateView& state,
    const GridMetrics& grid,
    const Field3D& nu_t,  // eddy viscosity
    Field3D& dudt_sgs,
    Field3D& dvdt_sgs,
    Field3D& dwdt_sgs
);

/*This function computes the scalar diffusion tendencies.
Takes in the state, grid, the eddy diffusivities, and the scalar and computes the scalar diffusion tendencies.*/
void compute_scalar_diffusion_tendencies(
    const TurbulenceStateView& state,
    const GridMetrics& grid,
    const Field3D& K_field,  // diffusivity
    const Field3D& phi,      // scalar field
    Field3D& dphi_dt_sgs     // tendency output
);

/*This function computes the stability correction.
Takes in the Richardson number and the critical Richardson number and computes the stability correction.*/
double stability_correction_ri(
    double Ri,           // Richardson number
    double Ri_crit = 0.25  // critical Richardson number
);

/*This function computes the TKE mixing length.
Takes in the filter width, TKE, Brunt-Väisälä frequency, and the stability factor and computes the TKE mixing length.*/
double compute_tke_mixing_length(
    double Delta,        // filter width [m]
    double e,           // TKE [m²/s²]
    double N,           // Brunt-Väisälä frequency [1/s]
    double c_l = 0.15   // mixing length coefficient
);

/*This function computes the Brunt-Väisälä frequency.
Takes in the state, and the row, column, and level and computes the Brunt-Väisälä frequency.*/
double compute_brunt_vaisala_frequency(
    const TurbulenceStateView& state,
    int i, int j, int k  // grid indices
);

/*This function applies the positivity limits.
Takes in the field, minimum value, and maximum value and applies the positivity limits.*/
void apply_positivity_limits(
    Field3D& field,
    double min_value = 0.0,
    double max_value = 1e6
);

} // namespace eddy_viscosity
