#pragma once
#include <vector>
#include "../../../include/radiation_base.hpp"

/*This file contains the implementation of the radiative transfer scheme.*/
namespace radiative_transfer 
{

// Compute heating rate from flux divergence
void heating_rate_from_flux_divergence(
    const std::vector<double>& rho,     // density [kg/m³] (nz)
    const std::vector<double>& dz,      // layer thickness [m] (nz)
    const std::vector<double>& Fnet,    // net flux [W/m²] (nz+1)
    std::vector<double>& dTdt          // output: heating rate [K/s] (nz)
);


/*This function computes the planck function.
Takes in the temperature and computes the planck function.*/
inline double planck_function(double T) 
{
    // Grey body emission (integrated over wavelengths)
    // B(T) = (σ/π) T⁴, so πB = σT⁴
    using namespace radiation_constants;
    return sigma * pow(T, 4.0) / pi;
}

/*This function computes the stefan-boltzmann law.
Takes in the temperature and computes the stefan-boltzmann law.*/
inline double stefan_boltzmann(double T) 
{
    using namespace radiation_constants;
    return sigma * pow(T, 4.0);
}

/*This function computes the grey longwave two stream solution.
Takes in the temperature, the optical depth layers, the surface temperature, 
the emissivity, the upward flux, and the downward flux and computes the
grey longwave two stream solution.*/
void grey_lw_two_stream(
    const std::vector<double>& T,           // temperature [K] (nz)
    const std::vector<double>& tau_layers,  // optical depth per layer (nz)
    double T_sfc,                          // surface temperature [K]
    double emissivity,                     // surface emissivity
    std::vector<double>& Fup,              // output: upward flux [W/m²] (nz+1)
    std::vector<double>& Fdn               // output: downward flux [W/m²] (nz+1)
);

/*This function computes the grey shortwave beer-lambert solution.
Takes in the optical depth layers, the cosine of the zenith angle, 
the solar constant, the surface albedo, the upward flux, and the 
downward flux and computes thegrey shortwave beer-lambert solution.*/
void grey_sw_beer_lambert(
    const std::vector<double>& tau_layers,  // optical depth per layer (nz)
    double mu0,                             // cos(solar zenith angle)
    double S0,                              // solar constant [W/m²]
    double albedo_sfc,                      // surface albedo
    std::vector<double>& Fup,               // output: upward flux [W/m²] (nz+1)
    std::vector<double>& Fdn                // output: downward flux [W/m²] (nz+1)
);

/*This function computes the optical depth profile.
Takes in the pressure, the reference pressure, and 
the optical depth and computes the optical depth profile.*/
std::vector<double> compute_optical_depth_profile(
    const std::vector<double>& p,      // pressure [Pa] (nz)
    double tau_ref,                    // reference optical depth
    double n                           // exponent
);

/*This function computes the layer optical depths.
Takes in the optical depth profile, and the vertical
grid spacing and computes the layer optical depths.*/
std::vector<double> compute_layer_optical_depths(
    const std::vector<double>& tau_profile,  // optical depth at interfaces (nz+1)
    const std::vector<double>& dz            // layer thickness [m] (nz)
);

/*This function checks the energy conservation.
Takes in the density, the vertical grid spacing, the heating rate, 
and the net flux and checks the energy conservation.*/
double check_energy_conservation(
    const std::vector<double>& rho,     // density [kg/m³] (nz)
    const std::vector<double>& dz,      // layer thickness [m] (nz)
    const std::vector<double>& dTdt,    // heating rate [K/s] (nz)
    const std::vector<double>& Fnet     // net flux [W/m²] (nz+1)
);

} // namespace radiative_transfer
