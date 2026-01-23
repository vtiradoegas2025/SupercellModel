#include "radiative_transfer.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>

/*This file contains the implementation of the radiative transfer scheme.*/

namespace radiative_transfer 
{

/*This function computes the heating rate from the flux divergence.
Takes in the density, the vertical grid spacing, the net flux, and the heating rate and computes the heating rate from the flux divergence.*/
void heating_rate_from_flux_divergence(
    const std::vector<double>& rho,
    const std::vector<double>& dz,
    const std::vector<double>& Fnet,
    std::vector<double>& dTdt
) 
{
    const size_t nz = rho.size();
    dTdt.resize(nz);

    // Use dry air specific heat capacity (simplified)
    const double cp = 1004.0;  // J/kg/K

    // Iterate over the vertical levels and compute the heating rate from the flux divergence.
    for (size_t k = 0; k < nz; ++k) 
    {
        // Central difference for flux divergence
        double dFdz = (Fnet[k+1] - Fnet[k]) / dz[k];
        // Heating rate: -1/(ρ cp) * dF/dz
        dTdt[k] = -dFdz / (rho[k] * cp);
    }
}

/*This function computes the grey longwave two stream solution.
Takes in the temperature, the optical depth layers, the surface temperature, 
the emissivity, the upward flux, and the downward flux and computes the
grey longwave two stream solution.*/
void grey_lw_two_stream(
    const std::vector<double>& T,
    const std::vector<double>& tau_layers,
    double T_sfc,
    double emissivity,
    std::vector<double>& Fup,
    std::vector<double>& Fdn
) 
{
    const size_t nz = T.size();
    Fup.resize(nz + 1);
    Fdn.resize(nz + 1);

    // Boundary conditions at top (TOA)
    Fdn[nz] = 0.0;  // no downward LW from space
    Fup[nz] = 0.0;  // initially zero, will be computed

    // Boundary conditions at surface
    double B_sfc = planck_function(T_sfc);
    Fdn[0] = 0.0;  // no downward flux at surface
    Fup[0] = emissivity * stefan_boltzmann(T_sfc) +
             (1.0 - emissivity) * Fdn[0];

    // Two-stream solution (simplified grey gas)
    // Work from top down for downward flux, bottom up for upward flux

  
    // Iterate over the vertical levels and compute the downward flux.
    for (int k = static_cast<int>(nz) - 1; k >= 0; --k) 
    {
        double tau_k = tau_layers[k];
        double B_k = planck_function(T[k]);

        // Simple two-stream approximation
        // dF↓/dτ = -F↓ + πB
        // Solution: F↓(τ) = F↓(0) * exp(-τ) + ∫ πB exp(-(τ-τ')) dτ'
        // Approximate as: F↓(k) ≈ F↓(k+1) * exp(-τ_k) + πB_k * (1 - exp(-τ_k))

        double exp_tau = exp(-tau_k);
        Fdn[k] = Fdn[k+1] * exp_tau + radiation_constants::pi * B_k * (1.0 - exp_tau);
    }

    // Iterate over the vertical levels and compute the upward flux.
    for (size_t k = 0; k < nz; ++k) 
    {
        double tau_k = tau_layers[k];
        double B_k = planck_function(T[k]);

        // F↑(k+1) = F↑(k) * exp(-τ_k) + πB_k * (1 - exp(-τ_k))
        double exp_tau = exp(-tau_k);
        Fup[k+1] = Fup[k] * exp_tau + radiation_constants::pi * B_k * (1.0 - exp_tau);
    }
}

/*This function computes the grey shortwave two stream solution.
Takes in the optical depth layers, the cosine of the zenith angle, the solar constant, 
the surface albedo, the upward flux, and the downward flux and computes the
grey shortwave two stream solution.*/
void grey_sw_beer_lambert(
    const std::vector<double>& tau_layers,
    double mu0,
    double S0,
    double albedo_sfc,
    std::vector<double>& Fup,
    std::vector<double>& Fdn
) 
{
    const size_t nz = tau_layers.size();
    Fup.resize(nz + 1);
    Fdn.resize(nz + 1);

    // Ensure mu0 > 0 (avoid division by zero)
    mu0 = std::max(mu0, 1e-6);

    // TOA downward flux
    Fdn[nz] = mu0 * S0;
    Fup[nz] = 0.0;

    // Surface boundary (reflection)
    double cum_tau = 0.0;  // cumulative optical depth from top

    // Iterate over the vertical levels and compute the downward flux.
    for (int k = static_cast<int>(nz) - 1; k >= 0; --k) 
    {
        cum_tau += tau_layers[k];
        // Beer-Lambert: F↓(z) = F↓(TOA) * exp(-τ/μ₀)
        Fdn[k] = Fdn[nz] * exp(-cum_tau / mu0);
    }

    // Surface reflection
    Fup[0] = albedo_sfc * Fdn[0];

    // Upward attenuation (simplified - no scattering)
    cum_tau = 0.0;

    // Iterate over the vertical levels and compute the upward flux.
    for (size_t k = 0; k < nz; ++k) 
    {
        cum_tau += tau_layers[k];
        // Simple upward attenuation (no absorption in upward path)
        Fup[k+1] = Fup[k] * exp(-cum_tau / mu0);
    }
}

/*This function computes the optical depth profile.
Takes in the pressure, the reference pressure, and the optical depth and computes the optical depth profile.*/
std::vector<double> compute_optical_depth_profile(
    const std::vector<double>& p,
    double tau_ref,
    double n
) 
{
    const size_t nz = p.size();
    std::vector<double> tau_profile(nz + 1);

    const double p_ref = 100000.0;  // reference pressure (Pa)

    // Optical depth at interfaces (assume linear in pressure)
    tau_profile[0] = 0.0;  // surface

    // Iterate over the vertical levels and compute the optical depth profile.
    for (size_t k = 1; k <= nz; ++k) 
    {
        double p_avg = (k < nz) ? 0.5 * (p[k-1] + p[k]) : p[nz-1];
        tau_profile[k] = tau_ref * pow(p_avg / p_ref, n);
    }

    return tau_profile;
}

/*This function computes the layer optical depths.
Takes in the optical depth profile, and the vertical grid spacing and computes the layer optical depths.*/
std::vector<double> compute_layer_optical_depths(
    const std::vector<double>& tau_profile,
    const std::vector<double>& dz
) 
{
    const size_t nz = dz.size();
    std::vector<double> dtau(nz);

    // Iterate over the vertical levels and compute the layer optical depths.
    for (size_t k = 0; k < nz; ++k) 
    {
        // Average optical depth across layer
        dtau[k] = 0.5 * (tau_profile[k] + tau_profile[k+1]);
    }

    return dtau;
}

/*This function checks the energy conservation.
Takes in the density, the vertical grid spacing, the heating rate, 
and the net flux and checks the energy conservation.*/
double check_energy_conservation(
    const std::vector<double>& rho,
    const std::vector<double>& dz,
    const std::vector<double>& dTdt,
    const std::vector<double>& Fnet
) 

{
    const size_t nz = rho.size();
    const double cp = 1004.0;  // J/kg/K

    // Column-integrated heating
    double integrated_heating = 0.0;

    // Iterate over the vertical levels and compute the column-integrated heating.
    for (size_t k = 0; k < nz; ++k) 
    {
        integrated_heating += rho[k] * cp * dTdt[k] * dz[k];
    }

    // Net flux difference (should equal integrated heating)
    double net_flux_diff = Fnet[nz] - Fnet[0];  // TOA - surface

    return integrated_heating - net_flux_diff;  // should be ~0
}

} // namespace radiative_transfer
