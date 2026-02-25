/**
 * @file radiative_transfer.cpp
 * @brief Implementation for the radiation module.
 *
 * Provides executable logic for the radiation runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/radiation subsystem.
 */

#include "radiative_transfer.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>


namespace radiative_transfer 
{

/**
 * @brief Computes the heating rate from the flux divergence.
 */
void heating_rate_from_flux_divergence(
    const std::vector<double>& rho,
    const std::vector<double>& dz,
    const std::vector<double>& Fnet,
    std::vector<double>& dTdt
) 
{
    const size_t nz = rho.size();
    dTdt.resize(nz);
    if (dz.size() != nz || Fnet.size() != nz + 1)
    {
        std::fill(dTdt.begin(), dTdt.end(), 0.0);
        return;
    }

    const double cp = 1004.0;

    for (size_t k = 0; k < nz; ++k) 
    {
        if (rho[k] <= 0.0 || !std::isfinite(rho[k]) || dz[k] <= 0.0 || !std::isfinite(dz[k]))
        {
            dTdt[k] = 0.0;
            continue;
        }
        double dFdz = (Fnet[k+1] - Fnet[k]) / dz[k];
        dTdt[k] = -dFdz / (rho[k] * cp);
        if (!std::isfinite(dTdt[k]))
        {
            dTdt[k] = 0.0;
        }
    }
}

/**
 * @brief Computes the grey longwave two stream solution.
 */
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
    if (tau_layers.size() != nz)
    {
        std::fill(Fup.begin(), Fup.end(), 0.0);
        std::fill(Fdn.begin(), Fdn.end(), 0.0);
        return;
    }

    Fdn[nz] = 0.0;
    Fup[nz] = 0.0;

    Fdn[0] = 0.0;
    Fup[0] = emissivity * stefan_boltzmann(T_sfc) +
             (1.0 - emissivity) * Fdn[0];


  
    for (int k = static_cast<int>(nz) - 1; k >= 0; --k) 
    {
        double tau_k = std::max(0.0, tau_layers[k]);
        double B_k = planck_function(T[k]);


        double exp_tau = exp(-tau_k);
        Fdn[k] = Fdn[k+1] * exp_tau + radiation_constants::pi * B_k * (1.0 - exp_tau);
    }

    for (size_t k = 0; k < nz; ++k) 
    {
        double tau_k = std::max(0.0, tau_layers[k]);
        double B_k = planck_function(T[k]);

        double exp_tau = exp(-tau_k);
        Fup[k+1] = Fup[k] * exp_tau + radiation_constants::pi * B_k * (1.0 - exp_tau);
    }
}

/**
 * @brief Computes the grey shortwave two stream solution.
 */
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
    if (nz == 0)
    {
        Fup[0] = 0.0;
        Fdn[0] = 0.0;
        return;
    }

    mu0 = std::max(mu0, 1e-6);

    Fdn[nz] = mu0 * S0;
    Fup[nz] = 0.0;

    for (int k = static_cast<int>(nz) - 1; k >= 0; --k) 
    {
        const double tau_k = std::max(0.0, tau_layers[k]);
        Fdn[k] = Fdn[k + 1] * exp(-tau_k / mu0);
    }

    Fup[0] = albedo_sfc * Fdn[0];

    for (size_t k = 0; k < nz; ++k) 
    {
        const double tau_k = std::max(0.0, tau_layers[k]);
        Fup[k + 1] = Fup[k] * exp(-tau_k / mu0);
    }
}

/**
 * @brief Computes the optical depth profile.
 */
std::vector<double> compute_optical_depth_profile(
    const std::vector<double>& p,
    double tau_ref,
    double n
) 
{
    const size_t nz = p.size();
    std::vector<double> tau_profile(nz + 1, 0.0);
    if (nz == 0 || !std::isfinite(tau_ref) || tau_ref <= 0.0 || !std::isfinite(n) || n <= 0.0)
    {
        return tau_profile;
    }

    const double p_ref = 100000.0;
    std::vector<double> p_interfaces(nz + 1, 0.0);
    p_interfaces[0] = std::max(0.0, p[0]);

    for (size_t k = 1; k < nz; ++k)
    {
        p_interfaces[k] = std::max(0.0, 0.5 * (p[k - 1] + p[k]));
    }
    p_interfaces[nz] = 0.0;

    for (size_t k = 0; k <= nz; ++k)
    {
        const double pressure_ratio = std::max(0.0, p_interfaces[k] / p_ref);
        tau_profile[k] = tau_ref * std::pow(pressure_ratio, n);
    }
    for (size_t k = 1; k <= nz; ++k)
    {
        tau_profile[k] = std::min(tau_profile[k], tau_profile[k - 1]);
    }

    return tau_profile;
}

/**
 * @brief Computes the layer optical depths.
 */
std::vector<double> compute_layer_optical_depths(
    const std::vector<double>& tau_profile,
    const std::vector<double>& dz
) 
{
    const size_t nz = dz.size();
    std::vector<double> dtau(nz, 0.0);
    if (tau_profile.size() != nz + 1)
    {
        return dtau;
    }

    for (size_t k = 0; k < nz; ++k) 
    {
        dtau[k] = std::max(0.0, tau_profile[k] - tau_profile[k + 1]);
    }

    return dtau;
}

/**
 * @brief Checks the energy conservation.
 */
double check_energy_conservation(
    const std::vector<double>& rho,
    const std::vector<double>& dz,
    const std::vector<double>& dTdt,
    const std::vector<double>& Fnet
) 

{
    const size_t nz = rho.size();
    if (dz.size() != nz || dTdt.size() != nz || Fnet.size() != nz + 1)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double cp = 1004.0;

    double integrated_heating = 0.0;

    for (size_t k = 0; k < nz; ++k) 
    {
        integrated_heating += rho[k] * cp * dTdt[k] * dz[k];
    }

    double net_flux_diff = Fnet[nz] - Fnet[0];

    return integrated_heating - net_flux_diff;
}

}
