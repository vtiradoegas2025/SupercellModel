#pragma once
#include <vector>
#include "../../../include/boundary_layer_base.hpp"

/*This header file contains the base classes and structures for the surface fluxes module.
The surface fluxes module is responsible for the surface fluxes of the simulation.
The surface fluxes scheme is chosen by the user in the configuration file.
This module is used to compute the surface fluxes of the simulation.*/

// Surface layer and flux utilities for boundary layer schemes
namespace surface_fluxes 
{

// Bulk aerodynamic surface flux computation
struct BulkFluxes 
{
    double ustar = 0.0;     // friction velocity [m/s]
    double tau_u = 0.0;     // zonal momentum stress [Pa]
    double tau_v = 0.0;     // meridional momentum stress [Pa]
    double H = 0.0;         // sensible heat flux [W/m²]
    double E = 0.0;         // moisture flux [kg/m²/s]
    double Cd = 0.0;        // drag coefficient
    double Ch = 0.0;        // heat transfer coefficient
    double Ce = 0.0;        // moisture transfer coefficient
};

// Compute bulk aerodynamic surface fluxes
BulkFluxes compute_bulk_fluxes(
    const SurfaceConfig& sfc_cfg,
    double u_sfc,           // surface zonal wind [m/s]
    double v_sfc,           // surface meridional wind [m/s]
    double theta_sfc,       // surface potential temperature [K]
    double qv_sfc,          // surface water vapor [kg/kg]
    double theta_air,       // air potential temperature [K]
    double qv_air,          // air water vapor [kg/kg]
    double p_air,           // air pressure [Pa]
    double rho_air,         // air density [kg/m³]
    double z_sfc            // height of measurement [m]
);

// Monin-Obukhov similarity theory surface fluxes (WRF-style)
// Based on Jiménez et al. (2012, MWR)
struct MoninObukhovFluxes 
{
    double ustar = 0.0;     // friction velocity [m/s]
    double tau_u = 0.0;     // zonal momentum stress [Pa]
    double tau_v = 0.0;     // meridional momentum stress [Pa]
    double H = 0.0;         // sensible heat flux [W/m²]
    double E = 0.0;         // moisture flux [kg/m²/s]
    double L = 0.0;         // Obukhov length [m]
    double zeta = 0.0;      // stability parameter z/L
    double Cd = 0.0;        // drag coefficient
    double Ch = 0.0;        // heat transfer coefficient
    double Ce = 0.0;        // moisture transfer coefficient
};

MoninObukhovFluxes compute_monin_obukhov_fluxes(
    const SurfaceConfig& sfc_cfg,
    double u_sfc,           // surface zonal wind [m/s]
    double v_sfc,           // surface meridional wind [m/s]
    double theta_sfc,       // surface potential temperature [K]
    double qv_sfc,          // surface water vapor [kg/kg]
    double theta_air,       // air potential temperature [K]
    double qv_air,          // air water vapor [kg/kg]
    double p_air,           // air pressure [Pa]
    double rho_air,         // air density [kg/m³]
    double z_sfc            // height of measurement [m]
);

// Stability functions for Monin-Obukhov theory
double psi_m(double zeta);  // momentum stability function
double psi_h(double zeta);  // heat stability function

// Bulk Richardson number computation
double bulk_richardson_number(
    double theta_sfc,       // surface potential temperature [K]
    double theta_air,       // air potential temperature [K]
    double u_sfc,           // surface zonal wind [m/s]
    double v_sfc,           // surface meridional wind [m/s]
    double z_sfc            // height [m]
);

// Virtual potential temperature
double virtual_potential_temperature(double theta, double qv);

// Convert between theta and T tendencies (for coupling)
void convert_theta_tendency_to_temperature(
    const std::vector<double>& dthetadt,
    const std::vector<double>& theta,
    const std::vector<double>& p,
    std::vector<double>& dTdt
);
} // namespace surface_fluxes
