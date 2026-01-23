#include "surface_fluxes.hpp"
#include <cmath>
#include <algorithm>

/*This file contains the implementation of the surface fluxes module.
This module is used to compute the surface fluxes of the simulation.
Takes in the surface configuration, wind speed, virtual potential temperature, 
and air temperature and computes the surface fluxes.
Passes out the surface fluxes, momentum stresses, sensible heat flux, latent heat flux, and Obukhov length.
to the calling function for use in the simulation.*/


namespace surface_fluxes 
{

/*This function computes the bulk fluxes of the simulation. Takes in the surface configuration, wind speed, virtual potential temperature, 
and air temperature and computes the surface fluxes.
Passes out the surface fluxes, momentum stresses, sensible heat flux, latent heat flux, and Obukhov length.
to the calling function for use in the simulation.*/
BulkFluxes compute_bulk_fluxes(
    const SurfaceConfig& sfc_cfg,
    double u_sfc, double v_sfc,
    double theta_sfc, double qv_sfc,
    double theta_air, double qv_air,
    double p_air, double rho_air,
    double z_sfc
) 
{
    // Initialize the bulk fluxes structure
    BulkFluxes fluxes;

    // Wind speed
    double U = std::sqrt(u_sfc * u_sfc + v_sfc * v_sfc);
    U = std::max(U, 0.1);  // minimum wind speed to avoid division by zero

    // Virtual potential temperatures for stability
    double theta_v_sfc = virtual_potential_temperature(theta_sfc, qv_sfc);
    double theta_v_air = virtual_potential_temperature(theta_air, qv_air);

    // Bulk Richardson number for stability correction
    double Ri_b = bulk_richardson_number(theta_v_sfc, theta_v_air, u_sfc, v_sfc, z_sfc);

    // Stability-corrected transfer coefficients
    // Simple stability function (neutral: 1.0, stable: <1.0, unstable: >1.0)
    double stability_factor = 1.0;

    // If the bulk Richardson number is greater than 0, the stability is stable and the transfer coefficients are reduced.
    if (Ri_b > 0.0) 
    {
        // Stable: reduce transfer
        stability_factor = std::max(0.1, 1.0 - 5.0 * Ri_b);
    } 
    else 
    {
        // Unstable: enhance transfer
        stability_factor = 1.0 - 10.0 * Ri_b;
    }

    // Neutral drag coefficient (roughness-dependent)
    double Cd_neutral = std::pow(boundary_layer_constants::kappa /
                                 std::log(z_sfc / sfc_cfg.z0m), 2);

    // Apply stability correction
    fluxes.Cd = Cd_neutral * stability_factor;
    fluxes.Ch = 0.75 * Cd_neutral * stability_factor;  // typically Ch < Cd
    fluxes.Ce = fluxes.Ch;  // assume Ce ≈ Ch for neutral

    // Friction velocity
    fluxes.ustar = std::sqrt(fluxes.Cd) * U;
    fluxes.ustar = std::max(fluxes.ustar, 1e-4);  // minimum value

    // Momentum stresses
    fluxes.tau_u = -rho_air * fluxes.Cd * U * u_sfc;
    fluxes.tau_v = -rho_air * fluxes.Cd * U * v_sfc;

    // Sensible heat flux
    double delta_theta = theta_sfc - theta_air;
    fluxes.H = rho_air * boundary_layer_constants::cp * fluxes.Ch * U * delta_theta;

    // Latent heat flux (moisture flux)
    double delta_qv = qv_sfc - qv_air;
    fluxes.E = rho_air * fluxes.Ce * U * delta_qv;

    return fluxes;
}

/*This function computes the Monin-Obukhov fluxes of the simulation. Takes in the surface configuration, wind speed, virtual potential temperature, 
and air temperature and computes the surface fluxes.
Passes out the surface fluxes, momentum stresses, sensible heat flux, latent heat flux, and Obukhov length.
to the calling function for use in the simulation.*/
MoninObukhovFluxes compute_monin_obukhov_fluxes(
    const SurfaceConfig& sfc_cfg,
    double u_sfc, double v_sfc,
    double theta_sfc, double qv_sfc,
    double theta_air, double qv_air,
    double p_air, double rho_air,
    double z_sfc
) 
{
    // Initialize the Monin-Obukhov fluxes structure
    MoninObukhovFluxes fluxes;

    // Wind speed
    double U = std::sqrt(u_sfc * u_sfc + v_sfc * v_sfc);
    U = std::max(U, 0.1);

    // Virtual potential temperatures
    double theta_v_sfc = virtual_potential_temperature(theta_sfc, qv_sfc);
    double theta_v_air = virtual_potential_temperature(theta_air, qv_air);

    // Initial guess for Obukhov length (iterate)
    double L_prev = 1e6;  // large value for neutral
    double L = L_prev;
    const int max_iter = 10;

    // Iterate to compute the Obukhov length
    for (int iter = 0; iter < max_iter; ++iter) 
    {
        // Stability parameter
        fluxes.zeta = z_sfc / L;

        // Stability functions (WRF formulation) - psi_m and psi_h are the stability functions for momentum and heat, respectively.
        double psi_m_val = psi_m(fluxes.zeta);
        double psi_h_val = psi_h(fluxes.zeta);

        // Transfer coefficients with stability correction
        double ln_z_z0m = std::log(z_sfc / sfc_cfg.z0m);
        double ln_z_z0h = std::log(z_sfc / sfc_cfg.z0h);

        fluxes.Cd = std::pow(boundary_layer_constants::kappa / (ln_z_z0m - psi_m_val), 2);
        fluxes.Ch = boundary_layer_constants::kappa * boundary_layer_constants::kappa /
                   ((ln_z_z0m - psi_m_val) * (ln_z_z0h - psi_h_val));

        // Friction velocity
        fluxes.ustar = std::sqrt(fluxes.Cd) * U;

        // Heat flux (iterative)
        double delta_theta_v = theta_v_sfc - theta_v_air;
        fluxes.H = rho_air * boundary_layer_constants::cp * fluxes.Ch * U * delta_theta_v;

        // Moisture flux
        double delta_qv = qv_sfc - qv_air;
        fluxes.Ce = fluxes.Ch;  // assume same as heat
        fluxes.E = rho_air * fluxes.Ce * U * delta_qv;

        // Update Obukhov length
        double theta_star = -fluxes.H / (rho_air * boundary_layer_constants::cp * fluxes.ustar);
        double qv_star = -fluxes.E / (rho_air * fluxes.ustar);

        double buoyancy_flux = boundary_layer_constants::g / theta_v_air *
                              (theta_star + 0.61 * theta_air * qv_star);

        // If the buoyancy flux is greater than 1e-10, compute the Obukhov length.
        if (std::abs(buoyancy_flux) > 1e-10) 
        {
            L = -fluxes.ustar * fluxes.ustar * fluxes.ustar *
                theta_v_air / (boundary_layer_constants::kappa * boundary_layer_constants::g * buoyancy_flux);
        } 
        else 
        {
            L = 1e6;  // neutral
        }

        // If convergence is reached, break the loop.
        if (std::abs(L - L_prev) / std::abs(L_prev) < 1e-3) 
        {
            break;
        }
        L_prev = L;
    }

    fluxes.L = L;

    // Momentum stresses
    fluxes.tau_u = -rho_air * fluxes.Cd * U * u_sfc;
    fluxes.tau_v = -rho_air * fluxes.Cd * U * v_sfc;

    return fluxes;
}

/*This function computes the stability function for momentum. Done WRF-style. Takes in the stability parameter 
and passes out the stability function.*/
double psi_m(double zeta) 
{
    // If the stability is stable, use the stable function.
    if (zeta >= 0.0) 
    {
        // Stable
        double x = (1.0 + 5.0 * zeta);
        return -5.0 * zeta - 5.0 * zeta * zeta / x;
    } 
    else 
    {
        // Unstable
        double x = std::pow(1.0 - 16.0 * zeta, 0.25);
        return 2.0 * std::log((1.0 + x) / 2.0) + std::log((1.0 + x * x) / 2.0) -
               2.0 * std::atan(x) + boundary_layer_constants::pi / 2.0;
    }
}

/*This function computes the stability function for heat. Done WRF-style. Takes in the stability parameter 
and passes out the stability function.*/
double psi_h(double zeta) 
{
    // If the stability is stable, use the stable function.
    if (zeta >= 0.0) 
    {
        // Stable
        return -5.0 * zeta;
    } 
    else 
    {
        // Unstable
        double x = std::pow(1.0 - 16.0 * zeta, 0.25);
        return 2.0 * std::log((1.0 + x * x) / 2.0);
    }
}

/*This function computes the bulk Richardson number. Takes in the virtual potential temperature at the surface, 
the virtual potential temperature of the air, the wind speed, and the height of the surface.
Passes out the bulk Richardson number.*/
double bulk_richardson_number(double theta_v_sfc, double theta_v_air,double u_sfc, double v_sfc, double z_sfc) 
{
    double delta_theta_v = theta_v_sfc - theta_v_air;
    double U = std::sqrt(u_sfc * u_sfc + v_sfc * v_sfc);
    U = std::max(U, 0.1);

    return boundary_layer_constants::g * delta_theta_v * z_sfc /
           (theta_v_air * U * U);
}

/*This function computes the virtual potential temperature. Takes in the potential temperature and 
the water vapor mixing ratioand passes out the virtual potential temperature.*/
double virtual_potential_temperature(double theta, double qv) 
{
    return theta * (1.0 + 0.61 * qv);
}

/*This function converts the potential temperature tendency to the temperature tendency. Takes in the potential temperature tendency, 
the potential temperature, the pressure, and passes out the temperature tendency.*/
void convert_theta_tendency_to_temperature(
    const std::vector<double>& dthetadt,
    const std::vector<double>& theta,
    const std::vector<double>& p,
    std::vector<double>& dTdt
) 
{
    // Initialize the temperature tendency structure
    const size_t nz = dthetadt.size();
    dTdt.resize(nz);

    const double R_d = 287.0;
    const double cp = boundary_layer_constants::cp;
    const double kappa = R_d / cp;
    const double p0 = 100000.0;

    // Convert the potential temperature tendency to the tem    perature tendency
    for (size_t k = 0; k < nz; ++k) 
    {
        // dT/dt = (p/p0)^(-kappa) * dtheta/dt
        dTdt[k] = dthetadt[k] * std::pow(p[k] / p0, -kappa);
    }
}

} // namespace surface_fluxes
