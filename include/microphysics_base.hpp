#pragma once
#include <vector>
#include <memory>
#include <string>
#include "field3d.hpp"

/*This header file contains the base classes and structures for the microphysics module.
The microphysics module is responsible for the microphysics of the simulation.
The microphysics scheme is chosen by the user in the configuration file.
This module is used to compute the microphysics of the simulation.*/

class MicrophysicsScheme 
{
public:
    virtual ~MicrophysicsScheme() = default;

    // Core interface - compute tendencies for all hydrometeor categories
    virtual void compute_tendencies(
        // Input state
        const Field3D& p,        // pressure
        const Field3D& theta,    // potential temperature
        const Field3D& qv,       // water vapor mixing ratio
        const Field3D& qc,       // cloud water mixing ratio
        const Field3D& qr,       // rain mixing ratio
        const Field3D& qi,       // cloud ice mixing ratio
        const Field3D& qs,       // snow mixing ratio
        const Field3D& qg,       // graupel mixing ratio
        const Field3D& qh,       // hail mixing ratio
        double dt,                                                    // timestep
        // Output tendencies (will be added to existing fields)
        Field3D& dtheta_dt,      // potential temperature tendency
        Field3D& dqv_dt,         // water vapor tendency
        Field3D& dqc_dt,         // cloud water tendency
        Field3D& dqr_dt,         // rain tendency
        Field3D& dqi_dt,         // cloud ice tendency
        Field3D& dqs_dt,         // snow tendency
        Field3D& dqg_dt,         // graupel tendency
        Field3D& dqh_dt          // hail tendency
    ) = 0;

    // Optional: radar reflectivity calculation
    virtual void compute_radar_reflectivity(
        const Field3D& qc,
        const Field3D& qr,
        const Field3D& qi,
        const Field3D& qs,
        const Field3D& qg,
        const Field3D& qh,
        Field3D& reflectivity_dbz
    ) {}

    // Optional: precipitation rates at surface
    virtual void compute_precipitation_rates(
        const Field3D& qr,
        const Field3D& qs,
        const Field3D& qg,
        const Field3D& qh,
        Field3D& precip_rate_rain,  // mm/hr
        Field3D& precip_rate_snow,  // mm/hr
        Field3D& precip_rate_grau,  // mm/hr
        Field3D& precip_rate_hail   // mm/hr
    ) {}

    // Get scheme name for diagnostics
    virtual std::string get_scheme_name() const = 0;

    // Get number of prognostic variables (for metadata)
    virtual int get_num_prognostic_vars() const = 0;
};

// Factory function declaration
std::unique_ptr<MicrophysicsScheme> create_microphysics_scheme(const std::string& scheme_name);

// Common physical constants used by microphysics schemes
namespace microphysics_constants 
{
    constexpr double T0 = 273.15;      // freezing temperature (K)
    constexpr double L_v = 2.5e6;      // latent heat of vaporization (J/kg)
    constexpr double L_f = 3.34e5;     // latent heat of fusion (J/kg)
    constexpr double L_s = L_v + L_f;  // latent heat of sublimation (J/kg)
    constexpr double cp = 1004.0;      // specific heat at constant pressure (J/kg·K)
    constexpr double R_v = 461.5;      // gas constant for water vapor (J/kg·K)
    constexpr double R_d = 287.0;      // gas constant for dry air (J/kg·K)
    constexpr double p0 = 100000.0;    // reference pressure (Pa)
    constexpr double rho_w = 1000.0;   // water density (kg/m³)
    constexpr double rho_i = 917.0;    // ice density (kg/m³)
}
