#pragma once
#include <vector>
#include <cmath>
#include "../../../include/microphysics_base.hpp"

/*This namespace contains the thermodynamics utilities for the microphysics schemes.
This namespace contains the theta_to_temperature, temperature_to_theta, 
saturation_vapor_pressure_water, saturation_vapor_pressure_ice, saturation_mixing_ratio_water, 
saturation_mixing_ratio_ice, relative_humidity, latent_heat_vaporization, virtual_temperature, 
air_density, temperature_tendency_to_theta, and saturation_adjustment 
functions.*/
namespace thermodynamics 
{

/*This function converts the potential temperature to the temperature.
Takes in the potential temperature and the pressure
and converts the potential temperature to the temperature.*/
inline double theta_to_temperature(double theta, double p) 
{
    using namespace microphysics_constants;
    double kappa = R_d / cp;
    return theta * pow(p / p0, kappa);
}


/*This function converts the temperature to the potential temperature.
Takes in the temperature and the pressure
and converts the temperature to the potential temperature.*/
inline double temperature_to_theta(double T, double p) 
{
    using namespace microphysics_constants;
    double kappa = R_d / cp;
    return T * pow(p0 / p, kappa);
}

/*This function computes the saturation vapor pressure over water.
Takes in the temperature and computes the saturation vapor pressure
over water.*/
inline double saturation_vapor_pressure_water(double T) 
{
    using namespace microphysics_constants;
    double T_c = T - T0;  // Celsius

    // If the temperature is greater than or equal to the freezing temperature, compute the saturation vapor pressure over water.
    if (T >= T0) 
    {
        // Over water
        return 611.21 * exp((18.678 - T_c/234.5) * T_c / (257.14 + T_c));
    } 
    else 
    {
        // Over ice (use ice formula below freezing)
        return 611.15 * exp((23.036 - T_c/333.7) * T_c / (279.82 + T_c));
    }
}

/*This function computes the saturation vapor pressure over ice.
Takes in the temperature and computes the saturation vapor pressure over ice.*/
inline double saturation_vapor_pressure_ice(double T) 
{
    using namespace microphysics_constants;
    double T_c = T - T0;  // Celsius
    return 611.15 * exp((23.036 - T_c/333.7) * T_c / (279.82 + T_c));
}

/*This function computes the saturation mixing ratio over water.
Takes in the temperature and the pressure
and computes the saturation mixing ratio over water.*/
inline double saturation_mixing_ratio_water(double T, double p) 
{
    using namespace microphysics_constants;
    double e_sat = saturation_vapor_pressure_water(T);
    return 0.62198 * e_sat / (p - e_sat);
}


/*This function computes the saturation mixing ratio over ice.
Takes in the temperature and the pressure
and computes the saturation mixing ratio over ice.*/
inline double saturation_mixing_ratio_ice(double T, double p) 
{
    using namespace microphysics_constants;
    double e_sat = saturation_vapor_pressure_ice(T);
    return 0.62198 * e_sat / (p - e_sat);
}

/*This function computes the relative humidity.
Takes in the vapor mixing ratio and the saturation vapor mixing ratio
and computes the relative humidity.*/
inline double relative_humidity(double qv, double qvs) 
{
    return (qvs > 0.0) ? qv / qvs : 0.0;
}

/*This function computes the latent heat of vaporization.
Takes in the temperature and computes the latent heat of vaporization.*/
inline double latent_heat_vaporization(double T) 
{
    using namespace microphysics_constants;
    // Simple temperature dependence
    return L_v - 2370.0 * (T - T0);
}

/*This function computes the virtual temperature.
Takes in the temperature, the vapor mixing ratio, the cloud water mixing ratio, and the rainwater mixing ratio
and computes the virtual temperature.*/
inline double virtual_temperature(double T, double qv, double qc = 0.0, double qr = 0.0) 
{
    using namespace microphysics_constants;
    double q_total = qv + qc + qr;  // approximate total condensate
    return T * (1.0 + 0.608 * qv - q_total);
}

/*This function computes the air density.
Takes in the pressure, the temperature, and the vapor mixing ratio
and computes the air density.*/
inline double air_density(double p, double T, double qv = 0.0) 
{
    using namespace microphysics_constants;
    double Tv = virtual_temperature(T, qv);
    return p / (R_d * Tv);
}

/*This function computes the temperature tendency to the potential temperature tendency.
Takes in the temperature tendency, the potential temperature, and the pressure
and computes the temperature tendency to the potential temperature tendency.*/
inline double temperature_tendency_to_theta(double dT_dt, double theta, double p) 
{
    using namespace microphysics_constants;
    double kappa = R_d / cp;
    return dT_dt * pow(p0 / p, kappa);
}

/*This function computes the saturation adjustment.
Takes in the temperature, the pressure, the vapor mixing ratio, 
and the cloud water mixing ratioand computes the saturation adjustment.*/
inline double saturation_adjustment(double T_in, double p, double& qv, double& qc) 
{
    using namespace microphysics_constants;

    double T = T_in;
    double qvs = saturation_mixing_ratio_water(T, p);

    // If the vapor mixing ratio is greater than the saturation vapor mixing ratio, compute the saturation adjustment.
    if (qv > qvs) 
    {
        // Excess vapor becomes cloud water
        double excess = qv - qvs;
        qc += excess;
        qv = qvs;
    } 
    // If the vapor mixing ratio is less than the saturation vapor mixing ratio 
    // and the cloud water mixing ratio is greater than 0, compute the saturation adjustment.
    else if (qv < qvs && qc > 0.0) 
    {
        // Evaporate cloud water to reach saturation
        double deficit = qvs - qv;
        double evaporate = std::min(deficit, qc);
        qc -= evaporate;
        qv += evaporate;
    }

    // Ensure non-negative values
    qc = std::max(0.0, qc);
    qv = std::max(0.0, qv);

    return T;
}

/*This function converts the potential temperature field to the temperature field.
Takes in the potential temperature field, the pressure field, and the temperature field
and converts the potential temperature field to the temperature field.*/
void convert_theta_to_temperature_field(
    const std::vector<std::vector<std::vector<float>>>& theta,
    const std::vector<std::vector<std::vector<float>>>& p,
    std::vector<std::vector<std::vector<float>>>& temperature
);

/*This function converts the temperature field to the potential temperature field.
Takes in the temperature field, the pressure field, and the potential temperature field
and converts the temperature field to the potential temperature field.*/
void convert_temperature_to_theta_field(
    const std::vector<std::vector<std::vector<float>>>& temperature,
    const std::vector<std::vector<std::vector<float>>>& p,
    std::vector<std::vector<std::vector<float>>>& theta
);

} // namespace thermodynamics
