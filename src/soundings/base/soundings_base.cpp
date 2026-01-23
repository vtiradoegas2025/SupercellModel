#include "../../include/soundings_base.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <stdexcept>

/*This file contains the implementation of the soundings base class.*/

namespace 
{

// Physical constants for thermodynamic calculations
const double R_d = 287.0;     // Dry air gas constant (J/kg·K)
const double cp = 1004.0;     // Specific heat at constant pressure (J/kg·K)
const double L_v = 2.5e6;     // Latent heat of vaporization (J/kg)
const double p0 = 100000.0;   // Reference pressure (Pa)
const double T0 = 273.15;     // Freezing temperature (K)

/**
 * @brief Calculate potential temperature
 * @param temperature_k Temperature in Kelvin
 * @param pressure_hpa Pressure in hPa
 * @return Potential temperature in Kelvin
 */
double calculate_potential_temperature(double temperature_k, double pressure_hpa) {
    const double kappa = R_d / cp;
    const double p_pa = pressure_hpa * 100.0;  // Convert to Pa
    return temperature_k * std::pow(p0 / p_pa, kappa);
}

/**
 * @brief Calculate mixing ratio from dewpoint
 * @param dewpoint_k Dewpoint temperature in Kelvin
 * @param pressure_hpa Pressure in hPa
 * @return Mixing ratio in kg/kg
 */
double calculate_mixing_ratio(double dewpoint_k, double pressure_hpa) {
    // Tetens formula for saturation vapor pressure
    const double e_s = 6.112 * std::exp(17.67 * (dewpoint_k - T0) / (dewpoint_k - 29.65));
    const double p_pa = pressure_hpa * 100.0;  // Convert to Pa
    return (0.62198 * e_s) / (p_pa/100.0 - e_s);  // kg/kg
}

/**
 * @brief Calculate equivalent potential temperature
 * @param temperature_k Temperature in Kelvin
 * @param dewpoint_k Dewpoint temperature in Kelvin
 * @param pressure_hpa Pressure in hPa
 * @return Equivalent potential temperature in Kelvin
 */
double calculate_equivalent_potential_temperature(double temperature_k, double dewpoint_k, double pressure_hpa) {
    const double r_v = calculate_mixing_ratio(dewpoint_k, pressure_hpa);
    const double T_L = temperature_k * std::pow(p0 / (pressure_hpa * 100.0), R_d/(cp + L_v * r_v / temperature_k));
    return T_L * std::exp((L_v * r_v) / (cp * T_L));
}

} // anonymous namespace

/**
 * @brief Linear interpolation between two points
 * @param x1 First x value
 * @param y1 First y value
 * @param x2 Second x value
 * @param y2 Second y value
 * @param x_target Target x value
 * @return Interpolated y value
 */
double linear_interpolate(double x1, double y1, double x2, double y2, double x_target) {
    if (x1 == x2) return y1;
    return y1 + (y2 - y1) * (x_target - x1) / (x2 - x1);
}

/**
 * @brief Quality control for sounding data
 * @param data Sounding data to validate
 * @param config Configuration with QC thresholds
 * @return true if data passes QC
 */
bool quality_control_sounding(SoundingData& data, const SoundingConfig& config) 
{
    // If the sounding data is invalid, return false.
    if (!data.is_valid()) 
    {
        std::cerr << "Sounding data is invalid - missing basic fields" << std::endl;
        return false;
    }

    size_t original_size = data.num_levels();
    std::vector<size_t> valid_indices;

    // Iterate over the vertical levels and perform the quality control.
    for (size_t i = 0; i < data.num_levels(); ++i) 
    {
        bool valid_level = true;

        // If the pressure is outside the range, set the valid level to false.
        if (data.pressure_hpa[i] < config.min_pressure_hpa ||
            data.pressure_hpa[i] > config.max_pressure_hpa) 
        {
            valid_level = false;
        }

        // If the temperature is outside the range, set the valid level to false.
        if (data.temperature_k[i] < config.min_temperature_k ||
            data.temperature_k[i] > config.max_temperature_k) 
        {
            valid_level = false;
        }

        // If the height is not monotonic, set the valid level to false.
        if (i > 0 && data.height_m[i] <= data.height_m[i-1]) 
        {
            valid_level = false;
        }

        // If the level is valid, add the index to the valid indices.
        if (valid_level) 
        {
            valid_indices.push_back(i);
        }
    }

    // If we have too few valid levels, reject the sounding
    if (valid_indices.size() < 5) 
    {
        std::cerr << "Insufficient valid levels in sounding: " << valid_indices.size()
                  << "/" << original_size << std::endl;
        return false;
    }

    // Filter data to only valid levels
    SoundingData filtered_data = data;
    filtered_data.clear();

    // Iterate over the valid indices and filter the data.
    for (size_t idx : valid_indices) 
    {
        filtered_data.height_m.push_back(data.height_m[idx]);
        filtered_data.pressure_hpa.push_back(data.pressure_hpa[idx]);
        filtered_data.temperature_k.push_back(data.temperature_k[idx]);

        // If the dewpoint is available, add the dewpoint to the filtered data.
        if (!data.dewpoint_k.empty()) 
        {
            filtered_data.dewpoint_k.push_back(data.dewpoint_k[idx]);
        }

        // If the wind speed is available, add the wind speed to the filtered data.
        if (!data.wind_speed_ms.empty()) 
        {
            filtered_data.wind_speed_ms.push_back(data.wind_speed_ms[idx]);
        }

        // If the wind direction is available, add the wind direction to the filtered data.
        if (!data.wind_direction_deg.empty()) 
        {
            filtered_data.wind_direction_deg.push_back(data.wind_direction_deg[idx]);
        }
    }

    data = filtered_data;

    std::cout << "Sounding quality control: " << data.num_levels()
              << " valid levels retained from " << original_size << std::endl;

    return true;
}

/**
 * @brief Calculate derived thermodynamic quantities
 * @param data Sounding data with basic fields
 */
void calculate_derived_quantities(SoundingData& data) 
{
    data.potential_temperature_k.resize(data.num_levels());
    data.equivalent_potential_temperature_k.resize(data.num_levels());
    data.mixing_ratio_kgkg.resize(data.num_levels());

    // Iterate over the vertical levels and calculate the derived thermodynamic quantities.
    for (size_t i = 0; i < data.num_levels(); ++i) 
    {
        // Potential temperature
        data.potential_temperature_k[i] = calculate_potential_temperature(
            data.temperature_k[i], data.pressure_hpa[i]);

        // If the dewpoint is available, calculate the mixing ratio and equivalent potential temperature.
        if (!data.dewpoint_k.empty() && i < data.dewpoint_k.size()) 
        {
            data.mixing_ratio_kgkg[i] = calculate_mixing_ratio(
                data.dewpoint_k[i], data.pressure_hpa[i]);

            // Equivalent potential temperature
            data.equivalent_potential_temperature_k[i] =
                calculate_equivalent_potential_temperature(
                    data.temperature_k[i], data.dewpoint_k[i], data.pressure_hpa[i]);
        }
    }
}

/**
 * @brief Interpolate sounding data to target heights using linear interpolation
 * @param source Source sounding data
 * @param target_heights Target heights in meters
 * @param config Sounding configuration
 * @return Interpolated sounding data
 */
SoundingData interpolate_sounding_linear(
    const SoundingData& source,
    const std::vector<double>& target_heights,
    const SoundingConfig& config) {

    SoundingData result;
    result.height_m = target_heights;

    // Initialize result arrays
    result.pressure_hpa.resize(target_heights.size());
    result.temperature_k.resize(target_heights.size());
    result.potential_temperature_k.resize(target_heights.size());

    // If the dewpoint is available, resize the result arrays.
    if (!source.dewpoint_k.empty()) 
    {
        result.dewpoint_k.resize(target_heights.size());
        result.mixing_ratio_kgkg.resize(target_heights.size());
        result.equivalent_potential_temperature_k.resize(target_heights.size());
    }

    // If the wind speed is available, resize the result arrays.
    if (!source.wind_speed_ms.empty()) {
        result.wind_speed_ms.resize(target_heights.size());
        result.wind_direction_deg.resize(target_heights.size());
    }

    // Iterate over the target heights and interpolate the sounding data.
    for (size_t i = 0; i < target_heights.size(); ++i) 
    {
        double target_h = target_heights[i];

        // Find bracketing levels in source data
        auto it = std::lower_bound(source.height_m.begin(), source.height_m.end(), target_h);

        // If the target height is below the lowest level, extrapolate using the bottom two levels.
        if (it == source.height_m.begin()) 
        {
            // If the extrapolation below the ground is requested, extrapolate using the bottom two levels.
            if (config.extrapolate_below_ground) 
            {
                // Extrapolate using bottom two levels
                size_t idx1 = 0, idx2 = 1;

                // If there are more than one level, interpolate using the bottom two levels.
                if (source.height_m.size() > 1)
                {
                    result.pressure_hpa[i] = linear_interpolate(
                        source.height_m[idx1], source.pressure_hpa[idx1],
                        source.height_m[idx2], source.pressure_hpa[idx2], target_h);
                    result.temperature_k[i] = linear_interpolate(
                        source.height_m[idx1], source.temperature_k[idx1],
                        source.height_m[idx2], source.temperature_k[idx2], target_h);
                } 
                else 
                {
                    // Single level - use as is
                    result.pressure_hpa[i] = source.pressure_hpa[0];
                    result.temperature_k[i] = source.temperature_k[0];
                }
            } 
            else 
            {
                // Use lowest level
                result.pressure_hpa[i] = source.pressure_hpa[0];
                result.temperature_k[i] = source.temperature_k[0];
            }
        } 
        // If the target height is above the highest level, extrapolate using the top two levels.
        else if (it == source.height_m.end()) 
        {
            // If the extrapolation above the top is requested, extrapolate using the top two levels.
            if (config.extrapolate_above_top) 
            {
                // Extrapolate using top two levels
                size_t idx1 = source.height_m.size() - 2;
                size_t idx2 = source.height_m.size() - 1;

                // If there are more than one level, interpolate using the top two levels.
                if (source.height_m.size() > 1) 
                {
                    result.pressure_hpa[i] = linear_interpolate(
                        source.height_m[idx1], source.pressure_hpa[idx1],
                        source.height_m[idx2], source.pressure_hpa[idx2], target_h);
                    result.temperature_k[i] = linear_interpolate(
                        source.height_m[idx1], source.temperature_k[idx1],
                        source.height_m[idx2], source.temperature_k[idx2], target_h);
                } 
                else 
                {
                    result.pressure_hpa[i] = source.pressure_hpa.back();
                    result.temperature_k[i] = source.temperature_k.back();
                }
            } 
            else 
            {
                // Use highest level
                result.pressure_hpa[i] = source.pressure_hpa.back();
                result.temperature_k[i] = source.temperature_k.back();
            }
        } 
        else 
        {
            // Interpolate between bracketing levels
            size_t idx2 = it - source.height_m.begin();
            size_t idx1 = idx2 - 1;

            result.pressure_hpa[i] = linear_interpolate(
                source.height_m[idx1], source.pressure_hpa[idx1],
                source.height_m[idx2], source.pressure_hpa[idx2], target_h);
            result.temperature_k[i] = linear_interpolate(
                source.height_m[idx1], source.temperature_k[idx1],
                source.height_m[idx2], source.temperature_k[idx2], target_h);

            // If the dewpoint is available, interpolate the dewpoint.
            if (!source.dewpoint_k.empty()) 
            {
                result.dewpoint_k[i] = linear_interpolate(
                    source.height_m[idx1], source.dewpoint_k[idx1],
                    source.height_m[idx2], source.dewpoint_k[idx2], target_h);
            }

            // If the wind speed is available, interpolate the wind speed.
            if (!source.wind_speed_ms.empty()) 
            {
                result.wind_speed_ms[i] = linear_interpolate(
                    source.height_m[idx1], source.wind_speed_ms[idx1],
                    source.height_m[idx2], source.wind_speed_ms[idx2], target_h);
            }

            // If the wind direction is available, interpolate the wind direction.
            if (!source.wind_direction_deg.empty()) 
            {
                result.wind_direction_deg[i] = linear_interpolate(
                    source.height_m[idx1], source.wind_direction_deg[idx1],
                    source.height_m[idx2], source.wind_direction_deg[idx2], target_h);
            }
        }
    }

    // Calculate derived quantities
    calculate_derived_quantities(result);

    return result;
}
