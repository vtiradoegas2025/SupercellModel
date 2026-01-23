#pragma once
#include <vector>
#include <string>
#include <memory>

/**
 * @brief Atmospheric sounding data structure
 * Contains vertical profiles of thermodynamic and kinematic variables
 */
struct SoundingData 
{
    std::vector<double> height_m;           // Height above ground (m)
    std::vector<double> pressure_hpa;       // Pressure (hPa)
    std::vector<double> temperature_k;      // Temperature (K)
    std::vector<double> dewpoint_k;         // Dewpoint temperature (K)
    std::vector<double> wind_speed_ms;      // Wind speed (m/s)
    std::vector<double> wind_direction_deg; // Wind direction (degrees from north)
    std::vector<double> relative_humidity;  // Relative humidity (%)
    std::vector<double> mixing_ratio_kgkg;  // Water vapor mixing ratio (kg/kg)

    // Derived quantities
    std::vector<double> potential_temperature_k; // Potential temperature (K)
    std::vector<double> equivalent_potential_temperature_k; // Equivalent potential temperature (K)

    // Metadata
    std::string station_id;
    std::string timestamp_utc;
    double latitude_deg;
    double longitude_deg;
    double elevation_m;

    // Quality flags
    std::vector<bool> valid_temperature;
    std::vector<bool> valid_moisture;
    std::vector<bool> valid_wind;

    /**
     * @brief Check if sounding data is valid
     * @return true if basic data is present
     */
    bool is_valid() const 
    {
        return !height_m.empty() &&
               height_m.size() == temperature_k.size() &&
               height_m.size() == pressure_hpa.size();
    }

    /**
     * @brief Get number of valid levels
     * @return Number of data points
     */
    size_t num_levels() const 
    {
        return height_m.size();
    }

    /**
     * @brief Clear all data
     */
    void clear() 
    {
        height_m.clear();
        pressure_hpa.clear();
        temperature_k.clear();
        dewpoint_k.clear();
        wind_speed_ms.clear();
        wind_direction_deg.clear();
        relative_humidity.clear();
        mixing_ratio_kgkg.clear();
        potential_temperature_k.clear();
        equivalent_potential_temperature_k.clear();
        valid_temperature.clear();
        valid_moisture.clear();
        valid_wind.clear();
    }
};

/**
 * @brief Configuration for sounding schemes
 */
struct SoundingConfig {
    std::string scheme_id = "none";           // Sounding scheme identifier
    std::string file_path = "";               // Path to sounding file
    bool use_fallback_profiles = true;        // Use procedural profiles if sounding fails
    double interpolation_method = 0;          // 0=linear, 1=spline, 2=log-linear for pressure
    bool extrapolate_below_ground = false;    // Extrapolate data below ground level
    bool extrapolate_above_top = false;       // Extrapolate data above sounding top

    // Quality control thresholds
    double min_pressure_hpa = 100.0;          // Minimum valid pressure
    double max_pressure_hpa = 1100.0;         // Maximum valid pressure
    double min_temperature_k = 200.0;         // Minimum valid temperature
    double max_temperature_k = 350.0;         // Maximum valid temperature
};

/**
 * @brief Base class for atmospheric sounding schemes
 * Provides interface for reading and processing sounding data
 */
class SoundingScheme 
{
public:
    virtual ~SoundingScheme() = default;

    /**
     * @brief Initialize the sounding scheme
     * @param config Sounding configuration
     */
    virtual void initialize(const SoundingConfig& config) = 0;

    /**
     * @brief Load sounding data from file
     * @param file_path Path to sounding file
     * @return SoundingData structure with loaded data
     */
    virtual SoundingData load_sounding(const std::string& file_path) = 0;

    /**
     * @brief Interpolate sounding data to target heights
     * @param sounding Original sounding data
     * @param target_heights_m Target heights for interpolation
     * @return Interpolated sounding data
     */
    virtual SoundingData interpolate_to_heights(
        const SoundingData& sounding,
        const std::vector<double>& target_heights_m) = 0;

    /**
     * @brief Get current configuration
     * @return Reference to current configuration
     */
    virtual const SoundingConfig& get_config() const = 0;

    /**
     * @brief Check if scheme is properly initialized
     * @return true if ready to use
     */
    virtual bool is_initialized() const = 0;
};

/**
 * @brief Factory function for creating sounding schemes
 * @param scheme_id Identifier for the desired scheme
 * @return Unique pointer to sounding scheme instance
 */
std::unique_ptr<SoundingScheme> create_sounding_scheme(const std::string& scheme_id);
