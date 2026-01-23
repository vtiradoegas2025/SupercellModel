#include "sharpy_sounding.hpp"
#include "../../base/soundings_base.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <filesystem>

// For now, we'll use standard C++ file I/O and create a framework
// that can be extended with HDF5/NetCDF libraries later
// #include <hdf5.h>      // Uncomment when HDF5 support is added
// #include <netcdf.h>    // Uncomment when NetCDF support is added

SharpySoundingScheme::SharpySoundingScheme()
    : initialized_(false) {
}

/*This function initializes the SHARPY sounding scheme.
Takes in the configuration and initializes the SHARPY sounding scheme.*/
void SharpySoundingScheme::initialize(const SoundingConfig& config) 
{
    config_ = config;
    initialized_ = true;

    std::cout << "Initialized SHARPY sounding scheme:" << std::endl;
    std::cout << "  File path: " << config_.file_path << std::endl;
    std::cout << "  Interpolation method: " <<
        (config_.interpolation_method == 0 ? "linear" :
         config_.interpolation_method == 1 ? "spline" : "log-linear") << std::endl;
    std::cout << "  Extrapolate below ground: " << (config_.extrapolate_below_ground ? "yes" : "no") << std::endl;
    std::cout << "  Extrapolate above top: " << (config_.extrapolate_above_top ? "yes" : "no") << std::endl;

    // If the file path is not empty and the file does not exist, throw an error.
    if (!config_.file_path.empty() && !std::filesystem::exists(config_.file_path)) 
    {
        std::cerr << "Warning: SHARPY file does not exist: " << config_.file_path << std::endl;

        // If the fallback profiles are not used, throw an error.
        if (!config_.use_fallback_profiles) 
        {
            throw std::runtime_error("SHARPY file not found and fallback profiles disabled");
        }
    }
}

/*This function loads the SHARPY sounding.
Takes in the file path and loads the SHARPY sounding.*/
SoundingData SharpySoundingScheme::load_sounding(const std::string& file_path) 
{
    // If the scheme is not initialized, throw an error.
    if (!initialized_) 
    {
        throw std::runtime_error("SHARPY sounding scheme not initialized");
    }

    std::cout << "Loading SHARPY sounding from: " << file_path << std::endl;

    // If the file does not exist, throw an error.
    if (!std::filesystem::exists(file_path)) 
    {
        throw std::runtime_error("SHARPY file not found: " + file_path);
    }

    // Detect file format
    std::string format = detect_file_format(file_path);
    std::cout << "Detected file format: " << format << std::endl;

    SoundingData data;

    try 
    {
        // If the file format is HDF5, read the SHARPY sounding.
        if (format == "hdf5") 
        {
            data = read_sharpy_hdf5(file_path);
        }

        // If the file format is NetCDF, read the SHARPY sounding.
         else if (format == "netcdf") 
        {
            data = read_sharpy_netcdf(file_path);
        } 
        else 
        {
            throw std::runtime_error("Unsupported file format: " + format);
        }

        // If the sounding fails quality control, throw an error.
        if (!quality_control_sounding(data, config_)) 
        {
            throw std::runtime_error("Sounding failed quality control");
        }

        // Calculate derived quantities
        calculate_derived_quantities(data);

        std::cout << "Successfully loaded sounding with " << data.num_levels() << " levels" << std::endl;
        std::cout << "  Height range: " << data.height_m.front() << " - " << data.height_m.back() << " m" << std::endl;
        std::cout << "  Pressure range: " << data.pressure_hpa.front() << " - " << data.pressure_hpa.back() << " hPa" << std::endl;

    }
     catch (const std::exception& e) 
     {
        std::cerr << "Error loading SHARPY sounding: " << e.what() << std::endl;

        // If the fallback profiles are not used, throw an error.
        if (!config_.use_fallback_profiles) 
        {
            throw;
        }
        std::cerr << "Using fallback profiles as configured" << std::endl;

        // Return empty data to signal fallback should be used
        data.clear();
    }

    return data;
}

/*This function interpolates the SHARPY sounding to the target heights.
Takes in the sounding data and the target heights and interpolates the SHARPY sounding to the target heights.*/
SoundingData SharpySoundingScheme::interpolate_to_heights(
    const SoundingData& sounding,
    const std::vector<double>& target_heights_m) 
    {

    // If the scheme is not initialized, throw an error.
    if (!initialized_) 
    {
        throw std::runtime_error("SHARPY sounding scheme not initialized");
    }

    // If the sounding data is invalid, throw an error.
    if (!sounding.is_valid()) 
    {
        throw std::runtime_error("Invalid sounding data provided for interpolation");
    }

    std::cout << "Interpolating sounding to " << target_heights_m.size() << " target heights" << std::endl;

    // Switch block to use the appropriate interpolation method.
    SoundingData interpolated;
    switch (static_cast<int>(config_.interpolation_method)) 
    {
        case 0: // Linear
        default:
            interpolated = interpolate_sounding_linear(sounding, target_heights_m, config_);
            break;
        case 1: // Spline (placeholder - would need spline library)
            std::cout << "Warning: Spline interpolation not yet implemented, using linear" << std::endl;
            interpolated = interpolate_sounding_linear(sounding, target_heights_m, config_);
            break;
        case 2: // Log-linear for pressure
            // For now, use linear - log-linear would need special handling
            std::cout << "Warning: Log-linear interpolation not yet implemented, using linear" << std::endl;
            interpolated = interpolate_sounding_linear(sounding, target_heights_m, config_);
            break;
    }

    return interpolated;
}

/*This function detects the file format of the SHARPY sounding.
Takes in the file path and detects the file format of the SHARPY sounding.*/
std::string SharpySoundingScheme::detect_file_format(const std::string& file_path) 
{
    // Check file extension first
    std::string extension = std::filesystem::path(file_path).extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

    // If the extension is HDF5, return "hdf5".
    if (extension == ".h5" || extension == ".hdf5") 
    {
        return "hdf5";
    } else if (extension == ".nc" || extension == ".nc4" || extension == ".netcdf") {
        return "netcdf";
    }

    // Try to detect from file content (basic check)
    std::ifstream file(file_path, std::ios::binary);

    // If the file is not open, return "unknown".
    if (!file.is_open()) 
    {
        return "unknown";
    }

    // Read first few bytes to check file signature
    char buffer[8];
    file.read(buffer, 8);

    // If the file has at least 8 bytes, check the file signature.
    if (file.gcount() >= 8) 
    {
        // HDF5 signature: \211HDF\r\n\032\n
        if (buffer[0] == '\211' && buffer[1] == 'H' && buffer[2] == 'D' &&
            buffer[3] == 'F' && buffer[4] == '\r' && buffer[5] == '\n' &&
            buffer[6] == '\032' && buffer[7] == '\n') {
            return "hdf5";
        }
        // NetCDF signature (CDF): CDF001, CDF002, etc.
        if (buffer[0] == 'C' && buffer[1] == 'D' && buffer[2] == 'F' &&
            buffer[3] >= '0' && buffer[3] <= '9') {
            return "netcdf";
        }
    }

    return "unknown";
}

/*This function validates the SHARPY file.
Takes in the file path and validates the SHARPY file.*/
bool SharpySoundingScheme::validate_sharpy_file(const std::string& file_path) 
{
    // Basic validation - check if file exists and has expected structure
    if (!std::filesystem::exists(file_path)) 
    {
        return false;
    }

    // For now, just check file size is reasonable (> 1KB)
    auto file_size = std::filesystem::file_size(file_path);
    if (file_size < 1024) 
    {
        std::cerr << "Warning: SHARPY file seems too small: " << file_size << " bytes" << std::endl;
        return false;
    }

    return true;
}

/*This function reads the SHARPY HDF5 file.
Takes in the file path and reads the SHARPY HDF5 file.*/
SoundingData SharpySoundingScheme::read_sharpy_hdf5(const std::string& file_path) 
{
    // Placeholder implementation - would need HDF5 library
    // This shows the structure that would be implemented

    std::cout << "Reading SHARPY HDF5 file (placeholder implementation)" << std::endl;
    std::cout << "To implement: add HDF5 library dependency and read SHARPY data structure" << std::endl;
    std::cout << "Expected SHARPY HDF5 structure:" << std::endl;
    std::cout << "  /profiles/height_m" << std::endl;
    std::cout << "  /profiles/pressure_hpa" << std::endl;
    std::cout << "  /profiles/temperature_k" << std::endl;
    std::cout << "  /profiles/dewpoint_k" << std::endl;
    std::cout << "  /profiles/wind_speed_ms" << std::endl;
    std::cout << "  /profiles/wind_direction_deg" << std::endl;
    std::cout << "  /metadata/station_id" << std::endl;
    std::cout << "  /metadata/timestamp_utc" << std::endl;

    // For now, create a sample sounding for testing
    return create_sample_sounding();
}

/*This function reads the SHARPY NetCDF file.
Takes in the file path and reads the SHARPY NetCDF file.*/
SoundingData SharpySoundingScheme::read_sharpy_netcdf(const std::string& file_path) 
{
    // Placeholder implementation - would need NetCDF library
    std::cout << "Reading SHARPY NetCDF file (placeholder implementation)" << std::endl;
    std::cout << "To implement: add NetCDF library dependency and read SHARPY data structure" << std::endl;

    // For now, create a sample sounding for testing
    return create_sample_sounding();
}


/*This function parses the SHARPY profile data.
Takes in the profile data and parses the SHARPY profile data.*/
SoundingData SharpySoundingScheme::parse_sharpy_profile(const std::vector<std::vector<double>>& profile_data) 
{
    // Placeholder for parsing SHARPY profile data structure
    SoundingData data;

    // This would parse the actual SHARPY data format
    // For now, return sample data

    return create_sample_sounding();
}

/*This function creates a sample SHARPY sounding.
Takes in the sample SHARPY sounding.*/
SoundingData SharpySoundingScheme::create_sample_sounding() 
{
    // Create a sample atmospheric sounding for testing
    // This represents a typical mid-latitude summer sounding

    SoundingData data;

    // Station metadata
    data.station_id = "KSAMPLE";
    data.timestamp_utc = "2024-01-15T12:00:00Z";
    data.latitude_deg = 35.0;
    data.longitude_deg = -97.0;
    data.elevation_m = 300.0;

    // Create profile data (surface to 20km)
    const size_t num_levels = 100;
    const double dz = 200.0; // 200m spacing

    // Iterate over the vertical levels and create the sample SHARPY sounding.
    for (size_t i = 0; i < num_levels; ++i) 
    {
        double height = i * dz;

        // Add the height to the sounding data.
        data.height_m.push_back(height);

        // Pressure (exponential decrease)
        double pressure = 1000.0 * std::exp(-height / 8000.0);
        data.pressure_hpa.push_back(pressure);

        // Temperature (dry adiabatic lapse rate near surface, stable aloft)
        double temp_k;

        // If the height is less than 2000 meters, add the temperature to the sounding data.
        if (height < 2000.0) 
        {
            // Boundary layer with slight inversion
            temp_k = 298.0 - 0.006 * height + 2.0 * std::exp(-height / 500.0);
        }
        
        // If the height is less than 12000 meters, add the temperature to the sounding data.
        else if (height < 12000.0) 
        {
            // Unstable layer
            temp_k = 298.0 - 0.006 * 2000.0 - 0.004 * (height - 2000.0);
        }
         else 
         {
            // Stable upper troposphere
            temp_k = 298.0 - 0.006 * 2000.0 - 0.004 * 10000.0 - 0.003 * (height - 12000.0);
        }
        data.temperature_k.push_back(temp_k);

        // Dewpoint (constant mixing ratio near surface, decreasing aloft)
        double dewpoint_k;

        // If the height is less than 3000 meters, add the dewpoint to the sounding data.
        if (height < 3000.0) 
        {
            dewpoint_k = 293.0 - 0.003 * height;
        } 
        else
        {
            dewpoint_k = 293.0 - 0.003 * 3000.0 - 0.002 * (height - 3000.0);
        }
        data.dewpoint_k.push_back(dewpoint_k);

        // Wind (simple hodograph)
        double wind_speed = 5.0 + 25.0 * std::tanh(height / 3000.0);
        double wind_dir = 180.0 + 30.0 * std::sin(height * 3.14159 / 10000.0);
        data.wind_speed_ms.push_back(wind_speed);
        data.wind_direction_deg.push_back(wind_dir);

        // Quality flags (all valid for sample data)
        data.valid_temperature.push_back(true);
        data.valid_moisture.push_back(true);
        data.valid_wind.push_back(true);
    }

    return data;
}
