/**
 * @file test_soundings.cpp
 * @brief Simple test program for the soundings module
 *
 * This test demonstrates the soundings module functionality without
 Placeholder file for the future and this was AI generated for a test.
 * requiring integration with the full SupercellModel.
 */

#include "../../include/soundings.hpp"
#include <iostream>
#include <vector>

int main() {
    std::cout << "Testing Soundings Module" << std::endl;
    std::cout << "========================" << std::endl;

    try {
        // Test 1: Initialize with sample sounding
        std::cout << "\nTest 1: Initialize sounding scheme" << std::endl;
        SoundingConfig config;
        config.scheme_id = "sharpy";
        config.file_path = "";  // Empty path will use sample data
        config.use_fallback_profiles = true;

        initialize_soundings(config);
        std::cout << "✓ Sounding scheme initialized successfully" << std::endl;

        // Test 2: Load sounding data
        std::cout << "\nTest 2: Load sounding data" << std::endl;
        SoundingData data = load_sounding_data();

        if (data.is_valid()) {
            std::cout << "✓ Sounding data loaded successfully" << std::endl;
            std::cout << "  Station: " << data.station_id << std::endl;
            std::cout << "  Levels: " << data.num_levels() << std::endl;
            std::cout << "  Height range: " << data.height_m.front()
                      << " - " << data.height_m.back() << " m" << std::endl;
            std::cout << "  Pressure range: " << data.pressure_hpa.front()
                      << " - " << data.pressure_hpa.back() << " hPa" << std::endl;
        } else {
            std::cout << "✗ Failed to load sounding data" << std::endl;
            return 1;
        }

        // Test 3: Interpolate to model grid
        std::cout << "\nTest 3: Interpolate to model grid" << std::endl;
        std::vector<double> model_heights;
        const double dz = 100.0;  // 100m spacing
        const int nz = 50;        // 50 levels
        for (int k = 0; k < nz; ++k) {
            model_heights.push_back(k * dz);
        }

        SoundingData interpolated = interpolate_sounding_to_grid(data, model_heights);

        if (interpolated.is_valid() && interpolated.num_levels() == nz) {
            std::cout << "✓ Interpolation successful" << std::endl;
            std::cout << "  Interpolated to " << interpolated.num_levels() << " levels" << std::endl;
            std::cout << "  Sample values at level 10:" << std::endl;
            std::cout << "    Height: " << interpolated.height_m[10] << " m" << std::endl;
            std::cout << "    Pressure: " << interpolated.pressure_hpa[10] << " hPa" << std::endl;
            std::cout << "    Temperature: " << interpolated.temperature_k[10] << " K" << std::endl;
            std::cout << "    Potential Temp: " << interpolated.potential_temperature_k[10] << " K" << std::endl;
            if (!interpolated.mixing_ratio_kgkg.empty()) {
                std::cout << "    Mixing Ratio: " << interpolated.mixing_ratio_kgkg[10] << " kg/kg" << std::endl;
            }
        } else {
            std::cout << "✗ Interpolation failed" << std::endl;
            return 1;
        }

        // Test 4: Test with invalid scheme
        std::cout << "\nTest 4: Test invalid scheme handling" << std::endl;
        reset_soundings();

        SoundingConfig invalid_config;
        invalid_config.scheme_id = "invalid_scheme";

        try {
            initialize_soundings(invalid_config);
            std::cout << "✗ Should have thrown exception for invalid scheme" << std::endl;
            return 1;
        } catch (const std::runtime_error& e) {
            std::cout << "✓ Correctly caught exception: " << e.what() << std::endl;
        }

        // Test 5: Test "none" scheme
        std::cout << "\nTest 5: Test 'none' scheme" << std::endl;
        SoundingConfig none_config;
        none_config.scheme_id = "none";

        initialize_soundings(none_config);
        if (!is_soundings_initialized()) {
            std::cout << "✓ 'none' scheme correctly returns null implementation" << std::endl;
        } else {
            std::cout << "✗ 'none' scheme should not be initialized" << std::endl;
            return 1;
        }

        std::cout << "\n🎉 All tests passed!" << std::endl;
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Test failed with exception: " << e.what() << std::endl;
        return 1;
    }
}
