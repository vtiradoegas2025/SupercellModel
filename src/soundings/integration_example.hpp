/**
 * @file integration_example.hpp
 * @brief Example showing how to integrate sounding system with main program
 *
 * This file demonstrates how the sounding components can be easily plugged
 * into the main SupercellModel program. It shows the minimal changes needed
 * in simulation.hpp and equations.cpp to support sounding-based initialization.
 */

#pragma once

/**
 * EXAMPLE: How to add sounding support to simulation.hpp
 *
 * Add these includes and declarations to simulation.hpp:
 *
 * #include "soundings.hpp"
 *
 * // Sounding configuration
 * extern SoundingConfig global_sounding_config;
 *
 * // Sounding API
 * void initialize_soundings(const SoundingConfig& config);
 * SoundingData load_sounding_data();
 * SoundingData interpolate_sounding_to_grid(const SoundingData& sounding, const std::vector<double>& model_heights);
 */

/**
 * EXAMPLE: How to modify equations.cpp to use soundings
 *
 * In the initialize() function, replace the current thermodynamic profile setup:
 *
 * void initialize() {
 *     // ... existing code ...
 *
 *     // NEW: Try to load sounding data
 *     SoundingData sounding_data = load_sounding_data();
 *
 *     if (sounding_data.is_valid()) {
 *         // Use sounding data for initialization
 *         std::cout << "Using sounding data for thermodynamic initialization" << std::endl;
 *
 *         // Create model height grid
 *         std::vector<double> model_heights;
 *         for (int k = 0; k < NZ; ++k) {
 *             model_heights.push_back(k * dz);
 *         }
 *
 *         // Interpolate sounding to model grid
 *         SoundingData interpolated = interpolate_sounding_to_grid(sounding_data, model_heights);
 *
 *         // Set thermodynamic fields from interpolated sounding
 *         for (int i = 0; i < NR; ++i) {
 *             for (int j = 0; j < NTH; ++j) {
 *                 for (int k = 0; k < NZ; ++k) {
 *                     theta[i][j][k] = interpolated.potential_temperature_k[k];
 *                     qv[i][j][k] = interpolated.mixing_ratio_kgkg[k];
 *                     // Set temperature from potential temperature (would need pressure field)
 *                     // temperature would be calculated from theta and p
 *                 }
 *             }
 *         }
 *     } else {
 *         // EXISTING CODE: Use procedural profiles as fallback
 *         std::cout << "Using procedural thermodynamic profiles" << std::endl;
 *
 *         // ... existing thermodynamic initialization code ...
 *     }
 *
 *     // ... rest of initialization ...
 * }
 */

/**
 * EXAMPLE: How to add sounding configuration to tornado_sim.cpp
 *
 * In the load_config function, add:
 *
 * // Load sounding configuration
 * if (config.count("environment.sounding.scheme_id")) {
 *     global_sounding_config.scheme_id = config["environment.sounding.scheme_id"];
 * }
 * if (config.count("environment.sounding.file_path")) {
 *     global_sounding_config.file_path = config["environment.sounding.file_path"];
 * }
 * if (config.count("environment.sounding.use_fallback_profiles")) {
 *     global_sounding_config.use_fallback_profiles = (config["environment.sounding.use_fallback_profiles"] == "true");
 * }
 *
 * // Initialize soundings
 * initialize_soundings(global_sounding_config);
 */

/**
 * EXAMPLE: YAML configuration additions
 *
 * Add to classic.yaml or other config files:
 *
 * environment:
 *   sounding:
 *     scheme_id: "sharpy"          # "sharpy", "none", or empty
 *     file_path: "data/sounding.h5" # Path to SHARPY file
 *     use_fallback_profiles: true   # Use procedural profiles if sounding fails
 *     interpolation_method: 0       # 0=linear, 1=spline, 2=log-linear
 *     extrapolate_below_ground: false
 *     extrapolate_above_top: false
 */

/**
 * EXAMPLE: Makefile additions
 *
 * Add to Makefile to include sounding sources:
 *
 * SRC += src/soundings/soundings.cpp \
 *        src/soundings/factory.cpp \
 *        src/soundings/base/soundings_base.cpp \
 *        src/soundings/schemes/sharpy/sharpy_sounding.cpp
 *
 * INCLUDE += -Iinclude
 *
 * # Add HDF5/NetCDF libraries when implemented:
 * # LIBS += -lhdf5 -lnetcdf
 */

/**
 * REQUIRED LIBRARIES (for full implementation):
 *
 * To read SHARPY HDF5 files:
 * - HDF5 library (libhdf5-dev on Ubuntu, hdf5 on macOS)
 *
 * To read SHARPY NetCDF files:
 * - NetCDF library (libnetcdf-dev on Ubuntu, netcdf on macOS)
 *
 * Installation commands:
 * Ubuntu: sudo apt-get install libhdf5-dev libnetcdf-dev
 * macOS: brew install hdf5 netcdf
 */
