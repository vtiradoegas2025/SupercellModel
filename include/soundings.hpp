#pragma once
#include "soundings_base.hpp"


/*This header file contains the functions for the soundings module.
The soundings module is responsible for the soundings of the simulation.
The soundings scheme is chosen by the user in the configuration file.
This module is used to get the soundings of the simulation.
This module is a placeholder for now and will be implemented in the future.*/

/**
 * @brief Initialize the atmospheric sounding system
 * @param config Sounding configuration
 */
void initialize_soundings(const SoundingConfig& config);

/**
 * @brief Load sounding data for model initialization
 * @return SoundingData structure, or empty if failed and fallback disabled
 */
SoundingData load_sounding_data();

/**
 * @brief Interpolate sounding data to model grid heights
 * @param sounding Original sounding data
 * @param model_heights Vector of model grid heights (m)
 * @return Interpolated SoundingData at model levels
 */
SoundingData interpolate_sounding_to_grid(
    const SoundingData& sounding,
    const std::vector<double>& model_heights);

/**
 * @brief Get current sounding configuration
 * @return Reference to current configuration
 */
const SoundingConfig& get_sounding_config();

/**
 * @brief Check if sounding system is initialized
 * @return true if ready to use
 */
bool is_soundings_initialized();

/**
 * @brief Reset sounding system
 */
void reset_soundings();
