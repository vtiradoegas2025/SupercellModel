/**
 * @file factory.hpp
 * @brief Declarations for the radiation module.
 *
 * Defines interfaces, data structures, and contracts used by
 * the radiation runtime and scheme implementations.
 * This file is part of the src/radiation subsystem.
 */

#pragma once
#include <memory>
#include <string>
#include "radiation_base.hpp"

class SimpleGreyScheme;

/**
 * @brief Creates a radiation scheme by configured name.
 */
std::unique_ptr<RadiationSchemeBase> create_radiation_scheme(const std::string& scheme_name);

/**
 * @brief Returns names of available radiation schemes.
 */
std::vector<std::string> get_available_radiation_schemes();
