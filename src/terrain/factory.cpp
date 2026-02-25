/**
 * @file factory.cpp
 * @brief Implementation for the terrain module.
 *
 * Provides executable logic for the terrain runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/terrain subsystem.
 */

#include "factory.hpp"
#include "schemes/bell/bell.hpp"
#include "schemes/schar/schar.hpp"
#include "schemes/none.hpp"
#include "string_utils.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>

namespace
{
/**
 * @brief Trims leading and trailing ASCII whitespace.
 */
std::string trim_copy(std::string value)
{
    const auto first = std::find_if_not(value.begin(), value.end(), [](unsigned char c) { return std::isspace(c) != 0; });
    if (first == value.end())
    {
        return "";
    }
    const auto last = std::find_if_not(value.rbegin(), value.rend(), [](unsigned char c) { return std::isspace(c) != 0; }).base();
    return std::string(first, last);
}

/**
 * @brief Normalizes terrain scheme names and aliases.
 */
std::string normalize_terrain_scheme_name(std::string scheme_name)
{
    scheme_name = tmv::strutil::lower_copy(trim_copy(std::move(scheme_name)));
    if (scheme_name.empty() || scheme_name == "flat")
    {
        return "none";
    }
    if (scheme_name == "schaer")
    {
        return "schar";
    }
    return scheme_name;
}
}

/**
 * @brief Creates a terrain scheme instance from configured id.
 */
std::unique_ptr<TerrainSchemeBase> create_terrain_scheme(const std::string& scheme_name) 
{
    const std::string normalized_name = normalize_terrain_scheme_name(scheme_name);

    if (normalized_name == "bell") 
    {
        return std::make_unique<BellScheme>();
    }
    else if (normalized_name == "schar") 
    {
        return std::make_unique<ScharScheme>();
    }
    else if (normalized_name == "none") {
        return std::make_unique<NoneScheme>();
    }
    else 
    {
        std::cerr << "Warning: Unknown terrain scheme '" << scheme_name
                  << "'. Falling back to 'none'." << std::endl;
        return std::make_unique<NoneScheme>();
    }
}

/**
 * @brief Gets the available terrain schemes.
 */
std::vector<std::string> get_available_terrain_schemes() 
{
    return {"none", "bell", "schar"};
}
