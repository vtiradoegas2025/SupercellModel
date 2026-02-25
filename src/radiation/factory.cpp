/**
 * @file factory.cpp
 * @brief Implementation for the radiation module.
 *
 * Provides executable logic for the radiation runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/radiation subsystem.
 */

#include "factory.hpp"
#include "schemes/simple_grey/simple_grey.hpp"

/**
 * @brief Creates the radiation scheme.
 */
std::unique_ptr<RadiationSchemeBase> create_radiation_scheme(const std::string& scheme_name) 
{
    if (scheme_name == "simple_grey") {
        return std::make_unique<SimpleGreyScheme>();
    }
    else 
    {
        throw std::runtime_error("Unknown radiation scheme: " + scheme_name);
    }
}

/**
 * @brief Gets the available radiation schemes.
 */
std::vector<std::string> get_available_radiation_schemes() 
{
    return {"simple_grey"};
}
