/**
 * @file factory.cpp
 * @brief Implementation for the turbulence module.
 *
 * Provides executable logic for the turbulence runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/turbulence subsystem.
 */

#include "factory.hpp"
#include "schemes/smagorinsky/smagorinsky.hpp"
#include "schemes/tke/tke.hpp"
#include <algorithm>
#include <cctype>

namespace
{
/**
 * @brief Normalizes turbulence scheme names and aliases.
 */
std::string normalize_scheme_name(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (value == "smag" || value == "smagorinsky_lilly" || value == "smagorinsky-lilly")
    {
        return "smagorinsky";
    }
    if (value == "1.5" || value == "1.5order" || value == "1.5-order")
    {
        return "tke";
    }
    return value;
}
}

/**
 * @brief Creates a turbulence scheme instance from configured id.
 */
std::unique_ptr<TurbulenceSchemeBase> create_turbulence_scheme(const std::string& scheme_name) 
{
    const std::string normalized = normalize_scheme_name(scheme_name);
    if (normalized == "smagorinsky") 
    {
        return std::make_unique<SmagorinskyScheme>();
    }
    else if (normalized == "tke") 
    {
        return std::make_unique<TKEScheme>();
    }
    else 
    {
        throw std::runtime_error("Unknown turbulence scheme: " + scheme_name);
    }
}

/**
 * @brief Gets the available turbulence schemes.
 */
std::vector<std::string> get_available_turbulence_schemes() 
{
    return {"smagorinsky", "tke"};
}
