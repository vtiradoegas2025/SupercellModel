/**
 * @file factory.cpp
 * @brief Implementation for the numerics module.
 *
 * Provides executable logic for the numerics runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/numerics subsystem.
 */

#include "factory.hpp"
#include "schemes/tvd/tvd.hpp"
#include "schemes/weno5/weno5.hpp"
#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace
{
std::string canonicalize_scheme_name(std::string value)
{
    std::string canonical;
    canonical.reserve(value.size());
    for (unsigned char c : value)
    {
        if (c == '_' || c == '-' || std::isspace(c))
        {
            continue;
        }
        canonical.push_back(static_cast<char>(std::tolower(c)));
    }
    return canonical;
}
}

/**
 * @brief Creates the advection scheme.
 */
std::unique_ptr<AdvectionSchemeBase> create_advection_scheme(const std::string& scheme_name) 
{
    const std::string canonical = canonicalize_scheme_name(scheme_name);

    if (canonical == "tvd") 
    {
        return std::make_unique<TVDScheme>();
    }
    else if (canonical == "weno5") 
    {
        return std::make_unique<WENO5Scheme>();
    }
    else 
    {
        throw std::runtime_error("Unknown advection scheme: " + scheme_name);
    }
}

std::vector<std::string> get_available_advection_schemes() 
{
    return {"tvd", "weno5"};
}
