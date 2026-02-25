/**
 * @file factory.cpp
 * @brief Implementation for the numerics module.
 *
 * Provides executable logic for the numerics runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/numerics subsystem.
 */

#include "factory.hpp"
#include "schemes/rk3/rk3.hpp"
#include "schemes/rk4/rk4.hpp"
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
 * @brief Creates the time stepping scheme.
 */
std::unique_ptr<TimeSteppingSchemeBase> create_time_stepping_scheme(const std::string& scheme_name) 
{
    const std::string canonical = canonicalize_scheme_name(scheme_name);

    if (canonical == "rk3" || canonical == "ssprk3") 
    {
        return std::make_unique<RK3Scheme>();
    }
    else if (canonical == "rk4") 
    {
        return std::make_unique<RK4Scheme>();
    }
    else 
    {
        throw std::runtime_error("Unknown time stepping scheme: " + scheme_name);
    }
}

/**
 * @brief Gets the available time stepping schemes.
 */
std::vector<std::string> get_available_time_stepping_schemes() 
{
    return {"rk3", "rk4"};
}
