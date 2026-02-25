/**
 * @file factory.cpp
 * @brief Implementation for the dynamics module.
 *
 * Provides executable logic for the dynamics runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/dynamics subsystem.
 */

#include "factory.hpp"
#include "string_utils.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>
#include <stdexcept>

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
 * @brief Normalizes user-provided dynamics scheme aliases.
 */
std::string normalize_dynamics_scheme_name(std::string scheme_name)
{
    scheme_name = tmv::strutil::lower_copy(trim_copy(std::move(scheme_name)));

    if (scheme_name == "mesocyclone")
    {
        return "supercell";
    }
    if (scheme_name == "axisymmetric")
    {
        return "tornado";
    }
    return scheme_name;
}

/**
 * @brief Builds a comma-separated list of available scheme ids.
 */
std::string available_schemes_csv()
{
    const auto schemes = get_available_dynamics_schemes();
    std::ostringstream oss;
    for (std::size_t i = 0; i < schemes.size(); ++i)
    {
        if (i > 0)
        {
            oss << ", ";
        }
        oss << schemes[i];
    }
    return oss.str();
}
}

/**
 * @brief Creates a dynamics scheme instance from its configured id.
 */
std::unique_ptr<DynamicsScheme> create_dynamics_scheme(const std::string& scheme_name) 
{
    const std::string normalized_name = normalize_dynamics_scheme_name(scheme_name);

    if (normalized_name == "supercell")
    {
        return std::make_unique<SupercellScheme>();
    } 
    else if (normalized_name == "tornado")
    {
        return std::make_unique<TornadoScheme>();
    } 
    else 
    {
        throw std::runtime_error("Unknown dynamics scheme: " + scheme_name +
                                 " (normalized: " + normalized_name +
                                 "). Available: " + available_schemes_csv());
    }
}

/**
 * @brief Returns supported dynamics scheme ids.
 */
std::vector<std::string> get_available_dynamics_schemes() 
{
    return {"supercell", "tornado"};
}
