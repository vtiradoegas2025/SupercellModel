/**
 * @file factory.cpp
 * @brief Implementation for the microphysics module.
 *
 * Provides executable logic for the microphysics runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/microphysics subsystem.
 */

#include "factory.hpp"
#include "schemes/kessler/kessler.hpp"
#include "schemes/lin/lin.hpp"
#include "schemes/thompson/thompson.hpp"
#include "schemes/milbrandt/milbrandt.hpp"

/**
 * @brief Creates a microphysics scheme instance from its configured id.
 */
std::unique_ptr<MicrophysicsScheme> create_microphysics_scheme(const std::string& scheme_name) 
{
    if (scheme_name == "kessler") 
    {
        return std::make_unique<KesslerScheme>();
    } else if (scheme_name == "lin") 
    {
        return std::make_unique<LinScheme>();
    } 
    else if (scheme_name == "thompson") 
    {
        return std::make_unique<ThompsonScheme>();
    } 
    else if (scheme_name == "milbrandt")
    {
        return std::make_unique<MilbrandtScheme>();
    } 
    else 
    {
        throw std::runtime_error("Unknown microphysics scheme: " + scheme_name);
    }
}

/**
 * @brief Returns supported microphysics scheme ids.
 */
std::vector<std::string> get_available_schemes() 
{
    return {"kessler", "lin", "thompson", "milbrandt"};
}
