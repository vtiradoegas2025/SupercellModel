/**
 * @file factory.cpp
 * @brief Implementation for the chaos module.
 *
 * Provides executable logic for the chaos runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/chaos subsystem.
 */

#include "factory.hpp"
#include <algorithm>
#include <cctype>

namespace
{
std::string normalize_scheme_name(std::string scheme_name)
{
    std::transform(scheme_name.begin(), scheme_name.end(), scheme_name.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (scheme_name == "pbl" ||
        scheme_name == "pbl_perturbation" ||
        scheme_name == "pbl_perturbations" ||
        scheme_name == "bl_perturbation" ||
        scheme_name == "bl_perturbations" ||
        scheme_name == "boundary_layer_perturbation" ||
        scheme_name == "boundary_layer_perturbations")
    {
        return "boundary_layer";
    }
    if (scheme_name == "ic")
    {
        return "initial_conditions";
    }
    if (scheme_name == "full")
    {
        return "full_stochastic";
    }
    return scheme_name;
}
}

namespace chaos 
{

std::unique_ptr<ChaosScheme> create_chaos_scheme(const std::string& scheme_name) 
{
    const std::string normalized_name = normalize_scheme_name(scheme_name);

    if (normalized_name == "none") 
    {
        return std::make_unique<NoneScheme>();
    } 

    else if (normalized_name == "initial_conditions") 
    {
        return std::make_unique<InitialConditionsScheme>();
    } 

    else if (normalized_name == "boundary_layer") 
    {
        return std::make_unique<BoundaryLayerScheme>();
    } 

    else if (normalized_name == "full_stochastic") 
    {
        return std::make_unique<FullStochasticScheme>();
    } 

    else 
    {
        throw std::runtime_error("Unknown chaos scheme: " + scheme_name +
                                 " (normalized: " + normalized_name + ")");
    }
}

}

/**
 * @brief Creates the chaos scheme.
 */
std::unique_ptr<chaos::ChaosScheme> create_chaos_scheme(const std::string& scheme_name) 
{
    return chaos::create_chaos_scheme(scheme_name);
}

/**
 * @brief Gets the available chaos schemes.
 */
std::vector<std::string> get_available_chaos_schemes() 
{
    return {"none", "initial_conditions", "boundary_layer", "full_stochastic"};
}
