/**
 * @file factory.cpp
 * @brief Implementation for the radar module.
 *
 * Provides executable logic for the radar runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/radar subsystem.
 */

#include "factory.hpp"
#include <stdexcept>

#include "schemes/reflectivity/reflectivity.hpp"
#include "schemes/velocity/velocity.hpp"
#include "schemes/zdr/zdr.hpp"
/**
 * @brief Creates the radar scheme.
 */

std::unique_ptr<RadarSchemeBase> RadarFactory::create(const std::string& scheme_id) 
{
    if (scheme_id == "reflectivity") 
    {
        return std::make_unique<ReflectivityScheme>();
    } 
    else if (scheme_id == "velocity") 
    {
        return std::make_unique<VelocityScheme>();
    } 
    else if (scheme_id == "zdr") 
    {
        return std::make_unique<ZDRScheme>();
    } 
    else 
    {
        throw std::runtime_error("Unknown radar scheme: " + scheme_id);
    }
}

/**
 * @brief Gets the available radar schemes.
 */
std::vector<std::string> RadarFactory::available_schemes() 
{
    return {"reflectivity", "velocity", "zdr"};
}
