#pragma once

#include <memory>
#include <string>
#include "radar_base.hpp"

/**
 * @brief Factory for creating radar schemes
 *
 * Follows the same pattern as microphysics, turbulence, etc.
 */

class RadarFactory 
{
public:
    /**
     * @brief Create a radar scheme by ID
     *
     * @param scheme_id Scheme identifier ("reflectivity", "velocity", "zdr")
     * @return Unique pointer to the radar scheme
     */
    static std::unique_ptr<RadarSchemeBase> create(const std::string& scheme_id);

    /**
     * @brief Get list of available radar schemes
     */
    static std::vector<std::string> available_schemes();
};
