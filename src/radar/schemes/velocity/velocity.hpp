#pragma once

#include "radar_base.hpp"

/**
 * @brief Doppler radial velocity radar scheme
 *
 * Computes radial velocity V_r = V · ê_r following AMS standards.
 *
 * Optional: scatterer motion correction using hydrometeor fall speeds.
 */

class VelocityScheme : public RadarSchemeBase 
{
private:
    // Grid dimensions
    int NR_, NTH_, NZ_;

    // Radar location
    double radar_x_, radar_y_, radar_z_;

public:
    void initialize(const RadarConfig& config, int NR, int NTH, int NZ) override;

    void compute(const RadarConfig& config, const RadarStateView& state, RadarOut& out) override;

    std::string name() const override { return "velocity"; }

private:
    /**
     * @brief Compute radial velocity at each grid point
     */
    void compute_radial_velocity(const RadarConfig& config, const RadarStateView& state, RadarOut& out);

    /**
     * @brief Optional scatterer motion correction
     *
     * V_r_scatterer = V_r_air - V_terminal · ê_r
     */
    void apply_scatterer_correction(const RadarConfig& config, const RadarStateView& state, RadarOut& out);

    /**
     * @brief Estimate terminal fall speed for hydrometeors (simplified)
     *
     * This is a placeholder - real implementation would use proper fall speed relations
     * based on particle size, density, etc.
     */
    float estimate_fall_speed(float q_rain, float q_snow, float q_graupel, float q_hail, float q_ice);
};
