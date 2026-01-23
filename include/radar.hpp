#pragma once

/**
 * @file radar.hpp
 * @brief Main radar module interface
 *
 * This header provides the main interface for the radar forward operator system.
 * It implements AMS-standard radar observables (Z_e, V_r, Z_DR) from model state.
 *
 * Key components:
 * - RadarSchemeBase: Abstract base class for all radar schemes
 * - RadarFactory: Factory for creating radar schemes
 * - RadarConfig: Configuration structure
 * - RadarStateView: Read-only view of model state
 * - RadarOut: Output structure for radar observables
 * - RadarGeometry: Shared geometry and sampling utilities
 *
 * Usage:
 * @code
 * #include "radar.hpp"
 *
 * // Create radar scheme
 * auto radar = RadarFactory::create("reflectivity");
 *
 * // Configure
 * RadarConfig config;
 * config.scheme_id = "reflectivity";
 * config.operator_tier = "fast_da";
 *
 * // Initialize
 * radar->initialize(config, NR, NTH, NZ);
 *
 * // Compute observables
 * RadarOut output;
 * output.initialize(NR, NTH, NZ);
 * radar->compute(config, state_view, output);
 * @endcode
 */

#include "radar_base.hpp"
#include <memory>

/**
 * @brief Radar module for computing synthetic radar observables
 *
 * This module provides AMS-standard forward operators for:
 * - Equivalent radar reflectivity factor (Z_e, dBZ)
 * - Radial Doppler velocity (V_r, m/s)
 * - Differential reflectivity (Z_DR, dB)
 *
 * The system is designed with multiple tiers for different fidelity levels:
 * - Basic/fast operators for real-time applications
 * - Research-grade operators following AMS literature
 * - Extensible framework for future enhancements
 */

// Forward declarations for external users
class RadarSchemeBase;
struct RadarConfig;
struct RadarStateView;
struct RadarOut;
class RadarGeometry;

/**
 * @brief Factory for creating radar schemes
 */
class RadarFactory;

/**
 * @brief High-level radar computation interface
 *
 * Convenience functions for common radar calculations.
 */
class RadarSystem 
{
public:
    /**
     * @brief Compute all radar observables at once
     *
     * This is a convenience function that creates and runs all three
     * radar schemes (reflectivity, velocity, zdr) with appropriate configurations.
     *
     * @param state Read-only view of model state
     * @param output Output structure (must be pre-initialized)
     * @param radar_x Radar location x-coordinate (m)
     * @param radar_y Radar location y-coordinate (m)
     * @param radar_z Radar location z-coordinate (m)
     */
    static void compute_all_observables(const RadarStateView& state, RadarOut& output,
                                       double radar_x = 0.0, double radar_y = 0.0, double radar_z = 0.0);

    /**
     * @brief Compute radar reflectivity only
     */
    static void compute_reflectivity(const RadarStateView& state, RadarOut& output);

    /**
     * @brief Compute Doppler velocity only
     */
    static void compute_velocity(const RadarStateView& state, RadarOut& output,
                                double radar_x = 0.0, double radar_y = 0.0, double radar_z = 0.0);

    /**
     * @brief Compute differential reflectivity only
     */
    static void compute_zdr(const RadarStateView& state, RadarOut& output);
};

/**
 * @brief Radar scheme identification and metadata
 */
namespace RadarSchemes 
{
    static const std::string REFLECTIVITY = "reflectivity";
    static const std::string VELOCITY = "velocity";
    static const std::string ZDR = "zdr";
}

/**
 * @brief Operator tiers for different fidelity levels
 */
namespace RadarTiers 
{
    // Reflectivity operator tiers
    namespace Reflectivity 
    {
        static const std::string FAST_DA = "fast_da";           // Fast data assimilation style
        static const std::string PSD_MOMENT = "psd_moment";     // PSD-based with moments
    }

    // ZDR operator tiers
    namespace ZDR 
    {
        static const std::string SIMPLE = "simple";             // Parameterized axis ratios
        static const std::string POLARIMETRIC_FO = "polarimetric_fo";  // Full forward operator
    }

    // Velocity operator tiers (future)
    namespace Velocity 
    {
        static const std::string POINT = "point";              // Point sampling
        static const std::string BEAM_VOLUME = "beam_volume";  // Volume averaging (future)
    }
}
