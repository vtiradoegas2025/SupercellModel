#pragma once

#include "radar_base.hpp"
#include <vector>

/**
 * @brief Reflectivity radar scheme
 *
 * Computes equivalent radar reflectivity factor Z_e and Z_dBZ
 * following AMS standards for Rayleigh scattering.
 *
 * Supports both fast DA-style operators and PSD-moment operators.
 */

class ReflectivityScheme : public RadarSchemeBase {
private:
    // Radar constants
    static constexpr double K_w_squared = 0.93;  // Dielectric factor for water at S-band
    static constexpr double Z_min_dBZ = -30.0;   // Minimum detectable reflectivity

    // Hydrometeor properties (S-band, from AMS literature)
    struct HydrometeorProps {
        double density;      // kg/m³
        double alpha;        // Size distribution parameter
        double c;           // Size distribution parameter
        double d;           // Size distribution parameter
    };

    // Species-specific properties
    static const HydrometeorProps rain_props;
    static const HydrometeorProps snow_props;
    static const HydrometeorProps graupel_props;
    static const HydrometeorProps hail_props;
    static const HydrometeorProps ice_props;

    // Grid dimensions
    int NR_, NTH_, NZ_;

public:
    void initialize(const RadarConfig& config, int NR, int NTH, int NZ) override;

    void compute(const RadarConfig& config, const RadarStateView& state, RadarOut& out) override;

    std::string name() const override { return "reflectivity"; }

private:
    /**
     * @brief Fast DA-style reflectivity operator (single/double moment aware)
     */
    void compute_fast_da(const RadarConfig& config, const RadarStateView& state, RadarOut& out);

    /**
     * @brief PSD-moment reflectivity operator (recommended with moments)
     */
    void compute_psd_moment(const RadarConfig& config, const RadarStateView& state, RadarOut& out);

    /**
     * @brief Compute reflectivity contribution from a single species
     *
     * @param q Mixing ratio field
     * @param Nt Number concentration field (optional, nullptr for single-moment)
     * @param props Hydrometeor properties
     * @param Ze_out Output reflectivity field (linear units)
     */
    void compute_species_reflectivity(
        const std::vector<std::vector<std::vector<float>>>& q,
        const std::vector<std::vector<std::vector<float>>>* Nt,
        const HydrometeorProps& props,
        std::vector<std::vector<std::vector<float>>>& Ze_out);

    /**
     * @brief Compute PSD-moment reflectivity using number concentrations
     *
     * @param q Mixing ratio field
     * @param Nt Number concentration field (required for moment calculations)
     * @param props Hydrometeor properties
     * @param Ze_out Output reflectivity field (linear units)
     * @param species_name Species name for diagnostics
     */
    void compute_moment_reflectivity(
        const std::vector<std::vector<std::vector<float>>>& q,
        const std::vector<std::vector<std::vector<float>>>* Nt,
        const HydrometeorProps& props,
        std::vector<std::vector<std::vector<float>>>& Ze_out,
        const std::string& species_name);

    /**
     * @brief Convert linear reflectivity to dBZ
     */
    static float linear_to_dbz(float Ze_linear) {
        if (Ze_linear <= 1e-10) {
            return Z_min_dBZ;
        }
        return 10.0f * std::log10(Ze_linear);
    }

    /**
     * @brief Convert dBZ to linear reflectivity
     */
    static float dbz_to_linear(float Z_dbz) {
        return std::pow(10.0f, Z_dbz / 10.0f);
    }
};
