#pragma once

#include <vector>
#include <string>
#include <memory>
#include "simulation.hpp"  // For grid dimensions

/**
 * @brief Radar module for computing synthetic radar observables
 *
 * This module implements AMS-standard radar forward operators for:
 * - Reflectivity (Z_e, Z_dBZ)
 * - Doppler radial velocity (V_r)
 * - Differential reflectivity (Z_DR)
 *.
 */

namespace RadarConstants 
{
    // Radar wavelength - S-band (10 cm) typical for weather radar
    constexpr double lambda = 0.1;  // meters

    // Dielectric factor for water at S-band
    constexpr double K_w_squared = 0.93;

    // Speed of light
    constexpr double c = 2.99792458e8;  // m/s
}

// Forward declarations
struct RadarConfig;
struct RadarStateView;
struct RadarOut;

/**
 * @brief Radar configuration parameters
 */
struct RadarConfig 
{
    std::string scheme_id;  // "reflectivity", "velocity", "zdr"

    // Operator tier selection
    std::string operator_tier;  // reflectivity: "fast_da"|"psd_moment"
                                // zdr: "simple"|"polarimetric_fo"

    // Radar location (for radial velocity and beam geometry)
    double radar_x = 0.0;
    double radar_y = 0.0;
    double radar_z = 0.0;

    // Sampling model
    std::string sampling = "point";  // "point" or "beam_volume" (future)

    // Hydrometeor fields available
    bool has_qr = true;   // rain mixing ratio
    bool has_qs = false;  // snow mixing ratio
    bool has_qg = false;  // graupel mixing ratio
    bool has_qh = false;  // hail mixing ratio
    bool has_qi = false;  // ice mixing ratio

    // Optional number concentrations (for double-moment schemes)
    bool has_Nr = false;  // rain number concentration
    bool has_Ns = false;  // snow number concentration
    bool has_Ng = false;  // graupel number concentration
    bool has_Nh = false;  // hail number concentration
    bool has_Ni = false;  // ice number concentration
};

/**
 * @brief Read-only view of model state for radar calculations
 */
struct RadarStateView 
{
    // Grid dimensions
    int NR, NTH, NZ;

    // Winds (required for velocity)
    const std::vector<std::vector<std::vector<float>>>* u = nullptr;
    const std::vector<std::vector<std::vector<float>>>* v = nullptr;
    const std::vector<std::vector<std::vector<float>>>* w = nullptr;

    // Hydrometeor mixing ratios (required for reflectivity/ZDR)
    const std::vector<std::vector<std::vector<float>>>* qr = nullptr;  // rain
    const std::vector<std::vector<std::vector<float>>>* qs = nullptr;  // snow
    const std::vector<std::vector<std::vector<float>>>* qg = nullptr;  // graupel
    const std::vector<std::vector<std::vector<float>>>* qh = nullptr;  // hail
    const std::vector<std::vector<std::vector<float>>>* qi = nullptr;  // ice

    // Optional number concentrations (for double-moment schemes)
    const std::vector<std::vector<std::vector<float>>>* Nr = nullptr;  // rain
    const std::vector<std::vector<std::vector<float>>>* Ns = nullptr;  // snow
    const std::vector<std::vector<std::vector<float>>>* Ng = nullptr;  // graupel
    const std::vector<std::vector<std::vector<float>>>* Nh = nullptr;  // hail
    const std::vector<std::vector<std::vector<float>>>* Ni = nullptr;  // ice

    // Optional thermodynamics (for temperature-dependent assumptions)
    const std::vector<std::vector<std::vector<float>>>* theta = nullptr;  // potential temperature
    const std::vector<std::vector<std::vector<float>>>* p = nullptr;      // pressure
};

/**
 * @brief Output structure for radar observables
 */
struct RadarOut 
{
    // Grid dimensions
    int NR, NTH, NZ;

    // Primary radar observables
    std::vector<std::vector<std::vector<float>>> Z_dBZ;      // Reflectivity (dBZ)
    std::vector<std::vector<std::vector<float>>> Vr;         // Radial velocity (m/s)
    std::vector<std::vector<std::vector<float>>> ZDR_dB;     // Differential reflectivity (dB)

    // Debug/internal fields (optional but recommended)
    std::vector<std::vector<std::vector<float>>> Ze_linear;  // Linear reflectivity (mm^6/m^3)
    std::vector<std::vector<std::vector<float>>> ZH_dBZ;      // Horizontal reflectivity (dBZ)
    std::vector<std::vector<std::vector<float>>> ZV_dBZ;      // Vertical reflectivity (dBZ)

    // Species contributions (for reflectivity debugging)
    std::vector<std::vector<std::vector<float>>> Ze_rain;     // Rain contribution
    std::vector<std::vector<std::vector<float>>> Ze_snow;     // Snow contribution
    std::vector<std::vector<std::vector<float>>> Ze_graupel;  // Graupel contribution
    std::vector<std::vector<std::vector<float>>> Ze_hail;     // Hail contribution
    std::vector<std::vector<std::vector<float>>> Ze_ice;      // Ice contribution

    /**
     * @brief Initialize all arrays to zero
     */
    void initialize(int nr, int nth, int nz) 
    {
        NR = nr; NTH = nth; NZ = nz;

        auto init_field = [&](auto& field) {
            field.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
        };

        init_field(Z_dBZ);
        init_field(Vr);
        init_field(ZDR_dB);
        init_field(Ze_linear);
        init_field(ZH_dBZ);
        init_field(ZV_dBZ);
        init_field(Ze_rain);
        init_field(Ze_snow);
        init_field(Ze_graupel);
        init_field(Ze_hail);
        init_field(Ze_ice);
    }
};

/**
 * @brief Base class for all radar schemes
 *
 * Follows the same factory pattern as other physics modules.
 */
class RadarSchemeBase 
{
public:
    virtual ~RadarSchemeBase() = default;

    /**
     * @brief Initialize the radar scheme
     */
    virtual void initialize(const RadarConfig& config, int NR, int NTH, int NZ) = 0;

    /**
     * @brief Compute radar observables
     */
    virtual void compute(const RadarConfig& config, const RadarStateView& state, RadarOut& out) = 0;

    /**
     * @brief Get scheme name
     */
    virtual std::string name() const = 0;
};

/**
 * @brief Shared radar geometry and sampling utilities
 */
class RadarGeometry 
{
public:
    /**
     * @brief Compute radar line-of-sight unit vector
     *
     * @param radar_x, radar_y, radar_z Radar location
     * @param x, y, z Sample point
     * @param e_r_unit Output: unit vector from radar to point
     * @param R Output: distance from radar to point
     */
    static void compute_line_of_sight(double radar_x, double radar_y, double radar_z,
                                     double x, double y, double z,
                                     double& e_r_x, double& e_r_y, double& e_r_z,
                                     double& R);

    /**
     * @brief Sample model state at a point (v1: nearest neighbor)
     *
     * Future versions may implement beam-volume averaging.
     */
    static void sample_state_point(const RadarStateView& state, int i, int j, int k,
                                  RadarStateView& point_state);
};
