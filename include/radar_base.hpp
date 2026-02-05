#pragma once

#include <vector>
#include <string>
#include <memory>
#include "field3d.hpp"
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
    const Field3D* u = nullptr;
    const Field3D* v = nullptr;
    const Field3D* w = nullptr;

    // Hydrometeor mixing ratios (required for reflectivity/ZDR)
    const Field3D* qr = nullptr;  // rain
    const Field3D* qs = nullptr;  // snow
    const Field3D* qg = nullptr;  // graupel
    const Field3D* qh = nullptr;  // hail
    const Field3D* qi = nullptr;  // ice

    // Optional number concentrations (for double-moment schemes)
    const Field3D* Nr = nullptr;  // rain
    const Field3D* Ns = nullptr;  // snow
    const Field3D* Ng = nullptr;  // graupel
    const Field3D* Nh = nullptr;  // hail
    const Field3D* Ni = nullptr;  // ice

    // Optional thermodynamics (for temperature-dependent assumptions)
    const Field3D* theta = nullptr;  // potential temperature
    const Field3D* p = nullptr;      // pressure
};

/**
 * @brief Output structure for radar observables
 */
struct RadarOut 
{
    // Grid dimensions
    int NR, NTH, NZ;

    // Primary radar observables
    Field3D Z_dBZ;      // Reflectivity (dBZ)
    Field3D Vr;         // Radial velocity (m/s)
    Field3D ZDR_dB;     // Differential reflectivity (dB)

    // Debug/internal fields (optional but recommended)
    Field3D Ze_linear;  // Linear reflectivity (mm^6/m^3)
    Field3D ZH_dBZ;      // Horizontal reflectivity (dBZ)
    Field3D ZV_dBZ;      // Vertical reflectivity (dBZ)

    // Species contributions (for reflectivity debugging)
    Field3D Ze_rain;     // Rain contribution
    Field3D Ze_snow;     // Snow contribution
    Field3D Ze_graupel;  // Graupel contribution
    Field3D Ze_hail;     // Hail contribution
    Field3D Ze_ice;      // Ice contribution

    /**
     * @brief Initialize all arrays to zero
     */
    void initialize(int nr, int nth, int nz) 
    {
        NR = nr; NTH = nth; NZ = nz;

        Z_dBZ.resize(nr, nth, nz, 0.0f);
        Vr.resize(nr, nth, nz, 0.0f);
        ZDR_dB.resize(nr, nth, nz, 0.0f);
        Ze_linear.resize(nr, nth, nz, 0.0f);
        ZH_dBZ.resize(nr, nth, nz, 0.0f);
        ZV_dBZ.resize(nr, nth, nz, 0.0f);
        Ze_rain.resize(nr, nth, nz, 0.0f);
        Ze_snow.resize(nr, nth, nz, 0.0f);
        Ze_graupel.resize(nr, nth, nz, 0.0f);
        Ze_hail.resize(nr, nth, nz, 0.0f);
        Ze_ice.resize(nr, nth, nz, 0.0f);
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
