#pragma once
#include <vector>
#include <memory>
#include <string>

/*This header file contains the base classes and structures for the radiation module.
The radiation module is responsible for the radiation of the simulation.
The radiation scheme is chosen by the user in the configuration file.
This module is used to compute the radiation fluxes, column state, and column tendencies of 
the simulation.*/

// Radiation module constants
namespace radiation_constants 
{
    inline constexpr double sigma = 5.67e-8;      // Stefan-Boltzmann constant (W/m²/K⁴)
    inline constexpr double S0 = 1366.0;          // Solar constant (W/m²)
    inline constexpr double pi = 3.14159265358979323846;
}

// Configuration for radiation schemes
struct RadiationConfig 
{
    std::string scheme_id = "simple_grey";
    bool do_lw = true;
    bool do_sw = true;
    double dt_radiation = 300.0;  // seconds (cadence)
    bool clear_sky_only = true;
    bool use_subcycling = false;
    bool enable_diagnostics = false;

    // Simple grey parameters
    double tau_lw_ref = 6.0;      // reference LW optical depth
    double tau_sw_ref = 0.22;     // reference SW optical depth
    double n_lw = 4.0;            // LW optical depth exponent
    double n_sw = 2.0;            // SW optical depth exponent
};

// Required fields for radiation computation
enum class RadiationRequirements 
{
    TEMPERATURE = 1,      // needs T[nz]
    THETA_PRESSURE = 2,   // needs theta[nz] + p[nz]
    SOLAR_GEOMETRY = 4,   // needs mu0, S0
    CLOUD_OPTICS = 8,     // needs cloud properties
    AEROSOL_OPTICS = 16   // needs aerosol properties
};

// Single column state for radiation computation
struct RadiationColumnStateView 
{
    // Required thermodynamic state (nz = number of layers)
    const std::vector<double>* T = nullptr;        // temperature [K] (nz)
    const std::vector<double>* theta = nullptr;    // potential temperature [K] (nz)
    const std::vector<double>* p = nullptr;        // pressure [Pa] (nz)
    const std::vector<double>* rho = nullptr;      // density [kg/m³] (nz)
    const std::vector<double>* dz = nullptr;       // layer thickness [m] (nz)
    const std::vector<double>* z_int = nullptr;    // interface heights [m] (nz+1)

    // Solar geometry (if SW enabled)
    double mu0 = 0.0;    // cos(solar zenith angle)
    double S0 = 1366.0;  // solar constant [W/m²]

    // Surface boundary conditions
    double Tsfc = 288.0;      // surface temperature [K]
    double albedo_sw = 0.2;   // shortwave albedo [0-1]
    double emissivity_lw = 0.95; // longwave emissivity [0-1]

    // Optional: water vapor for more realistic optics
    const std::vector<double>* qv = nullptr;       // water vapor mixing ratio (nz)

    // Future: cloud/aerosol properties
    // const std::vector<double>* cldfra = nullptr;   // cloud fraction
    // const std::vector<double>* re_liq = nullptr;   // liquid effective radius
    // const std::vector<double>* re_ice = nullptr;   // ice effective radius
};

// Output tendencies and fluxes
struct RadiationColumnTendencies 
{
    std::vector<double> dTdt_rad;      // total radiative heating [K/s] (nz)
    std::vector<double> dTdt_lw;       // LW heating [K/s] (nz)
    std::vector<double> dTdt_sw;       // SW heating [K/s] (nz)
};

struct RadiationColumnFluxes {
    std::vector<double> Fup_lw;        // LW upward flux [W/m²] (nz+1)
    std::vector<double> Fdn_lw;        // LW downward flux [W/m²] (nz+1)
    std::vector<double> Fup_sw;        // SW upward flux [W/m²] (nz+1)
    std::vector<double> Fdn_sw;        // SW downward flux [W/m²] (nz+1)

    // Convenience computed fields
    std::vector<double> Fnet_lw;       // LW net flux [W/m²] (nz+1)
    std::vector<double> Fnet_sw;       // SW net flux [W/m²] (nz+1)
    std::vector<double> Fnet_total;    // total net flux [W/m²] (nz+1)
};

// Abstract base class for radiation schemes
class RadiationSchemeBase 
{
public:
    virtual ~RadiationSchemeBase() = default;

    // Scheme identification and requirements
    virtual std::string name() const = 0;
    virtual int required_fields() const = 0;

    // Initialization
    virtual void initialize(const RadiationConfig& cfg) = 0;

    // Core computation
    virtual void compute_column(
        const RadiationConfig& cfg,
        const RadiationColumnStateView& col,
        RadiationColumnTendencies& tend,
        RadiationColumnFluxes* fluxes_opt = nullptr) = 0;
};

// Factory function type
using RadiationSchemeFactory = std::unique_ptr<RadiationSchemeBase> (*)(const RadiationConfig&);

// Global radiation scheme instance
extern std::unique_ptr<RadiationSchemeBase> radiation_scheme;
extern RadiationConfig global_radiation_config;

// Factory and initialization functions
std::unique_ptr<RadiationSchemeBase> create_radiation_scheme(const std::string& scheme_name);
std::vector<std::string> get_available_radiation_schemes();
void initialize_radiation(const std::string& scheme_name, const RadiationConfig& cfg = RadiationConfig{});
