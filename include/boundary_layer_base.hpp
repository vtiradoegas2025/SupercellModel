#pragma once
#include <vector>
#include <memory>
#include <string>

/*This header file contains the base classes and structures for the boundary layer module.
The boundary layer module is responsible for the boundary layer physics in the simulation.
The boundary layer scheme is chosen by the user in the configuration file.
This module is used to compute the boundary layer physics of the simulation.*/

// Boundary layer module constants
namespace boundary_layer_constants 
{
    inline constexpr double kappa = 0.4;       // von Karman constant
    inline constexpr double g = 9.81;          // gravity [m/s²]
    inline constexpr double cp = 1004.0;      // specific heat [J/kg/K]
    inline constexpr double T0 = 273.15;      // freezing temperature [K]
    inline constexpr double pi = 3.14159265358979323846;
}

// Configuration for boundary layer schemes
struct BoundaryLayerConfig 
{
    std::string scheme_id = "slab";
    double dt_pbl = 60.0;  // seconds (cadence)
    bool apply_surface_fluxes = true;
    std::string surface_layer_id = "bulk";  // "bulk" or "monin_obukhov"
    bool clear_sky_only = true;
    bool enable_nonlocal_transport = false;  // YSU-style
    bool enable_tke = false;                 // MYNN-style
    std::string pbl_top_method = "ri";       // "ri", "theta", or "tke"

    // Scheme-specific parameters
    double pbl_max_height = 2000.0;  // maximum PBL height [m]
    double min_ustar = 1e-4;         // minimum friction velocity [m/s]
};

// Surface properties configuration
struct SurfaceConfig 
{
    double z0m = 0.1;      // roughness length for momentum [m]
    double z0h = 0.01;     // roughness length for heat [m]
    double albedo = 0.2;   // surface albedo [0-1]
    double emissivity = 0.95; // surface emissivity [0-1]

    // Surface state (can be prescribed or computed)
    double Tsfc = 288.0;   // surface temperature [K]
    double qsfc = 0.01;    // surface specific humidity [kg/kg]
};

// Required fields for boundary layer computation
enum class BoundaryLayerRequirements 
{
    BASIC = 1,           // u, v, theta, qv, p, rho, z
    MOISTURE = 2,        // includes moisture effects
    TKE = 4,            // requires TKE field
    CLOUDS = 8          // includes cloud effects
};

// Single column state for boundary layer computation
struct BoundaryLayerColumnStateView 
{
    // Required basic state (nz = number of layers)
    const std::vector<double>* u = nullptr;        // zonal wind [m/s] (nz)
    const std::vector<double>* v = nullptr;        // meridional wind [m/s] (nz)
    const std::vector<double>* theta = nullptr;    // potential temperature [K] (nz)
    const std::vector<double>* qv = nullptr;       // water vapor mixing ratio [kg/kg] (nz)
    const std::vector<double>* p = nullptr;        // pressure [Pa] (nz)
    const std::vector<double>* rho = nullptr;      // density [kg/m³] (nz)
    const std::vector<double>* z_int = nullptr;    // interface heights [m] (nz+1)

    // Optional for advanced schemes
    const std::vector<double>* tke = nullptr;      // turbulent kinetic energy [m²/s²] (nz)
    const std::vector<double>* qc = nullptr;       // cloud water mixing ratio (nz)
    const std::vector<double>* qi = nullptr;       // cloud ice mixing ratio (nz)

    // Surface layer inputs (lowest model level)
    double u_sfc = 0.0;     // surface zonal wind [m/s]
    double v_sfc = 0.0;     // surface meridional wind [m/s]
    double theta_sfc = 288.0; // surface potential temperature [K]
    double qv_sfc = 0.01;   // surface water vapor [kg/kg]
    double z_sfc = 10.0;    // height of lowest model level [m]
};

// Output tendencies and diagnostics
struct BoundaryLayerTendencies 
{
    std::vector<double> dudt_pbl;      // momentum tendency u [m/s²] (nz)
    std::vector<double> dvdt_pbl;      // momentum tendency v [m/s²] (nz)
    std::vector<double> dthetadt_pbl;  // temperature tendency [K/s] (nz)
    std::vector<double> dqvdt_pbl;     // moisture tendency [kg/kg/s] (nz)
    std::vector<double> dtkedt_pbl;    // TKE tendency [m²/s³] (nz) - for MYNN
};

struct BoundaryLayerDiagnostics 
{
    // Eddy diffusivities at interfaces (nz+1)
    std::vector<double> K_m;        // momentum diffusivity [m²/s]
    std::vector<double> K_h;        // heat/moisture diffusivity [m²/s]
    std::vector<double> K_e;        // TKE diffusivity [m²/s] - for MYNN

    // PBL characteristics
    double pbl_height = 0.0;        // PBL height [m]
    double ustar = 0.0;            // friction velocity [m/s]

    // Surface fluxes
    double tau_u = 0.0;            // zonal momentum stress [Pa]
    double tau_v = 0.0;            // meridional momentum stress [Pa]
    double H = 0.0;                // sensible heat flux [W/m²]
    double E = 0.0;                // moisture flux [kg/m²/s] (latent heat flux / L_v)

    // Stability diagnostics
    double L = 0.0;                // Obukhov length [m]
    double rib_bulk = 0.0;         // bulk Richardson number

    // Transfer coefficients (optional diagnostics)
    double Cd = 0.0;               // drag coefficient
    double Ch = 0.0;               // heat transfer coefficient
    double Ce = 0.0;               // moisture transfer coefficient
};

// Abstract base class for boundary layer schemes
class BoundaryLayerSchemeBase 
{
public:
    virtual ~BoundaryLayerSchemeBase() = default;

    // Scheme identification and requirements
    virtual std::string name() const = 0;
    virtual int required_fields() const = 0;

    // Initialization
    virtual void initialize(const BoundaryLayerConfig& cfg) = 0;

    // Core computation
    virtual void compute_column(
        const BoundaryLayerConfig& cfg,
        const SurfaceConfig& sfc,
        const BoundaryLayerColumnStateView& col,
        BoundaryLayerTendencies& tend,
        BoundaryLayerDiagnostics* diag_opt = nullptr) = 0;
};

// Factory function type
using BoundaryLayerSchemeFactory = std::unique_ptr<BoundaryLayerSchemeBase> (*)(const BoundaryLayerConfig&);

// Global boundary layer scheme instance and configuration
extern std::unique_ptr<BoundaryLayerSchemeBase> boundary_layer_scheme;
extern BoundaryLayerConfig global_boundary_layer_config;
extern SurfaceConfig global_surface_config;

// Factory and initialization functions
std::unique_ptr<BoundaryLayerSchemeBase> create_boundary_layer_scheme(const std::string& scheme_name);
std::vector<std::string> get_available_boundary_layer_schemes();
void initialize_boundary_layer(const std::string& scheme_name,
                              const BoundaryLayerConfig& cfg = BoundaryLayerConfig{},
                              const SurfaceConfig& sfc = SurfaceConfig{});
void step_boundary_layer(double current_time);
