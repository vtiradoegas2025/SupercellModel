#pragma once
#include <vector>
#include <memory>
#include <string>

// Turbulence module constants
namespace turbulence_constants {
    inline constexpr double g = 9.81;          // gravity [m/s²]
    inline constexpr double kappa = 0.4;       // von Karman constant
    inline constexpr double C_s_default = 0.18; // default Smagorinsky coefficient
    inline constexpr double Pr_t_default = 0.7; // default turbulent Prandtl number
    inline constexpr double Sc_t_default = 0.7; // default turbulent Schmidt number
}

// Configuration for turbulence schemes
struct TurbulenceConfig {
    std::string scheme_id = "smagorinsky";
    double dt_sgs = 1.0;  // seconds (cadence, often = dynamics dt)
    std::string mode = "3d";  // "3d", "horizontal_only", "vertical_only"
    std::string filter_width = "cubic_root";  // "dx", "cubic_root", "user"
    double Cs = turbulence_constants::C_s_default;  // Smagorinsky coefficient
    double Pr_t = turbulence_constants::Pr_t_default;  // turbulent Prandtl number
    double Sc_t = turbulence_constants::Sc_t_default;  // turbulent Schmidt number
    bool dynamic_Cs = false;  // dynamic Smagorinsky (future extension)
    std::string stability_correction = "none";  // "none", "ri", "tke_based"
    bool enable_moist_buoyancy = false;  // for TKE buoyancy production
    bool enable_diagnostics = false;
    double nu_t_max = 1000.0;  // maximum eddy viscosity [m²/s]
    double K_max = 1000.0;     // maximum eddy diffusivity [m²/s]
};

// Grid metrics for turbulence calculations
struct GridMetrics {
    double dx = 1000.0;  // horizontal grid spacing [m]
    double dy = 1000.0;  // horizontal grid spacing [m]
    std::vector<double> dz;  // vertical grid spacing [m] (nz)
    std::vector<double> z_int;  // interface heights [m] (nz+1)
};

// Required fields for turbulence computation
enum class TurbulenceRequirements {
    BASIC = 1,           // u, v, w, rho, theta
    MOISTURE = 2,        // includes qv for buoyancy
    TKE = 4,            // requires TKE field
    STABILITY = 8       // requires stability corrections
};

// Single column/3D state for turbulence computation
struct TurbulenceStateView {
    // Required basic state
    const std::vector<std::vector<std::vector<float>>>* u = nullptr;    // zonal wind [m/s]
    const std::vector<std::vector<std::vector<float>>>* v = nullptr;    // meridional wind [m/s]
    const std::vector<std::vector<std::vector<float>>>* w = nullptr;    // vertical wind [m/s]
    const std::vector<std::vector<std::vector<float>>>* rho = nullptr;  // density [kg/m³]
    const std::vector<std::vector<std::vector<float>>>* theta = nullptr; // potential temperature [K]

    // Optional for advanced schemes
    const std::vector<std::vector<std::vector<float>>>* qv = nullptr;   // water vapor mixing ratio [kg/kg]
    const std::vector<std::vector<std::vector<float>>>* tke = nullptr;  // TKE [m²/s²]

    // Grid information
    int NR = 0, NTH = 0, NZ = 0;  // grid dimensions
};

// Output tendencies and diagnostics
struct TurbulenceTendencies {
    std::vector<std::vector<std::vector<float>>> dudt_sgs;      // u-momentum tendency [m/s²]
    std::vector<std::vector<std::vector<float>>> dvdt_sgs;      // v-momentum tendency [m/s²]
    std::vector<std::vector<std::vector<float>>> dwdt_sgs;      // w-momentum tendency [m/s²]
    std::vector<std::vector<std::vector<float>>> dthetadt_sgs;  // temperature tendency [K/s]
    std::vector<std::vector<std::vector<float>>> dqvdt_sgs;     // moisture tendency [kg/kg/s] (optional)
    std::vector<std::vector<std::vector<float>>> dtkedt_sgs;    // TKE tendency [m²/s³] (TKE scheme)
};

struct TurbulenceDiagnostics {
    // Eddy coefficients (can be 3D or interface-based)
    std::vector<std::vector<std::vector<float>>> nu_t;      // eddy viscosity [m²/s]
    std::vector<std::vector<std::vector<float>>> K_theta;   // temperature diffusivity [m²/s]
    std::vector<std::vector<std::vector<float>>> K_q;       // moisture diffusivity [m²/s]
    std::vector<std::vector<std::vector<float>>> K_tke;     // TKE diffusivity [m²/s] (TKE scheme)

    // TKE budget terms (TKE scheme diagnostics)
    std::vector<std::vector<std::vector<float>>> shear_prod;    // shear production [m²/s³]
    std::vector<std::vector<std::vector<float>>> buoy_prod;     // buoyancy production [m²/s³]
    std::vector<std::vector<std::vector<float>>> dissipation;   // dissipation [m²/s³]

    // Dynamic coefficients (future extension)
    std::vector<std::vector<std::vector<float>>> Cs_eff;    // effective Smagorinsky coefficient
    std::vector<std::vector<std::vector<float>>> Pr_t_eff;  // effective Prandtl number
};

// Abstract base class for turbulence schemes
class TurbulenceSchemeBase {
public:
    virtual ~TurbulenceSchemeBase() = default;

    // Scheme identification and requirements
    virtual std::string name() const = 0;
    virtual int required_fields() const = 0;

    // Initialization
    virtual void initialize(const TurbulenceConfig& cfg) = 0;

    // Core computation
    virtual void compute(
        const TurbulenceConfig& cfg,
        const GridMetrics& grid,
        const TurbulenceStateView& state,
        TurbulenceTendencies& tend,
        TurbulenceDiagnostics* diag_opt = nullptr) = 0;
};

// Factory function type
using TurbulenceSchemeFactory = std::unique_ptr<TurbulenceSchemeBase> (*)(const TurbulenceConfig&);

// Global turbulence scheme instance and configuration
extern std::unique_ptr<TurbulenceSchemeBase> turbulence_scheme;
extern TurbulenceConfig global_turbulence_config;

// Factory and initialization functions
std::unique_ptr<TurbulenceSchemeBase> create_turbulence_scheme(const std::string& scheme_name);
std::vector<std::string> get_available_turbulence_schemes();
void initialize_turbulence(const std::string& scheme_name,
                          const TurbulenceConfig& cfg = TurbulenceConfig{});
void step_turbulence(double current_time, TurbulenceTendencies& tendencies);
