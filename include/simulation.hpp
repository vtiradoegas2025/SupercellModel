#pragma once
#include <vector>
#include <memory>
#include "radiation_base.hpp"
#include "boundary_layer_base.hpp"
#include "turbulence_base.hpp"
#include "advection_base.hpp"
#include "diffusion_base.hpp"
#include "time_stepping_base.hpp"
#include "terrain_base.hpp"
#include "chaos_base.hpp"

/*This header file contains the main simulation class and the main simulation functions.
The simulation class is responsible for the main simulation loop.
The simulation functions are responsible for the main simulation functions.
This module is used to initialize the simulation, run the simulation, and clean up the simulation.*/



// Forward declare dynamics scheme
class DynamicsScheme;

// Grid parameters (configurable at runtime)
extern int NR;         // radial points
extern int NTH;        // azimuthal points (theta)
extern int NZ;         // vertical points

// Grid resolution parameters (computed from grid size)
extern double dr;      // radial resolution (m)
extern double dz;      // vertical resolution (m)
extern double dt;      // time step (s)
extern double dtheta;  // azimuthal resolution (rad)

// Physical constants (shared)
inline constexpr double g = 9.81;       // gravity (m/s^2)
inline constexpr double R_d = 287.0;    // dry air gas constant (J/kg·K)
inline constexpr double cp = 1004.0;    // specific heat at constant pressure (J/kg·K)
inline constexpr double gamma = cp / (cp - R_d); // adiabatic index
inline constexpr double p0 = 100000.0;  // reference pressure (Pa)
inline constexpr double theta0 = 300.0; // base potential temperature (K)

// Anelastic approximation parameters
extern std::vector<double> rho0_base;  // base-state density profile (nz levels)

// Kessler microphysics parameters (warm rain)
inline constexpr double qc0 = 1.0e-3;   // autoconversion threshold (kg/kg)
inline constexpr double c_auto = 1.0e-3; // autoconversion rate (s⁻¹)
inline constexpr double c_accr = 2.2;   // accretion coefficient (s⁻¹)
inline constexpr double c_evap = 3.0e-3; // evaporation rate (s⁻¹)
inline constexpr double L_v = 2.5e6;    // latent heat of vaporization (J/kg)
inline constexpr double a_term = 65.0;  // terminal velocity coefficient
inline constexpr double b_term = 0.125; // terminal velocity exponent
inline constexpr double Vt_max = 20.0;  // maximum terminal velocity (m/s)

// Ice microphysics parameters
inline constexpr double T0 = 273.15;    // freezing temperature (K)
inline constexpr double L_f = 3.34e5;   // latent heat of fusion (J/kg)
inline constexpr double L_s = L_v + L_f; // latent heat of sublimation (J/kg)
inline constexpr double c_rime = 1.0;   // riming efficiency
inline constexpr double c_freeze = 1.0e-3; // homogeneous freezing rate (s⁻¹)
inline constexpr double c_melt = 1.0e-3; // melting rate (s⁻¹)
inline constexpr double c_subl = 1.0e-3; // sublimation rate (s⁻¹)
inline constexpr double a_hail = 114.0; // hail terminal velocity coefficient
inline constexpr double b_hail = 0.5;   // hail terminal velocity exponent
inline constexpr double a_grau = 40.0;  // graupel terminal velocity coefficient
inline constexpr double b_grau = 0.37;  // graupel terminal velocity exponent
inline constexpr double Vt_max_hail = 40.0; // max hail terminal velocity (m/s)
inline constexpr double Vt_max_grau = 15.0; // max graupel terminal velocity (m/s)

// Radiation constants
inline constexpr double sigma_sb = 5.67e-8;  // Stefan-Boltzmann (W/m²/K⁴)
inline constexpr double S0_solar = 1366.0;   // Solar constant (W/m²)

// Simulation time tracking
extern double simulation_time;  // current simulation time [s]

// Shared fields (defined in equations.cpp for integrated build; in GUI_STANDALONE, tracer is defined in gui.cpp)
extern std::vector<std::vector<std::vector<float>>> rho;      // [r][th][z]
extern std::vector<std::vector<std::vector<float>>> p;        // [r][th][z]
extern std::vector<std::vector<std::vector<float>>> u;        // radial wind [r][th][z]
extern std::vector<std::vector<std::vector<float>>> w;        // vertical wind [r][th][z]
extern std::vector<std::vector<std::vector<float>>> v_theta;  // azimuthal wind [r][th][z]
extern std::vector<std::vector<std::vector<float>>> tracer;   // passive scalar [r][th][z]

// Thermodynamic fields
extern std::vector<std::vector<std::vector<float>>> theta;    // potential temperature [r][th][z] (K)
extern std::vector<std::vector<std::vector<float>>> qv;       // water vapor mixing ratio [r][th][z] (kg/kg)
extern std::vector<std::vector<std::vector<float>>> qc;       // cloud water mixing ratio [r][th][z] (kg/kg)
extern std::vector<std::vector<std::vector<float>>> qr;       // rain water mixing ratio [r][th][z] (kg/kg)
extern std::vector<std::vector<std::vector<float>>> qi;       // cloud ice mixing ratio [r][th][z] (kg/kg)
extern std::vector<std::vector<std::vector<float>>> qs;       // snow mixing ratio [r][th][z] (kg/kg)
extern std::vector<std::vector<std::vector<float>>> qh;       // hail mixing ratio [r][th][z] (kg/kg)
extern std::vector<std::vector<std::vector<float>>> qg;       // graupel mixing ratio [r][th][z] (kg/kg)

// Wind profile for hodograph initialization
struct WindProfile 
{
    double u_sfc, v_sfc;     // surface winds (m/s)
    double u_1km, v_1km;     // 1km winds (m/s)
    double u_6km, v_6km;     // 6km winds (m/s)
};
extern WindProfile global_wind_profile;
extern double global_cape_target; // J/kg

// Nested grid configuration
struct NestedGridConfig 
{
    bool enabled = false;
    int center_i = NR/2;     // center of nested grid in parent coordinates
    int center_j = NTH/2;
    int center_k = NZ/2;
    int nest_size_r = 32;   // nested grid size
    int nest_size_th = 32;
    int nest_size_z = 32;
    double refinement = 3.0; // refinement ratio
};
extern NestedGridConfig nested_config;

// Nested grid fields (subset of main grid with higher resolution)
extern std::vector<std::vector<std::vector<float>>> nest_rho;
extern std::vector<std::vector<std::vector<float>>> nest_p;
extern std::vector<std::vector<std::vector<float>>> nest_u;
extern std::vector<std::vector<std::vector<float>>> nest_w;
extern std::vector<std::vector<std::vector<float>>> nest_v_theta;
extern std::vector<std::vector<std::vector<float>>> nest_theta;
extern std::vector<std::vector<std::vector<float>>> nest_qv;
extern std::vector<std::vector<std::vector<float>>> nest_qc;
extern std::vector<std::vector<std::vector<float>>> nest_qr;
extern std::vector<std::vector<std::vector<float>>> nest_qh;
extern std::vector<std::vector<std::vector<float>>> nest_qg;
extern std::vector<std::vector<std::vector<float>>> nest_tracer;

// Radar reflectivity field
extern std::vector<std::vector<std::vector<float>>> radar_reflectivity;

// Dynamics scheme
extern std::unique_ptr<DynamicsScheme> dynamics_scheme;

// Vorticity diagnostic fields
extern std::vector<std::vector<std::vector<float>>> vorticity_r;
extern std::vector<std::vector<std::vector<float>>> vorticity_theta;
extern std::vector<std::vector<std::vector<float>>> vorticity_z;
extern std::vector<std::vector<std::vector<float>>> stretching_term;
extern std::vector<std::vector<std::vector<float>>> tilting_term;
extern std::vector<std::vector<std::vector<float>>> baroclinic_term;

// Angular momentum diagnostics (for tornado dynamics)
extern std::vector<std::vector<std::vector<float>>> angular_momentum;
extern std::vector<std::vector<std::vector<float>>> angular_momentum_tendency;

// Pressure diagnostics
extern std::vector<std::vector<std::vector<float>>> p_prime;
extern std::vector<std::vector<std::vector<float>>> dynamic_pressure;
extern std::vector<std::vector<std::vector<float>>> buoyancy_pressure;

// Radiation fields
extern std::vector<std::vector<std::vector<float>>> dtheta_dt_rad;  // radiation tendency

// Boundary layer fields
extern std::vector<std::vector<std::vector<float>>> dtheta_dt_pbl;  // PBL theta tendency
extern std::vector<std::vector<std::vector<float>>> dqv_dt_pbl;     // PBL moisture tendency
extern std::vector<std::vector<std::vector<float>>> du_dt_pbl;      // PBL u-wind tendency
extern std::vector<std::vector<std::vector<float>>> dv_dt_pbl;      // PBL v-wind tendency
extern std::vector<std::vector<std::vector<float>>> dtke_dt_pbl;    // PBL TKE tendency

// TKE field (for MYNN scheme)
extern std::vector<std::vector<std::vector<float>>> tke;
extern RadiationConfig global_radiation_config;

// Boundary layer configuration
extern BoundaryLayerConfig global_boundary_layer_config;
extern SurfaceConfig global_surface_config;

// Turbulence configuration
extern TurbulenceConfig global_turbulence_config;

// Numerics configuration
extern GridMetrics global_grid_metrics;
extern AdvectionConfig global_advection_config;
extern DiffusionConfig global_diffusion_config;
extern TimeSteppingConfig global_time_stepping_config;

// Simulation API
void initialize();
void resize_fields();
void advect_tracer(double dt_advect);
void advect_thermodynamics(double dt_advect);
void step_microphysics(double dt_micro);
void step_dynamics();
void apply_boundary_conditions();
void calculate_radar_reflectivity();

// Microphysics API
void initialize_microphysics(const std::string& scheme_name = "kessler");

// Dynamics API
void initialize_dynamics(const std::string& scheme_name = "supercell");
void step_dynamics(double current_time = 0.0);
void step_dynamics_old(double current_time = 0.0);
void step_dynamics_new(double dt_dynamics, double current_time = 0.0);
void compute_dynamics_diagnostics();

// Radiation API
void step_radiation(double current_time);

// Boundary layer API
void initialize_boundary_layer(const std::string& scheme_name,
                              const BoundaryLayerConfig& cfg,
                              const SurfaceConfig& sfc);
void step_boundary_layer(double current_time);

// Turbulence API
void initialize_turbulence(const std::string& scheme_name,
                          const TurbulenceConfig& cfg);
void step_turbulence(double current_time, TurbulenceTendencies& tendencies);

// Chaos API (stochastic perturbations)
void initialize_chaos(const chaos::ChaosConfig& cfg);
void apply_chaos_initial_conditions();
void apply_chaos_tendencies();
void step_chaos_noise(double dt);
const chaos::ChaosDiagnostics& get_chaos_diagnostics();
void reset_chaos_diagnostics();

// Numerics API
void initialize_numerics();

// Terrain API
void build_terrain_fields();

// Nested grid API
void initialize_nested_grid();
void feedback_to_parent();

// GUI API (optional)
#ifdef ENABLE_GUI
void run_gui(int autoExportMs = 0);
#endif
