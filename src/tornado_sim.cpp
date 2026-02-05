#include <vector>
#include <string>
#include <algorithm>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <iostream>
#include <chrono>
#include <fstream>
#include <unordered_map>
#include "simulation.hpp"


/*This file contains the implementation of the tornado simulation.
It manages the initialization of the tornado simulation and the computation of the tornado simulation.*/

// Forward declarations for functions defined in equations.cpp
extern void initialize_radar(const std::string& scheme_name);


// Compute Weisman-Klemp style wind profile
/*This function computes the wind profile.
Takes in the wind profile and computes the wind profile.*/

void compute_wind_profile(const WindProfile& profile, double z, double& u, double& v) 
{
    // Define the heights in meters
    const double z_sfc = 0.0;
    const double z_1km = 1000.0;
    const double z_6km = 6000.0;

    // If the height is less than or equal to 1km, interpolate the wind profile from the surface to 1km.
    if (z <= z_1km) 
    {
        // Linear interpolation from surface to 1km
        double frac = (z - z_sfc) / (z_1km - z_sfc);
        u = profile.u_sfc + frac * (profile.u_1km - profile.u_sfc);
        v = profile.v_sfc + frac * (profile.v_1km - profile.v_sfc);
    } 

    // If the height is between 1km and 6km, interpolate the wind profile from 1km to 6km.
    else if (z <= z_6km) 
    {
        // Linear interpolation from 1km to 6km
        double frac = (z - z_1km) / (z_6km - z_1km);
        u = profile.u_1km + frac * (profile.u_6km - profile.u_1km);
        v = profile.v_1km + frac * (profile.v_6km - profile.v_1km);
    } 
    else 
    {
        // Constant above 6km
        u = profile.u_6km;
        v = profile.v_6km;
    }
}

/*This function parses a YAML file.
Takes in the filename and parses the YAML file.*/
std::unordered_map<std::string, std::string> parse_yaml_simple(const std::string& filename)
{
    std::unordered_map<std::string, std::string> config;
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Could not open config file: " << filename << std::endl;
        return config;
    }

    std::string line;
    std::vector<std::string> section_stack;

    // Iterate over the lines in the YAML file and parse the YAML file.
    while (std::getline(file, line))
    {
        // Remove comments
        size_t comment_pos = line.find('#');
        if (comment_pos != std::string::npos)
        {
            line = line.substr(0, comment_pos);
        }

        // Count leading spaces BEFORE trimming
        size_t indent = 0;
        while (indent < line.size() && line[indent] == ' ') indent++;

        // Convert spaces to indentation level (assume 2 spaces per level)
        size_t indent_level = indent / 2;

        // Trim whitespace from line
        std::string trimmed_line = line;
        trimmed_line.erase(trimmed_line.begin(), std::find_if(trimmed_line.begin(), trimmed_line.end(), [](unsigned char ch) {
            return !std::isspace(ch);
        }));
        trimmed_line.erase(std::find_if(trimmed_line.rbegin(), trimmed_line.rend(), [](unsigned char ch) {
            return !std::isspace(ch);
        }).base(), trimmed_line.end());

        if (trimmed_line.empty()) continue;

        // Use trimmed line for parsing
        line = trimmed_line;

        // If the line ends with a colon, it is a section header.
        if (line.back() == ':')
        {
            std::string section_name = line.substr(0, line.size() - 1);

            // Adjust section stack based on indentation
            while (section_stack.size() > indent_level) 
            {
                section_stack.pop_back();
            }

            // If the section stack size is equal to the indentation level, add the section name to the section stack.
            if (section_stack.size() == indent_level) {
                section_stack.push_back(section_name);
            } 
            else 
            {
                // Replace current level
                section_stack[indent_level] = section_name;
            }

            continue;
        }

        // Parse key-value pairs
        size_t colon_pos = line.find(':');

        // If the colon is not found, continue.
        if (colon_pos != std::string::npos)
        {
            std::string key = line.substr(0, colon_pos);
            std::string value = line.substr(colon_pos + 1);

            // Trim whitespace
            key.erase(key.begin(), std::find_if(key.begin(), key.end(), [](unsigned char ch) {
                return !std::isspace(ch);
            }));
            key.erase(std::find_if(key.rbegin(), key.rend(), [](unsigned char ch) {
                return !std::isspace(ch);
            }).base(), key.end());

            value.erase(value.begin(), std::find_if(value.begin(), value.end(), [](unsigned char ch) {
                return !std::isspace(ch);
            }));
            value.erase(std::find_if(value.rbegin(), value.rend(), [](unsigned char ch) {
                return !std::isspace(ch);
            }).base(), value.end());

            // Build full key from section stack
            std::string full_key;
            for (const auto& section : section_stack) 
            {
                if (!full_key.empty()) full_key += ".";
                full_key += section;
            }
            if (!full_key.empty()) full_key += ".";
            full_key += key;
            config[full_key] = value;
        }
    }

    return config;
}

// Global wind profile for initialization
WindProfile global_wind_profile = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

// Global CAPE target for thermodynamic profile scaling
double global_cape_target = 2500.0; // J/kg

// Global microphysics scheme name
std::string global_microphysics_scheme = "kessler";

// Global dynamics scheme name
std::string global_dynamics_scheme_name = "";

/*This function loads the configuration from a YAML file.
Takes in the configuration path, duration, write every, and output directory and loads the configuration from the YAML file.*/
void load_config(const std::string& config_path, int& duration_s, int& write_every_s, std::string& /*outdir*/)
{
    if (config_path.empty()) {
        return;
    }

    auto config = parse_yaml_simple(config_path);
    std::cout << "Loaded config with " << config.size() << " keys" << std::endl;

    //if the duration is specified in the config, load the duration.
    if (config.count("duration_s"))
    {
        duration_s = std::stoi(config["duration_s"]);
    }

    // if the write every is specified in the config, load the write every.
    if (config.count("output.write_every_s"))
    {
        write_every_s = std::stoi(config["output.write_every_s"]);
    }

    // if the CAPE target is specified in the config, load the CAPE target.
    if (config.count("environment.cape_target_jkg"))
    {
        global_cape_target = std::stod(config["environment.cape_target_jkg"]);
    }

    // If the nested grid is enabled in the config, load the nested grid.
    if (config.count("nested.enabled"))
    {
        nested_config.enabled = (config["nested.enabled"] == "true");
    }

    // If the nested grid refinement is specified in the config, load the nested grid refinement.
    if (config.count("nested.refinement"))
    {
        nested_config.refinement = std::stod(config["nested.refinement"]);
    }

    // If the nested grid size is specified in the config, load the nested grid size.
    if (config.count("nested.size_r"))
    {
        nested_config.nest_size_r = std::stoi(config["nested.size_r"]);
    }

    // If the nested grid size theta is specified in the config, load the nested grid size theta.
    if (config.count("nested.size_th"))
    {
        nested_config.nest_size_th = std::stoi(config["nested.size_th"]);
    }

    // If the nested grid size z is specified in the config, load the nested grid size z.
    if (config.count("nested.size_z"))
    {
        nested_config.nest_size_z = std::stoi(config["nested.size_z"]);
    }

    // Load grid parameters
    bool grid_changed = false;

    // If the grid size x is specified in the config, load the grid size x.
    if (config.count("grid.nx"))
    {
        int new_NR = std::stoi(config["grid.nx"]);

        // If the grid size x is not equal to the current grid size x, load the grid size x.
        if (new_NR != NR) 
        {
            NR = new_NR;
            grid_changed = true;
        }
    }

    // If the grid size y is specified in the config, load the grid size y.
    if (config.count("grid.ny"))
    {
        int new_NTH = std::stoi(config["grid.ny"]);

        // If the grid size y is not equal to the current grid size y, load the grid size y.
        if (new_NTH != NTH) 
        {
            NTH = new_NTH;
            grid_changed = true;
        }
    }

    // If the grid size z is specified in the config, load the grid size z.
    if (config.count("grid.nz"))
    {
        int new_NZ = std::stoi(config["grid.nz"]);

        // If the grid size z is not equal to the current grid size z, load the grid size z.
        if (new_NZ != NZ) 
        {
            NZ = new_NZ;
            grid_changed = true;
        }
    }

    // Load grid resolution parameters
    if (config.count("grid.dx"))
    {
        dr = std::stod(config["grid.dx"]);
    }

    // If the grid resolution y is specified in the config, load the grid resolution y.
    // Note: For cylindrical coordinates, dtheta is computed from NTH, but we can load dy for reference
    if (config.count("grid.dy"))
    {
        // Store for potential future use (currently dtheta is computed from NTH)
        // This parameter is loaded but may not be actively used in cylindrical coordinate system
    }

    // If the grid resolution z is specified in the config, load the grid resolution z.
    if (config.count("grid.dz"))
    {
        dz = std::stod(config["grid.dz"]);
    }

    // If the grid resolution t is specified in the config, load the grid resolution t.
    if (config.count("grid.dt"))
    {
        dt = std::stod(config["grid.dt"]);
    }

    // Always resize fields to ensure they're allocated (even if grid didn't change from defaults)
    {
        extern void resize_fields();
        resize_fields();
    }

    // If the hodograph u surface is specified in the config, load the hodograph u surface.
    if (config.count("environment.hodograph.u_sfc_ms"))
    {
        global_wind_profile.u_sfc = std::stod(config["environment.hodograph.u_sfc_ms"]);
    }

    // If the hodograph v surface is specified in the config, load the hodograph v surface.
    if (config.count("environment.hodograph.v_sfc_ms"))
    {
        global_wind_profile.v_sfc = std::stod(config["environment.hodograph.v_sfc_ms"]);
    }

    // If the hodograph u 1km is specified in the config, load the hodograph u 1km.
    if (config.count("environment.hodograph.u_1km_ms"))
    {
        global_wind_profile.u_1km = std::stod(config["environment.hodograph.u_1km_ms"]);
    }

    // If the hodograph v 1km is specified in the config, load the hodograph v 1km.
    if (config.count("environment.hodograph.v_1km_ms"))
    {
        global_wind_profile.v_1km = std::stod(config["environment.hodograph.v_1km_ms"]);
    }

    // If the hodograph u 6km is specified in the config, load the hodograph u 6km.
    if (config.count("environment.hodograph.u_6km_ms"))
    {
        global_wind_profile.u_6km = std::stod(config["environment.hodograph.u_6km_ms"]);
    }

    // If the hodograph v 6km is specified in the config, load the hodograph v 6km.
    if (config.count("environment.hodograph.v_6km_ms"))
    {
        global_wind_profile.v_6km = std::stod(config["environment.hodograph.v_6km_ms"]);
    }

    // Load microphysics scheme
    if (config.count("microphysics.scheme"))
    {
        std::string requested_scheme = config["microphysics.scheme"];
        // Validate scheme name
        std::vector<std::string> valid_schemes = {"kessler", "lin", "thompson", "milbrandt"};
        bool valid = false;

        // Iterate over the valid schemes and check if the requested scheme is valid.
        for (const auto& scheme : valid_schemes)
        {
            // If the requested scheme is equal to the current scheme, set the valid flag to true and break.
            if (requested_scheme == scheme) 
            {
                valid = true;
                break;
            }
        }

        // If the requested scheme is valid, set the microphysics scheme to the requested scheme.
        if (valid) 
        {
            global_microphysics_scheme = requested_scheme;
        } 
        else
        {
            std::cout << "Warning: Invalid microphysics scheme '" << requested_scheme << "'. Valid options: ";

            // Iterate over the valid schemes and print the valid schemes.
            for (size_t i = 0; i < valid_schemes.size(); ++i) 
            {
                std::cout << valid_schemes[i];
                if (i < valid_schemes.size() - 1) std::cout << ", ";
            }
            std::cout << ". Using default: kessler" << std::endl;
            global_microphysics_scheme = "kessler";
        }
    }

    // If the radiation scheme is specified in the config, load the radiation scheme.
    if (config.count("radiation.scheme"))
    {
        global_radiation_config.scheme_id = config["radiation.scheme"];
    }

    // If the radiation longwave is specified in the config, load the radiation longwave.
    if (config.count("radiation.do_lw"))
    {
        global_radiation_config.do_lw = (config["radiation.do_lw"] == "true");
    }

    // If the radiation shortwave is specified in the config, load the radiation shortwave.
    if (config.count("radiation.do_sw"))
    {
        global_radiation_config.do_sw = (config["radiation.do_sw"] == "true");
    }

    // If the radiation timestep is specified in the config, load the radiation timestep.
    if (config.count("radiation.dt"))
    {
        global_radiation_config.dt_radiation = std::stod(config["radiation.dt"]);
    }

    // If the radiation longwave reference is specified in the config, load the radiation longwave reference.
    if (config.count("radiation.tau_lw_ref"))
    {
        global_radiation_config.tau_lw_ref = std::stod(config["radiation.tau_lw_ref"]);
    }

    // If the radiation shortwave reference is specified in the config, load the radiation shortwave reference.
    if (config.count("radiation.tau_sw_ref"))
    {
        global_radiation_config.tau_sw_ref = std::stod(config["radiation.tau_sw_ref"]);
    }

    // If the boundary layer scheme is specified in the config, load the boundary layer scheme.
    if (config.count("boundary_layer.scheme"))
    {
        global_boundary_layer_config.scheme_id = config["boundary_layer.scheme"];
    }

    // If the boundary layer apply surface fluxes is specified in the config, load the boundary layer apply surface fluxes.
    if (config.count("boundary_layer.apply_surface_fluxes"))
    {
        global_boundary_layer_config.apply_surface_fluxes = (config["boundary_layer.apply_surface_fluxes"] == "true");
    }

    // If the boundary layer timestep is specified in the config, load the boundary layer timestep.
    if (config.count("boundary_layer.dt_pbl"))
    {
        global_boundary_layer_config.dt_pbl = std::stod(config["boundary_layer.dt_pbl"]);
    }

    // If the surface z0m is specified in the config, load the surface z0m.
    if (config.count("surface.z0m"))
    {
        global_surface_config.z0m = std::stod(config["surface.z0m"]);
    }

    // If the surface z0h is specified in the config, load the surface z0h.
    if (config.count("surface.z0h"))
    {
        global_surface_config.z0h = std::stod(config["surface.z0h"]);
    }

    // If the turbulence scheme is specified in the config, load the turbulence scheme.
    if (config.count("turbulence.scheme"))
    {
        global_turbulence_config.scheme_id = config["turbulence.scheme"];
    }

    // If the turbulence timestep is specified in the config, load the turbulence timestep.
    if (config.count("turbulence.dt_sgs"))
    {
        global_turbulence_config.dt_sgs = std::stod(config["turbulence.dt_sgs"]);
    }

    // If the turbulence Cs is specified in the config, load the turbulence Cs.
    if (config.count("turbulence.Cs"))
    {
        global_turbulence_config.Cs = std::stod(config["turbulence.Cs"]);
    }

    // If the turbulence Pr_t is specified in the config, load the turbulence Pr_t.
    if (config.count("turbulence.Pr_t"))
    {
        global_turbulence_config.Pr_t = std::stod(config["turbulence.Pr_t"]);
    }

    // If the dynamics scheme is specified in the config, load the dynamics scheme.
    if (config.count("dynamics.scheme"))
    {
        global_dynamics_scheme_name = config["dynamics.scheme"];
    }

    // Load numerics scheme configurations
    if (config.count("numerics.advection"))
    {
        global_advection_config.scheme_id = config["numerics.advection"];
    }
    if (config.count("numerics.diffusion"))
    {
        global_diffusion_config.scheme_id = config["numerics.diffusion"];
    }
    if (config.count("numerics.time_stepping"))
    {
        global_time_stepping_config.scheme_id = config["numerics.time_stepping"];
    }

    // Load chaos configuration parameters (will be used when initialize_chaos is called)
    if (config.count("chaos.scheme"))
    {
        global_chaos_config.scheme_id = config["chaos.scheme"];
    }
    if (config.count("chaos.seed"))
    {
        global_chaos_config.seed = std::stoull(config["chaos.seed"]);
    }
    if (config.count("chaos.member_id"))
    {
        global_chaos_config.member_id = std::stoi(config["chaos.member_id"]);
    }
    if (config.count("chaos.dt_chaos"))
    {
        global_chaos_config.dt_chaos = std::stod(config["chaos.dt_chaos"]);
    }
    if (config.count("chaos.tau_t"))
    {
        global_chaos_config.tau_t = std::stod(config["chaos.tau_t"]);
    }
    if (config.count("chaos.Lx"))
    {
        global_chaos_config.Lx = std::stod(config["chaos.Lx"]);
    }
    if (config.count("chaos.Ly"))
    {
        global_chaos_config.Ly = std::stod(config["chaos.Ly"]);
    }

    // Load terrain configuration parameters
    if (config.count("terrain.scheme"))
    {
        global_terrain_config.scheme_id = config["terrain.scheme"];
    }
    if (config.count("terrain.ztop"))
    {
        global_terrain_config.ztop = std::stod(config["terrain.ztop"]);
    }
    if (config.count("terrain.coord_id"))
    {
        global_terrain_config.coord_id = config["terrain.coord_id"];
    }

    // For now, we'll keep the grid parameters fixed
    // Future enhancement: make grid configurable
    std::cout << "\n" << "=" << std::string(60, '=') << std::endl;
    std::cout << "CONFIGURATION LOADED" << std::endl;
    std::cout << "=" << std::string(60, '=') << std::endl;
    std::cout << "Config file: " << config_path << std::endl;
    std::cout << "Duration: " << duration_s << "s" << std::endl;
    std::cout << "Write every: " << write_every_s << "s" << std::endl;
    std::cout << "\nGrid parameters:" << std::endl;
    std::cout << "  NR=" << NR << ", NTH=" << NTH << ", NZ=" << NZ << std::endl;
    std::cout << "  dr=" << dr << "m, dz=" << dz << "m, dt=" << dt << "s" << std::endl;
    std::cout << "\nEnvironment:" << std::endl;
    std::cout << "  CAPE target: " << global_cape_target << " J/kg" << std::endl;
    std::cout << "  Wind profile: SFC(" << global_wind_profile.u_sfc << "," << global_wind_profile.v_sfc
              << ") 1km(" << global_wind_profile.u_1km << "," << global_wind_profile.v_1km
              << ") 6km(" << global_wind_profile.u_6km << "," << global_wind_profile.v_6km << ")" << std::endl;
    std::cout << "\nPhysics schemes:" << std::endl;
    std::cout << "  Dynamics: " << global_dynamics_scheme_name << std::endl;
    std::cout << "  Microphysics: " << global_microphysics_scheme << std::endl;
    std::cout << "  Radiation: " << global_radiation_config.scheme_id
              << " (LW: " << (global_radiation_config.do_lw ? "on" : "off")
              << ", SW: " << (global_radiation_config.do_sw ? "on" : "off")
              << ", dt: " << global_radiation_config.dt_radiation << "s)" << std::endl;
    std::cout << "  Boundary Layer: " << global_boundary_layer_config.scheme_id
              << " (surface fluxes: " << (global_boundary_layer_config.apply_surface_fluxes ? "on" : "off")
              << ", dt: " << global_boundary_layer_config.dt_pbl << "s)" << std::endl;
    std::cout << "  Turbulence: " << global_turbulence_config.scheme_id
              << " (Cs: " << global_turbulence_config.Cs
              << ", dt: " << global_turbulence_config.dt_sgs << "s)" << std::endl;
    std::cout << "  Chaos: " << global_chaos_config.scheme_id << std::endl;
    std::cout << "=" << std::string(60, '=') << "\n" << std::endl;
}

/*This function is the main function.
Takes in the arguments and runs the simulation.*/
int main(int argc, char** argv) 
{
    // Basic CLI: --headless, --export-ms=N, --duration=SEC, --write-every=SEC, --outdir=PATH, --config=PATH, --export-all
    bool headless = false;
    int export_ms = 0;          // GUI auto-export cadence (ms)
    int duration_s = -1;        // headless run wall-clock duration (s), negative = infinite
    int write_every_s = 0;      // headless export cadence (s); if >0, write all theta slices periodically
    std::string outdir = "data/exports";
    std::string config_path = ""; // path to YAML config file
    // All fields are now always exported

    // Iterate over the arguments and parse the arguments.
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "--headless") headless = true;
        else if (arg.rfind("--export-ms=", 0) == 0)
        {
            export_ms = std::max(0, std::stoi(arg.substr(12)));
        }
        else if (arg.rfind("--duration=", 0) == 0)
        {
            duration_s = std::stoi(arg.substr(11));
        }
        else if (arg == "--duration" && i + 1 < argc)
        {
            duration_s = std::stoi(argv[++i]);
        }
        else if (arg.rfind("--write-every=", 0) == 0)
        {
            write_every_s = std::max(0, std::stoi(arg.substr(14)));
        }
        else if (arg.rfind("--outdir=", 0) == 0)
        {
            outdir = arg.substr(9);
        }
        else if (arg.rfind("--config=", 0) == 0)
        {
            config_path = arg.substr(9);
        }
        else if (arg == "--config" && i + 1 < argc)
        {
            config_path = argv[++i];
        }
        // --export-all flag removed - all fields are now always exported
    }

    // Load configuration from YAML file if specified
    load_config(config_path, duration_s, write_every_s, outdir);

    // Initialize microphysics scheme from config
    initialize_microphysics(global_microphysics_scheme);

    // Initialize radar scheme (reflectivity by default)
    initialize_radar("reflectivity");

    // Initialize the simulation state once (after grid is resized)
    initialize();

    // Initialize numerics framework (must be before chaos to set grid metrics)
    initialize_numerics();

    // Initialize chaos perturbation scheme (use loaded config if available)
    // Must be after numerics initialization to have correct grid metrics
    if (global_chaos_config.scheme_id.empty()) {
        global_chaos_config.scheme_id = "none";  // default
    }
    initialize_chaos(global_chaos_config);

    // Apply chaos initial condition perturbations (must be done after initialize() and chaos scheme initialization)
    apply_chaos_initial_conditions();

    // Initialize dynamics scheme (use config if available, default to tornado)
    std::string dynamics_scheme_name = "tornado";  // default
    if (!global_dynamics_scheme_name.empty()) {
        dynamics_scheme_name = global_dynamics_scheme_name;
    }
    initialize_dynamics(dynamics_scheme_name);

    // Initialize radiation scheme (use loaded config)
    if (global_radiation_config.scheme_id.empty()) {
        global_radiation_config.scheme_id = "simple_grey";  // default
    }
    initialize_radiation(global_radiation_config.scheme_id, global_radiation_config);

    // Initialize boundary layer scheme (use loaded config)
    if (global_boundary_layer_config.scheme_id.empty()) {
        global_boundary_layer_config.scheme_id = "slab";  // default
    }
    initialize_boundary_layer(global_boundary_layer_config.scheme_id, global_boundary_layer_config, global_surface_config);

    // Initialize turbulence scheme (use loaded config)
    if (global_turbulence_config.scheme_id.empty()) {
        global_turbulence_config.scheme_id = "smagorinsky";  // default
    }
    initialize_turbulence(global_turbulence_config.scheme_id, global_turbulence_config);

    // Terrain system is available but not initialized by default
    // To enable terrain, call: initialize_terrain("bell", TerrainConfig{}) or similar

    if (headless)
    {
        // Run a headless loop with optional periodic export
        const int thetaIndex = 0;
#ifdef EXPORT_NPY
        auto save_field_slice_npy = [&](const Field3D& field, int theta, const std::string& filename)
        {
            std::vector<float> buf;
            buf.resize(static_cast<size_t>(NR) * static_cast<size_t>(NZ));
            size_t idx = 0;
            
            // Debug: Check values before saving
            static bool debug_save = true;
            if (debug_save && theta == 0) {
                float field_min = 1e10, field_max = -1e10;
                int nan_count = 0;
                for (int i = 0; i < NR && i < 10; ++i) {
                    for (int k = 0; k < NZ && k < 10; ++k) {
                        float val = field[i][theta][k];
                        if (std::isnan(val)) nan_count++;
                        if (val < field_min) field_min = val;
                        if (val > field_max) field_max = val;
                    }
                }
                std::cout << "[SAVE DEBUG] Saving field slice theta=" << theta << ":" << std::endl;
                std::cout << "  Sample range: [" << field_min << ", " << field_max << "]" << std::endl;
                std::cout << "  NaN count (sample): " << nan_count << std::endl;
                debug_save = false;
            }
            
            for (int k = 0; k < NZ; ++k)
            {
                for (int i = 0; i < NR; ++i)
                {
                    float val = static_cast<float>(field[i][theta][k]);
                    buf[idx++] = val;
                }
            }
            std::string header_dict = "{'descr': '<f4', 'fortran_order': False, 'shape': (" + std::to_string(NZ) + ", " + std::to_string(NR) + "), }";
            size_t header_len = header_dict.size() + 1;
            const size_t preamble = 6 + 2 + 2;
            size_t total = preamble + header_len;
            size_t padding = (16 - (total % 16)) % 16;
            header_len += padding;
            std::ofstream out(filename, std::ios::binary);
            if (out)
            {
                out.write("\x93NUMPY", 6);
                out.put(static_cast<char>(1));
                out.put(static_cast<char>(0));
                uint16_t hl = static_cast<uint16_t>(header_len);
                char lenb[2]; lenb[0] = static_cast<char>(hl & 0xFF); lenb[1] = static_cast<char>((hl >> 8) & 0xFF);
                out.write(lenb, 2);
                out.write(header_dict.c_str(), static_cast<std::streamsize>(header_dict.size()));
                for (size_t i = 0; i < header_len - (header_dict.size() + 1); ++i) out.put(' ');
                out.put('\n');
                out.write(reinterpret_cast<const char*>(buf.data()), static_cast<std::streamsize>(buf.size() * sizeof(float)));
                out.close();
            }
        };

        /*This function writes all the fields to the output directory.
        Takes in the export index and writes all the fields to the output directory.*/
        auto write_all_fields = [&](int export_index)
        {
            // Debug: Check values before writing
            if (export_index == 0) {
                std::cout << "\n[EXPORT DEBUG] Writing timestep " << export_index << " (t=" << simulation_time << "s)" << std::endl;
                float theta_min = 1e10, theta_max = -1e10;
                float u_min = 1e10, u_max = -1e10;
                int nan_count = 0;
                for (int i = 0; i < NR && i < 10; ++i) {
                    for (int j = 0; j < NTH && j < 3; ++j) {
                        for (int k = 0; k < NZ && k < 3; ++k) {
                            if (std::isnan(theta[i][j][k])) nan_count++;
                            if (theta[i][j][k] < theta_min) theta_min = theta[i][j][k];
                            if (theta[i][j][k] > theta_max) theta_max = theta[i][j][k];
                            if (u[i][j][k] < u_min) u_min = u[i][j][k];
                            if (u[i][j][k] > u_max) u_max = u[i][j][k];
                        }
                    }
                }
                std::cout << "  Theta sample: [" << theta_min << ", " << theta_max << "] K" << std::endl;
                std::cout << "  Wind (u) sample: [" << u_min << ", " << u_max << "] m/s" << std::endl;
                std::cout << "  NaN count (sample): " << nan_count << std::endl;
                if (theta_min < 0 || theta_max > 500) {
                    std::cerr << "  ⚠️  ERROR: Theta values are wrong before export!" << std::endl;
                }
                std::cout << std::endl;
            }
            
            // Create outdir/step_XXXXXX
            std::ostringstream stepdir;
            stepdir << outdir << "/step_" << std::setfill('0') << std::setw(6) << export_index;
            std::filesystem::create_directories(stepdir.str());

            for (int th = 0; th < NTH; ++th)
            {
                std::string base_path = stepdir.str() + "/th" + std::to_string(th);

                // Export all fields
                save_field_slice_npy(u, th, base_path + "_u.npy");
                save_field_slice_npy(v_theta, th, base_path + "_v.npy");
                save_field_slice_npy(w, th, base_path + "_w.npy");
                save_field_slice_npy(rho, th, base_path + "_rho.npy");
                save_field_slice_npy(p, th, base_path + "_p.npy");
                save_field_slice_npy(theta, th, base_path + "_theta.npy");
                save_field_slice_npy(qv, th, base_path + "_qv.npy");
                save_field_slice_npy(qc, th, base_path + "_qc.npy");
                save_field_slice_npy(qr, th, base_path + "_qr.npy");
                save_field_slice_npy(qh, th, base_path + "_qh.npy");
                save_field_slice_npy(qg, th, base_path + "_qg.npy");
                save_field_slice_npy(radar_reflectivity, th, base_path + "_radar.npy");
                save_field_slice_npy(tracer, th, base_path + "_tracer.npy");
            }
        };

#endif

        // Ensure output directory exists if periodic writing is requested
        if (write_every_s > 0)
        {
            std::filesystem::create_directories(outdir);
        }

        auto lastExport = std::chrono::steady_clock::now();
        auto startRun = std::chrono::steady_clock::now();
        int steps = 0;
        int export_index = 0;
        double simulation_time = 0.0;  // simulation time in seconds
        
        // Debug: Check values before time stepping
        {
            float theta_min = 1e10, theta_max = -1e10;
            float u_min = 1e10, u_max = -1e10;
            int nan_count = 0;
            for (int i = 0; i < NR && i < 10; ++i) {
                for (int j = 0; j < NTH && j < 5; ++j) {
                    for (int k = 0; k < NZ && k < 5; ++k) {
                        if (std::isnan(theta[i][j][k])) nan_count++;
                        if (theta[i][j][k] < theta_min) theta_min = theta[i][j][k];
                        if (theta[i][j][k] > theta_max) theta_max = theta[i][j][k];
                        if (u[i][j][k] < u_min) u_min = u[i][j][k];
                        if (u[i][j][k] > u_max) u_max = u[i][j][k];
                    }
                }
            }
            std::cout << "\n[TIME STEP DEBUG] Before time stepping (t=0):" << std::endl;
            std::cout << "  Theta sample: min=" << theta_min << "K, max=" << theta_max << "K" << std::endl;
            std::cout << "  Wind (u) sample: min=" << u_min << "m/s, max=" << u_max << "m/s" << std::endl;
            std::cout << "  NaN count (sample): " << nan_count << std::endl;
            if (theta_min < 0 || theta_max > 500) {
                std::cerr << "  ⚠️  ERROR: Theta already corrupted before time stepping!" << std::endl;
            }
            std::cout << std::endl;
        }
        
        while (true)
        {
            // Apply physics processes
            step_radiation(simulation_time);
            step_boundary_layer(simulation_time);
            step_dynamics(simulation_time);
            simulation_time += dt;  // increment simulation time
            
            // Debug: Check values periodically during time stepping
            if (steps % 100 == 0 && steps < 1000) {
                float theta_min = 1e10, theta_max = -1e10;
                float u_min = 1e10, u_max = -1e10;
                int nan_count = 0, inf_count = 0;
                for (int i = 0; i < NR && i < 10; ++i) {
                    for (int j = 0; j < NTH && j < 5; ++j) {
                        for (int k = 0; k < NZ && k < 5; ++k) {
                            if (std::isnan(theta[i][j][k])) nan_count++;
                            if (std::isinf(theta[i][j][k])) inf_count++;
                            if (theta[i][j][k] < theta_min) theta_min = theta[i][j][k];
                            if (theta[i][j][k] > theta_max) theta_max = theta[i][j][k];
                            if (u[i][j][k] < u_min) u_min = u[i][j][k];
                            if (u[i][j][k] > u_max) u_max = u[i][j][k];
                        }
                    }
                }
                std::cout << "[TIME STEP DEBUG] Step " << steps << " (t=" << simulation_time << "s):" << std::endl;
                std::cout << "  Theta sample: min=" << theta_min << "K, max=" << theta_max << "K" << std::endl;
                std::cout << "  Wind (u) sample: min=" << u_min << "m/s, max=" << u_max << "m/s" << std::endl;
                std::cout << "  NaN/Inf count (sample): " << nan_count << "/" << inf_count << std::endl;
                if (theta_min < 0 || theta_max > 500 || std::abs(u_min) > 1000 || std::abs(u_max) > 1000) {
                    std::cerr << "  ⚠️  WARNING: Values going wrong at step " << steps << "!" << std::endl;
                }
                std::cout << std::endl;
            }
#ifdef EXPORT_NPY
            // GUI-style ms export (kept for compatibility when headless + --export-ms)
            if (export_ms > 0)
            {
                auto now = std::chrono::steady_clock::now();
                if (std::chrono::duration_cast<std::chrono::milliseconds>(now - lastExport).count() >= export_ms)
                {
                    lastExport = now;
                    // Preserve previous behavior: write a single theta slice (thetaIndex=0) to data/
                    std::string tmp = std::string("data/.tracer_slice_th0.npy.tmp");
                    std::string fin = std::string("data/tracer_slice_th0.npy");
                    save_field_slice_npy(tracer, thetaIndex, tmp);
                    std::rename(tmp.c_str(), fin.c_str());
                }
            }
            // New: periodic full-theta export in seconds to outdir/step_xxxxxx
            if (write_every_s > 0)
            {
                auto now = std::chrono::steady_clock::now();
                if (std::chrono::duration_cast<std::chrono::seconds>(now - lastExport).count() >= write_every_s)
                {
                    lastExport = now;
                    write_all_fields(export_index++);
                }
            }
#endif
            ++steps;
            if (steps % 1000 == 0) { /* simple heartbeat */ }
            if (duration_s >= 0)
            {
                auto now = std::chrono::steady_clock::now();
                if (std::chrono::duration_cast<std::chrono::seconds>(now - startRun).count() >= duration_s)
                {
                    break;
                }
            }
        }
    }
    else
    {
        // No GUI in default build; if GUI is enabled, delegate.
        #ifdef ENABLE_GUI
        run_gui(export_ms);
        #else
        (void)export_ms; // unused
        #endif
    }

    return 0;
}


