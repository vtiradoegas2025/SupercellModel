#include "simulation.hpp"
#include "turbulence/factory.hpp"
#include "turbulence/base/eddy_viscosity.hpp"
#include <iostream>
#include <vector>

// Global turbulence scheme instance and configuration
std::unique_ptr<TurbulenceSchemeBase> turbulence_scheme = nullptr;
TurbulenceConfig global_turbulence_config;

// Turbulence tendency fields
TurbulenceTendencies turbulence_tendencies;

// Initialize turbulence scheme
void initialize_turbulence(const std::string& scheme_name,
                          const TurbulenceConfig& cfg) {
    try {
        global_turbulence_config = cfg;
        turbulence_scheme = create_turbulence_scheme(scheme_name);
        turbulence_scheme->initialize(global_turbulence_config);

        // Initialize tendency fields (will be resized when needed)
        std::cout << "Initialized turbulence scheme: " << scheme_name << std::endl;
        std::cout << "  SGS cadence: " << global_turbulence_config.dt_sgs << " s" << std::endl;
        std::cout << "  Cs = " << global_turbulence_config.Cs << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error initializing turbulence: " << e.what() << std::endl;
        throw;
    }
}

// Step turbulence forward in time
void step_turbulence(double current_time, TurbulenceTendencies& tendencies) {
    if (!turbulence_scheme) return;

    // Only compute turbulence at specified cadence
    static double last_turbulence_time = -global_turbulence_config.dt_sgs;

    if (current_time - last_turbulence_time < global_turbulence_config.dt_sgs) {
        // Hold previous tendencies
        return;
    }

    last_turbulence_time = current_time;

    // Set up state view
    TurbulenceStateView state;
    state.u = &u;
    state.v = &v_theta;
    state.w = &w;
    state.rho = &rho;
    state.theta = &theta;
    state.qv = &qv;  // may be empty
    state.tke = &tke;  // may be empty for Smagorinsky
    state.NR = NR;
    state.NTH = NTH;
    state.NZ = NZ;

    // Set up grid metrics
    GridMetrics grid;
    grid.dx = dr;
    grid.dy = dtheta * 1000.0;  // approximate arc length in meters
    grid.dz.resize(NZ);
    for (int k = 0; k < NZ; ++k) {
        grid.dz[k] = dz;
    }
    grid.z_int.resize(NZ + 1);
    grid.z_int[0] = 0.0;
    for (int k = 1; k <= NZ; ++k) {
        grid.z_int[k] = grid.z_int[k-1] + dz;
    }

    // Resize tendency arrays
    eddy_viscosity::initialize_3d_field(tendencies.dudt_sgs, NR, NTH, NZ, 0.0f);
    eddy_viscosity::initialize_3d_field(tendencies.dvdt_sgs, NR, NTH, NZ, 0.0f);
    eddy_viscosity::initialize_3d_field(tendencies.dwdt_sgs, NR, NTH, NZ, 0.0f);
    eddy_viscosity::initialize_3d_field(tendencies.dthetadt_sgs, NR, NTH, NZ, 0.0f);
    eddy_viscosity::initialize_3d_field(tendencies.dqvdt_sgs, NR, NTH, NZ, 0.0f);
    eddy_viscosity::initialize_3d_field(tendencies.dtkedt_sgs, NR, NTH, NZ, 0.0f);

    // Compute turbulence
    turbulence_scheme->compute(global_turbulence_config, grid, state, tendencies);
}
