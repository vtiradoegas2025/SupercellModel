#include "simulation.hpp"
#include "numerics/advection/factory.hpp"
#include "numerics/diffusion/factory.hpp"
#include "numerics/time_stepping/factory.hpp"
#include <iostream>

// Define the global numerics variables declared in simulation.hpp
std::unique_ptr<AdvectionSchemeBase> advection_scheme;
AdvectionConfig global_advection_config;
std::unique_ptr<DiffusionSchemeBase> diffusion_scheme;
DiffusionConfig global_diffusion_config;
std::unique_ptr<TimeSteppingSchemeBase> time_stepping_scheme;
TimeSteppingConfig global_time_stepping_config;

// Define global grid metrics declared in turbulence_base.hpp
GridMetrics global_grid_metrics;

// Initialize numerics schemes
void initialize_numerics() {
    // Initialize grid metrics (these would be set from the actual grid)
    global_grid_metrics.dx = dr;
    global_grid_metrics.dy = dtheta * 1000.0;  // approximate arc length
    global_grid_metrics.dz.assign(NZ, dz);
    global_grid_metrics.z_int.resize(NZ + 1);
    global_grid_metrics.z_int[0] = 0.0;
    for (int k = 1; k <= NZ; ++k) {
        global_grid_metrics.z_int[k] = global_grid_metrics.z_int[k-1] + dz;
    }

    // Initialize advection scheme (use config if set, default TVD)
    if (global_advection_config.scheme_id.empty()) {
        global_advection_config.scheme_id = "tvd";
    }
    advection_scheme = create_advection_scheme(global_advection_config.scheme_id);
    advection_scheme->initialize(global_advection_config);

    // Initialize diffusion scheme (use config if set, default explicit)
    if (global_diffusion_config.scheme_id.empty()) {
        global_diffusion_config.scheme_id = "explicit";
    }
    diffusion_scheme = create_diffusion_scheme(global_diffusion_config.scheme_id);
    diffusion_scheme->initialize(global_diffusion_config);

    // Initialize time stepping scheme (use config if set, default RK3)
    if (global_time_stepping_config.scheme_id.empty()) {
        global_time_stepping_config.scheme_id = "rk3";
    }
    time_stepping_scheme = create_time_stepping_scheme(global_time_stepping_config.scheme_id);

    std::cout << "Initialized numerics framework:" << std::endl;
    std::cout << "  Advection: " << advection_scheme->name() << std::endl;
    std::cout << "  Diffusion: " << diffusion_scheme->name() << std::endl;
    std::cout << "  Time stepping: " << time_stepping_scheme->name() << std::endl;
}
