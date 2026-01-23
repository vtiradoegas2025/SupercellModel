#include "simulation.hpp"
#include "terrain/factory.hpp"
#include <iostream>

/*This file contains the implementation of the terrain system.
It manages the initialization of the terrain system and the computation of the terrain system.*/

// Global terrain scheme instance and configuration
std::unique_ptr<TerrainSchemeBase> terrain_scheme;
TerrainConfig global_terrain_config;
Topography2D global_topography;
TerrainMetrics3D global_terrain_metrics;

/*This function initializes the terrain scheme.
Takes in the scheme name and configuration and initializes the terrain scheme.*/
void initialize_terrain(const std::string& scheme_name, const TerrainConfig& cfg) {
    try {
        global_terrain_config = cfg;
        // Create a simple "none" terrain scheme that does nothing
        terrain_scheme = create_terrain_scheme("none");
        terrain_scheme->initialize(global_terrain_config);

        // Initialize empty topography and metrics fields
        build_terrain_fields();

        std::cout << "Initialized terrain scheme: " << scheme_name << " (stub implementation)" << std::endl;

    } 
    catch (const std::exception& e) 
    {
        std::cerr << "Error initializing terrain: " << e.what() << std::endl;
        throw;
    }
}

/*This function builds the terrain fields.
Takes in the topography and metrics and builds the terrain fields.*/
void build_terrain_fields() 
{
    // Initialize topography with zeros
    global_topography.h.assign(NR, std::vector<double>(NTH, 0.0));
    global_topography.hx.assign(NR, std::vector<double>(NTH, 0.0));
    global_topography.hy.assign(NR, std::vector<double>(NTH, 0.0));

    // Initialize metrics with identity transformations
    global_terrain_metrics.z.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));
    global_terrain_metrics.J.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 1.0)));
    global_terrain_metrics.mx.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));
    global_terrain_metrics.my.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));
    global_terrain_metrics.zeta.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));

    // Set zeta levels (uniform spacing)
    double dz_zeta = global_terrain_config.ztop / (NZ - 1);

    // Iterate over the levels and set the zeta levels.
    for (int k = 0; k < NZ; ++k) 
    {
        double zeta = k * dz_zeta;
        // Iterate over the rows and columns and set the zeta levels.
        for (int i = 0; i < NR; ++i) 
        {
            // Iterate over the columns and set the zeta levels.
            for (int j = 0; j < NTH; ++j) 
            {
                global_terrain_metrics.zeta[i][j][k] = zeta;
                global_terrain_metrics.z[i][j][k] = zeta;  // z = zeta for flat terrain
            }
        }
    }

    std::cout << "Terrain fields initialized (flat terrain)" << std::endl;
}
