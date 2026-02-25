/**
 * @file none.cpp
 * @brief Implementation for the terrain module.
 *
 * Provides executable logic for the terrain runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/terrain subsystem.
 */

#include "none.hpp"
#include <iostream>


NoneScheme::NoneScheme() {}


/**
 * @brief Initializes the none terrain scheme.
 */

void NoneScheme::initialize(const TerrainConfig& cfg) 
{
    std::cout << "Initialized None terrain scheme (flat terrain)" << std::endl;
}


/**
 * @brief Builds the topography.
 */
void NoneScheme::build_topography(const TerrainConfig& cfg, Topography2D& topo) 
{
    const int NR = topo.h.size();
    const int NTH = NR > 0 ? topo.h[0].size() : 0;

    for (int i = 0; i < NR; ++i) 
    {
        for (int j = 0; j < NTH; ++j) 
        {
            topo.h[i][j] = 0.0;
            if (!topo.hx.empty()) topo.hx[i][j] = 0.0;
            if (!topo.hy.empty()) topo.hy[i][j] = 0.0;
        }
    }
}

/**
 * @brief Builds the metrics.
 */
void NoneScheme::build_metrics(const TerrainConfig& cfg,
                              const Topography2D& topo,
                              TerrainMetrics3D& metrics,
                              TerrainDiagnostics* diag_opt) 
{
    const int NR = metrics.z.size_r();
    const int NTH = metrics.z.size_th();
    const int NZ = metrics.z.size_z();

    auto zeta_levels = topography::build_zeta_levels(NZ, cfg.ztop);

    for (int i = 0; i < NR; ++i) 
    {
        for (int j = 0; j < NTH; ++j) 
        {
            for (int k = 0; k < NZ; ++k) {
                double zeta = zeta_levels[k];
                metrics.z(i, j, k) = zeta;
                metrics.J(i, j, k) = 1.0;
                metrics.mx(i, j, k) = 0.0;
                metrics.my(i, j, k) = 0.0;
                metrics.zeta(i, j, k) = zeta;
            }
        }
    }

    if (diag_opt) 
    {
        diag_opt->max_height = 0.0;
        diag_opt->max_slope_x = 0.0;
        diag_opt->max_slope_y = 0.0;
        diag_opt->min_jacobian = 1.0;
        diag_opt->max_jacobian = 1.0;
        diag_opt->coordinate_folding = false;
    }
}
