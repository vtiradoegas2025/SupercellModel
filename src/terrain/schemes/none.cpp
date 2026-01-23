#include "none.hpp"
#include <iostream>

/*This file contains the implementation of the none terrain scheme.
It manages the initialization of the none terrain scheme and the building of the topography.*/

NoneScheme::NoneScheme() {}


/*This function initializes the none terrain scheme.
Takes in the configuration and initializes the none terrain scheme.*/

void NoneScheme::initialize(const TerrainConfig& cfg) 
{
    std::cout << "Initialized None terrain scheme (flat terrain)" << std::endl;
}


/*This function builds the topography.
Takes in the configuration and the topography and builds the topography.*/
void NoneScheme::build_topography(const TerrainConfig& cfg, Topography2D& topo) 
{
    const int NR = topo.h.size();
    const int NTH = NR > 0 ? topo.h[0].size() : 0;

    // Iterate over the rows and columns and set the heights to zero.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over the columns and set the heights to zero.
        for (int j = 0; j < NTH; ++j) 
        {
            // Set the heights to zero.
            topo.h[i][j] = 0.0;
            if (!topo.hx.empty()) topo.hx[i][j] = 0.0;
            if (!topo.hy.empty()) topo.hy[i][j] = 0.0;
        }
    }
}

/*This function builds the metrics.
Takes in the configuration and the topography and the metrics and builds the metrics.*/
void NoneScheme::build_metrics(const TerrainConfig& cfg,
                              const Topography2D& topo,
                              TerrainMetrics3D& metrics,
                              TerrainDiagnostics* diag_opt) 
{
    const int NR = metrics.z.size();
    const int NTH = NR > 0 ? metrics.z[0].size() : 0;
    const int NZ = NTH > 0 ? metrics.z[0][0].size() : 0;

    // Build zeta levels
    auto zeta_levels = topography::build_zeta_levels(NZ, cfg.ztop);

    // Cartesian coordinates (z = ζ, J = 1, m_x = m_y = 0)
    // Iterate over the rows, columns, and levels and build the metrics.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over the columns and build the metrics.
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over the levels and build the metrics.
            for (int k = 0; k < NZ; ++k) {
                double zeta = zeta_levels[k];
                metrics.z[i][j][k] = zeta;
                metrics.J[i][j][k] = 1.0;
                metrics.mx[i][j][k] = 0.0;
                metrics.my[i][j][k] = 0.0;
                metrics.zeta[i][j][k] = zeta;
            }
        }
    }

    // If the diagnostics are requested, fill the diagnostics.
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
