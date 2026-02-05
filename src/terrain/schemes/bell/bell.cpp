#include "bell.hpp"
#include "simulation.hpp"
#include <iostream>
#include <cmath>

/*This file contains the implementation of the bell terrain scheme.
It manages the initialization of the bell terrain scheme and the building of the topography.*/

BellScheme::BellScheme() 
{
}

/*This function initializes the bell terrain scheme.
Takes in the configuration and initializes the bell terrain scheme.*/
void BellScheme::initialize(const TerrainConfig& cfg) 
{
    config_ = cfg;

    std::cout << "Initialized Bell terrain scheme" << std::endl;
    std::cout << "  h0 = " << cfg.bell.h0 << " m" << std::endl;
    std::cout << "  a = " << cfg.bell.a << " m" << std::endl;
    std::cout << "  axisymmetric = " << (cfg.bell.axisymmetric ? "true" : "false") << std::endl;
    std::cout << "  coordinate type = " << cfg.coord_id << std::endl;
}

/*This function builds the topography.
Takes in the configuration and the topography and builds the topography.*/
void BellScheme::build_topography(const TerrainConfig& cfg, Topography2D& topo) 
{
    const int NR = topo.h.size();
    const int NTH = NR > 0 ? topo.h[0].size() : 0;

    // Center the bell in the domain
    const double center_i = (NR - 1) / 2.0;
    const double center_j = (NTH - 1) / 2.0;


    // Iterate over the rows and columns and build the topography.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over the columns and build the topography.
        for (int j = 0; j < NTH; ++j) 
        {
            double x, y;
            get_coordinates(i, j, x, y);

            // Center the bell
            double x_centered = x - center_i * dr;
            double y_centered = y - center_j * dtheta * 1000.0;  // approximate

            auto result = topography::eval_bell(x_centered, y_centered, cfg.bell);

            topo.h[i][j] = result.h;

            // If the derivatives are computed, add the derivatives to the topography.
            if (cfg.compute_derivatives) 
            {
                topo.hx[i][j] = result.hx;
                topo.hy[i][j] = result.hy;
            }
        }
    }
}

/*This function builds the metrics.
Takes in the configuration and the topography and the metrics and builds the metrics.*/
void BellScheme::build_metrics(const TerrainConfig& cfg,
                              const Topography2D& topo,
                              TerrainMetrics3D& metrics,
                              TerrainDiagnostics* diag_opt) {
    const int NR = metrics.z.size();
    const int NTH = NR > 0 ? metrics.z[0].size() : 0;
    const int NZ = NTH > 0 ? metrics.z[0][0].size() : 0;

    // Build zeta levels
    auto zeta_levels = topography::build_zeta_levels(NZ, cfg.ztop);

    // Create coordinate transformation
    topography::TerrainFollowingCoordinate coord(cfg, topo);

    // Iterate over the rows, columns, and levels and compute the metrics.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over the columns and compute the metrics.
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over the levels and compute the metrics.
            for (int k = 0; k < NZ; ++k) 
            {
                double zeta = zeta_levels[k];
                double z, J, mx, my;

                coord.compute_metrics(zeta, i, j, z, J, mx, my);

                metrics.z[i][j][k] = z;
                metrics.J[i][j][k] = J;
                metrics.mx[i][j][k] = mx;
                metrics.my[i][j][k] = my;
                metrics.zeta[i][j][k] = zeta;
            }
        }
    }

    // If the diagnostics are requested, fill the diagnostics.
    if (diag_opt) 
    {
        diag_opt->max_height = cfg.bell.h0;
        topography::check_coordinate_validity(metrics, *diag_opt);
    }
}

/*This function gets the coordinates.
Takes in the row and column and the coordinates and gets the coordinates.*/
void BellScheme::get_coordinates(int i, int j, double& x, double& y) const 
{
    // Convert grid indices to physical coordinates
    // Center at (0,0) for simplicity
    x = (i - NR/2) * dr;
    y = (j - NTH/2) * dtheta * 1000.0;  // approximate arc length
}
