#pragma once
#include "terrain_base.hpp"

/*This file contains the declaration of the topography module.
It manages the evaluation of the topography and the initialization of the topography.
Commented out code is for the future and was having issues with the compilation. with 
binary operator issues.*/

// Topography generation and coordinate utilities
namespace topography {

// Evaluate bell-shaped mountain/ridge
struct BellResult {
    double h;   // height [m]
    double hx;  // dh/dx [m/m]
    double hy;  // dh/dy [m/m]
};

BellResult eval_bell(double x, double y, const TerrainConfig::BellParams& params);

// Evaluate Schär mountain (MWR 2002 test case)
struct ScharResult {
    double h;   // height [m]
    double hx;  // dh/dx [m/m]
};

ScharResult eval_schar(double x, const TerrainConfig::ScharParams& params);

// Terrain-following coordinate transformations
// TODO: Implement terrain coordinate transformations - commented out due to compilation issues
/*
class TerrainFollowingCoordinate {
private:
    TerrainConfig config_;
    Topography2D topography_;

public:
    TerrainFollowingCoordinate(const TerrainConfig& cfg, const Topography2D& topo);

    // Transform from terrain-following ζ to physical z
    double zeta_to_z(double zeta, int i, int j) const;

    // Transform from physical z to terrain-following ζ
    double z_to_zeta(double z, int i, int j) const;

    // Compute metric terms at a grid point
    void compute_metrics(double zeta, int i, int j,
                        double& z, double& J, double& mx, double& my) const;

    // Smoothed coordinate surfaces (optional)
    double apply_coordinate_smoothing(double h, double zeta) const;
};
*/

// Initialize 2D topography field
void initialize_topography(Topography2D& topo, int NR, int NTH);

// Initialize 3D metrics field
void initialize_metrics(TerrainMetrics3D& metrics, int NR, int NTH, int NZ);

// Build vertical coordinate levels (ζ)
std::vector<double> build_zeta_levels(int NZ, double ztop);

// Check coordinate validity (no folding)
bool check_coordinate_validity(const TerrainMetrics3D& metrics, TerrainDiagnostics& diag);
