#include "topography.hpp"
#include <cmath>
#include <algorithm>

/*This file contains the implementation of the topography module.
It manages the evaluation of the topography and the initialization of the topography.
Commented out code is for the future and was having issues with the compilation. with 
binary operator issues.*/

namespace topography 
{

/*This function evaluates the bell topography.
Takes in the x and y coordinates and the parameters and evaluates the bell topography.*/
BellResult eval_bell(double x, double y, const TerrainConfig::BellParams& params) 
{
    BellResult result;

    // If the topography is axisymmetric, evaluate the axisymmetric bell topography.
    if (params.axisymmetric) 
    {
        // 3D axisymmetric bell: h(r) = h0 * a^3 / (a^2 + r^2)^{3/2}
        double r2 = x*x + y*y;
        double denom = params.a*params.a + r2;
        double denom_32 = std::pow(denom, 1.5);

        result.h = params.h0 * params.a*params.a*params.a / denom_32;

        // If the denominator is not zero, evaluate the derivatives.
        if (denom_32 > 1e-12) 
        {
            double factor = -3.0 * params.h0 * params.a*params.a*params.a / std::pow(denom, 2.5);
            result.hx = factor * x;
            result.hy = factor * y;
        }
    } 
    else 
    {
        // 2D ridge: h(x) = h0 / (1 + (x/a)^2)
        double xa_ratio = x / params.a;
        double denom = 1.0 + xa_ratio*xa_ratio;

        result.h = params.h0 / denom;

        // Derivative dh/dx
        result.hx = -2.0 * params.h0 * xa_ratio / (params.a * denom*denom);
        result.hy = 0.0;  // 2D ridge independent of y
    }

    return result;
}

/*This function evaluates the Schär topography.
Takes in the x coordinate and the parameters and evaluates the Schär topography.*/
ScharResult eval_schar(double x, const TerrainConfig::ScharParams& params) 
{
    ScharResult result;

    // Schär mountain: h(x) = H * exp(-(x/a)^2) * cos^2(π*x/ℓ)
    double xa_ratio = x / params.a;
    double envelope = std::exp(-xa_ratio*xa_ratio);
    double k = M_PI / params.ell;
    double cos_term = std::cos(k * x);
    double cos2_term = cos_term * cos_term;

    result.h = params.H * envelope * cos2_term;

    // Derivative dh/dx
    double d_envelope_dx = -2.0 * xa_ratio / params.a * envelope;
    double d_cos2_dx = -2.0 * k * cos_term * std::sin(k * x);

    result.hx = params.H * (d_envelope_dx * cos2_term + envelope * d_cos2_dx);

    return result;
}

// TerrainFollowingCoordinate methods
TerrainFollowingCoordinate::TerrainFollowingCoordinate(const TerrainConfig& cfg, const Topography2D& topo)
    : config_(cfg), topography_(topo) 
{
}

double TerrainFollowingCoordinate::zeta_to_z(double zeta, int i, int j) const 
{
    double h = topography_.h[i][j];
    double ztop = config_.ztop;

    // If the coordinate id is "btf", evaluate the basic terrain-following.
    if (config_.coord_id == "btf") 
    {
        // Basic terrain-following: z = ζ + h * (1 - ζ/ztop)
        return zeta + h * (1.0 - zeta / ztop);
    }
    // If the coordinate id is "smoothed", evaluate the smoothed coordinate surfaces.
    else if (config_.coord_id == "smoothed") 
    {
        // Smoothed coordinate surfaces
        double h_eff = apply_coordinate_smoothing(h, zeta);
        return zeta + h_eff * (1.0 - zeta / ztop);
    } 
    else 
    {
        // Default to Cartesian
        return zeta;
    }
}

/*This function converts the z coordinate to the zeta coordinate.*/
double TerrainFollowingCoordinate::z_to_zeta(double z, int i, int j) const 
{
    double h = topography_.h[i][j];
    double ztop = config_.ztop;

    if (config_.coord_id == "btf") {
        // Basic terrain-following: ζ = ztop * (z - h) / (ztop - h)
        if (std::abs(ztop - h) > terrain_constants::epsilon) {
            return ztop * (z - h) / (ztop - h);
        } else {
            return z;  // degenerate case
        }
    } else if (config_.coord_id == "smoothed") {
        // Smoothed coordinate surfaces (inverse mapping)
        // This is approximate - full inverse would need iterative solution
        double h_eff_avg = h * 0.5;  // simplified
        if (std::abs(ztop - h_eff_avg) > terrain_constants::epsilon) {
            return ztop * (z - h_eff_avg) / (ztop - h_eff_avg);
        } else {
            return z;
        }
    } 
    else {
        return z;  // Cartesian
    }
}

void TerrainFollowingCoordinate::compute_metrics(double zeta, int i, int j,
                                                double& z_out, double& J, double& mx, double& my) const 
{
    double h = topography_.h[i][j];
    double hx = topography_.hx.empty() ? 0.0 : topography_.hx[i][j];
    double hy = topography_.hy.empty() ? 0.0 : topography_.hy[i][j];
    double ztop = config_.ztop;

    // If the coordinate id is "btf", evaluate the basic terrain-following metrics. 
    if (config_.coord_id == "btf")
    {
        // Basic terrain-following metrics
        z_out = zeta + h * (1.0 - zeta/ztop);

        // Jacobian dz/dζ
        J = 1.0 - h/ztop;

        // Metric terms
        double zeta_factor = 1.0 - zeta/ztop;

        // If the Jacobian is not zero, evaluate the metric terms.
        if (std::abs(J) > terrain_constants::epsilon) 
        {
            mx = -zeta_factor * hx / J;
            my = -zeta_factor * hy / J;
        }
        else 
        {
            mx = 0.0;
            my = 0.0;
        }
    }
    // If the coordinate id is "smoothed", evaluate the smoothed coordinate surfaces.   
    else if (config_.coord_id == "smoothed") {
        // Smoothed coordinate surfaces (simplified)
        double h_eff = apply_coordinate_smoothing(h, zeta);
        double hx_eff = hx;  // simplified - would need proper smoothing
        double hy_eff = hy;

        z_out = zeta + h_eff * (1.0 - zeta/ztop);

        J = 1.0 - h_eff/ztop;

        double zeta_factor = 1.0 - zeta/ztop;
        if (std::abs(J) > terrain_constants::epsilon) {
            mx = -zeta_factor * hx_eff / J;
            my = -zeta_factor * hy_eff / J;
        } else {
            mx = 0.0;
            my = 0.0;
        }
    } else {
        // Cartesian coordinates
        z_out = zeta;
        J = 1.0;
        mx = 0.0;
        my = 0.0;
    }
}

double TerrainFollowingCoordinate::apply_coordinate_smoothing(double h, double zeta) const 
{
    // Exponential decay of terrain influence with height
    double decay_factor = std::exp(-config_.smoothing_decay * zeta / config_.ztop);
    return h * decay_factor;
}

/*This function initializes the topography.
Takes in the topography and the number of rows and columns and initializes the topography.*/

void initialize_topography(Topography2D& topo, int NR, int NTH) 
{
    topo.h.assign(NR, std::vector<double>(NTH, 0.0));
    topo.hx.assign(NR, std::vector<double>(NTH, 0.0));
    topo.hy.assign(NR, std::vector<double>(NTH, 0.0));
}

/*This function initializes the metrics.
Takes in the metrics and the number of rows and columns and initializes the metrics.*/
void initialize_metrics(TerrainMetrics3D& metrics, int NR, int NTH, int NZ) 
{
    metrics.z.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));
    metrics.J.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 1.0)));
    metrics.mx.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));
    metrics.my.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));
    metrics.zeta.assign(NR, std::vector<std::vector<double>>(NTH, std::vector<double>(NZ, 0.0)));
}

/*This function builds the zeta levels.
Takes in the number of levels and the top of the domain and builds the zeta levels.*/
std::vector<double> build_zeta_levels(int NZ, double ztop) 
{
    std::vector<double> zeta(NZ);

    // Iterate over the vertical levels and build the zeta levels.
    for (int k = 0; k < NZ; ++k) 
    {
        zeta[k] = k * ztop / (NZ - 1);  // uniform spacing
    }
    return zeta;
}

/*This function checks the coordinate validity.
Takes in the metrics and the diagnostics and checks the coordinate validity.*/
bool check_coordinate_validity(const TerrainMetrics3D& metrics, TerrainDiagnostics& diag) 
{
    const int NR = metrics.J.size();
    const int NTH = NR > 0 ? metrics.J[0].size() : 0;
    const int NZ = NTH > 0 ? metrics.J[0][0].size() : 0;

    diag.min_jacobian = 1e6;
    diag.max_jacobian = -1e6;
    diag.coordinate_folding = false;
    diag.warnings.clear();

    // Iterate over the rows, columns, and levels and check the coordinate validity.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over the columns and check the coordinate validity.
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over the levels and check the coordinate validity.
            for (int k = 0; k < NZ; ++k) 
            {
                double J = metrics.J[i][j][k];
                diag.min_jacobian = std::min(diag.min_jacobian, J);
                diag.max_jacobian = std::max(diag.max_jacobian, J);

                // If the Jacobian is less than or equal to the epsilon, set the coordinate folding to true.
                if (J <= terrain_constants::epsilon) 
                {
                    diag.coordinate_folding = true;
                    diag.warnings.push_back("Coordinate folding detected at (i,j,k)=(" +
                                          std::to_string(i) + "," + std::to_string(j) + "," +
                                          std::to_string(k) + "), J=" + std::to_string(J));
                }
            }
        }
    }

    return !diag.coordinate_folding;
}

} // namespace topography
