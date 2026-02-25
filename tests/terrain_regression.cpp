#include "simulation.hpp"
#include "terrain/base/topography.hpp"
#include "terrain/schemes/bell/bell.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

// Minimal grid globals used by terrain components in this standalone regression.
int NR = 17;
int NTH = 17;
int NZ = 12;
double dr = 200.0;
double dz = 100.0;
double dt = 0.1;
double dtheta = 2.0 * 3.14159265358979323846 / static_cast<double>(NTH);

// Terrain module calls this hook after terrain rebuilds; numerics is not linked in this test.
void refresh_grid_metrics_from_terrain() {}

namespace
{
bool nearly_equal(double a, double b, double tol = 1.0e-12)
{
    return std::abs(a - b) <= tol;
}

int expect_true(bool cond, const std::string& message)
{
    if (!cond)
    {
        std::cerr << "[terrain-regression] FAIL: " << message << std::endl;
        return 1;
    }
    return 0;
}

int test_zeta_levels_single_layer()
{
    int failures = 0;
    const auto zeta = topography::build_zeta_levels(1, 5000.0);
    failures += expect_true(zeta.size() == 1, "build_zeta_levels(NZ=1) must return one level");
    failures += expect_true(nearly_equal(zeta[0], 0.0), "single zeta level must be 0");
    return failures;
}

int test_bell_is_centered_in_domain()
{
    int failures = 0;
    NR = 17;
    NTH = 17;
    NZ = 12;
    dr = 250.0;
    dtheta = 2.0 * 3.14159265358979323846 / static_cast<double>(NTH);

    TerrainConfig cfg;
    cfg.scheme_id = "bell";
    cfg.coord_id = "btf";
    cfg.bell.h0 = 1000.0;
    cfg.bell.a = 2200.0;
    cfg.bell.axisymmetric = true;

    Topography2D topo;
    topography::initialize_topography(topo, NR, NTH);

    BellScheme scheme;
    scheme.initialize(cfg);
    scheme.build_topography(cfg, topo);

    int max_i = 0;
    int max_j = 0;
    double max_h = -1.0;
    for (int i = 0; i < NR; ++i)
    {
        for (int j = 0; j < NTH; ++j)
        {
            if (topo.h[i][j] > max_h)
            {
                max_h = topo.h[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    const int center_i = (NR - 1) / 2;
    const int center_j = (NTH - 1) / 2;
    failures += expect_true(std::abs(max_i - center_i) <= 1, "bell peak i-index must be near domain center");
    failures += expect_true(std::abs(max_j - center_j) <= 1, "bell peak j-index must be near domain center");
    failures += expect_true(max_h > 0.0, "bell peak height must be positive");
    return failures;
}

int test_low_ztop_is_sanitized_for_metrics()
{
    int failures = 0;
    NR = 25;
    NTH = 25;
    NZ = 16;
    dr = 200.0;
    dtheta = 2.0 * 3.14159265358979323846 / static_cast<double>(NTH);

    TerrainConfig cfg;
    cfg.scheme_id = "bell";
    cfg.coord_id = "terrain_following";  // alias should normalize to btf
    cfg.ztop = 500.0;                    // intentionally below terrain peak
    cfg.compute_derivatives = true;
    cfg.compute_metrics = true;
    cfg.bell.h0 = 1400.0;
    cfg.bell.a = 2000.0;
    cfg.bell.axisymmetric = true;

    initialize_terrain(cfg.scheme_id, cfg);

    double max_height = 0.0;
    for (int i = 0; i < NR; ++i)
    {
        for (int j = 0; j < NTH; ++j)
        {
            max_height = std::max(max_height, global_topography.h[i][j]);
        }
    }
    failures += expect_true(
        global_terrain_config.ztop > max_height,
        "terrain.ztop must be adjusted to exceed maximum terrain height");

    bool all_jacobians_positive = true;
    for (int i = 0; i < NR && all_jacobians_positive; ++i)
    {
        for (int j = 0; j < NTH && all_jacobians_positive; ++j)
        {
            for (int k = 0; k < NZ; ++k)
            {
                const double J = global_terrain_metrics.J(i, j, k);
                if (!std::isfinite(J) || J <= terrain_constants::epsilon)
                {
                    all_jacobians_positive = false;
                    break;
                }
            }
        }
    }
    failures += expect_true(all_jacobians_positive, "terrain Jacobian must remain positive and finite");
    failures += expect_true(global_terrain_config.coord_id == "btf", "terrain coord alias must normalize to btf");
    return failures;
}
} // namespace

int main()
{
    int failures = 0;
    failures += test_zeta_levels_single_layer();
    failures += test_bell_is_centered_in_domain();
    failures += test_low_ztop_is_sanitized_for_metrics();

    if (failures > 0)
    {
        std::cerr << "[terrain-regression] FAILED with " << failures << " check(s)." << std::endl;
        return 1;
    }

    std::cout << "[terrain-regression] all checks passed" << std::endl;
    return 0;
}
