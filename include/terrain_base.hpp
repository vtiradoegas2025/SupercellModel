#pragma once
#include <vector>
#include <memory>
#include <string>

/*This header file contains the base classes and structures for the terrain module.
The terrain module is responsible for the terrain of the simulation.
The terrain scheme is chosen by the user in the configuration file.
This module is used to build the terrain of the simulation.*/

// Terrain module constants
namespace terrain_constants 
{
    inline constexpr double epsilon = 1e-12;        // small number for divisions
}

// Configuration for terrain schemes
struct TerrainConfig 
{
    std::string scheme_id = "none";
    std::string coord_id = "btf";  // "btf" (basic terrain-following) or "smoothed"
    double ztop = 20000.0;         // model top height [m]
    bool compute_derivatives = true;  // compute h_x, h_y
    bool compute_metrics = true;      // compute J, m_x, m_y, z

    // Terrain-specific parameters
    struct BellParams 
    {
        double h0 = 1000.0;         // maximum height [m]
        double a = 5000.0;          // half-width [m]
        bool axisymmetric = false;  // 3D bell vs 2D ridge
    } bell;

    struct ScharParams 
    {
        double H = 1000.0;          // maximum height [m]
        double a = 5000.0;          // envelope width [m]
        double ell = 4000.0;        // small-scale wavelength [m]
    } schar;

    // Coordinate smoothing (for coord_id="smoothed")
    double smoothing_decay = 0.1;  // decay rate with height
};

// 2D topography field
struct Topography2D 
{
    std::vector<std::vector<double>> h;    // terrain height [m] (NR x NTH)
    std::vector<std::vector<double>> hx;   // dh/dx [m/m] (optional)
    std::vector<std::vector<double>> hy;   // dh/dy [m/m] (optional)
};

// 3D terrain-following coordinate metrics
struct TerrainMetrics3D 
{
    std::vector<std::vector<std::vector<double>>> z;     // physical height [m] (NR x NTH x NZ)
    std::vector<std::vector<std::vector<double>>> J;     // Jacobian dz/dζ (NR x NTH x NZ)
    std::vector<std::vector<std::vector<double>>> mx;    // metric term -dz_x/dz_ζ (NR x NTH x NZ)
    std::vector<std::vector<std::vector<double>>> my;    // metric term -dz_y/dz_ζ (NR x NTH x NZ)
    std::vector<std::vector<std::vector<double>>> zeta;  // terrain-following coordinate [m] (NR x NTH x NZ)
};

// Diagnostics for terrain
struct TerrainDiagnostics 
{
    double max_height = 0.0;        // maximum terrain height [m]
    double max_slope_x = 0.0;       // maximum |dh/dx|
    double max_slope_y = 0.0;       // maximum |dh/dy|
    double min_jacobian = 1.0;      // minimum Jacobian (should be > 0)
    double max_jacobian = 1.0;      // maximum Jacobian
    bool coordinate_folding = false; // true if J <= 0 anywhere
    std::vector<std::string> warnings; // warning messages
};

// Abstract base class for terrain schemes
class TerrainSchemeBase 
{
public:
    virtual ~TerrainSchemeBase() = default;

    virtual std::string name() const = 0;

    // Initialize with configuration
    virtual void initialize(const TerrainConfig& cfg) = 0;

    // Build 2D topography field
    virtual void build_topography(const TerrainConfig& cfg, Topography2D& topo) = 0;

    // Build 3D coordinate metrics from topography
    virtual void build_metrics(const TerrainConfig& cfg,
                              const Topography2D& topo,
                              TerrainMetrics3D& metrics,
                              TerrainDiagnostics* diag_opt = nullptr) = 0;
};

// Factory function type
using TerrainSchemeFactory = std::unique_ptr<TerrainSchemeBase> (*)(const TerrainConfig&);

// Global terrain scheme instance and configuration
extern std::unique_ptr<TerrainSchemeBase> terrain_scheme;
extern TerrainConfig global_terrain_config;
extern Topography2D global_topography;
extern TerrainMetrics3D global_terrain_metrics;

// Factory and initialization functions
std::unique_ptr<TerrainSchemeBase> create_terrain_scheme(const std::string& scheme_name);
std::vector<std::string> get_available_terrain_schemes();
void initialize_terrain(const std::string& scheme_name, const TerrainConfig& cfg = TerrainConfig{});
void build_terrain_fields();
