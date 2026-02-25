/**
 * @file bell.hpp
 * @brief Declarations for the terrain module.
 *
 * Defines interfaces, data structures, and contracts used by
 * the terrain runtime and scheme implementations.
 * This file is part of the src/terrain subsystem.
 */

#pragma once
#include "terrain/base/topography.hpp"

/**
 * @brief Bell-mountain terrain scheme for idealized topography tests.
 */
class BellScheme : public TerrainSchemeBase 
{
private:
    TerrainConfig config_;

public:
    /**
     * @brief Constructs the bell terrain scheme.
     */
    BellScheme();
    /**
 * @brief Gets the name of the bell terrain scheme.
 */
    std::string name() const override { return "bell"; }

    /**
     * @brief Initializes scheme parameters from terrain configuration.
     */
    void initialize(const TerrainConfig& cfg) override;

    /**
 * @brief Builds the topography.
 */
    void build_topography(const TerrainConfig& cfg, Topography2D& topo) override;

    /**
 * @brief Builds the metrics.
 */
    void build_metrics(const TerrainConfig& cfg,
                      const Topography2D& topo,
                      TerrainMetrics3D& metrics,
                      TerrainDiagnostics* diag_opt = nullptr) override;

private:
    /**
     * @brief Computes physical x/y coordinates for a grid index pair.
     */
    void get_coordinates(int i, int j, double& x, double& y) const;
};
