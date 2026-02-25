/**
 * @file none.hpp
 * @brief Declarations for the terrain module.
 *
 * Defines interfaces, data structures, and contracts used by
 * the terrain runtime and scheme implementations.
 * This file is part of the src/terrain subsystem.
 */

#pragma once
#include "terrain/base/topography.hpp"

  

/**
 * @brief Flat-terrain scheme that emits zero topography and identity metrics.
 */
class NoneScheme : public TerrainSchemeBase {
public:
    /**
 * @brief Gets the name of the none terrain scheme.
 */
    NoneScheme();

    /**
 * @brief Gets the name of the none terrain scheme.
 */
    std::string name() const override { return "none"; }

    /**
 * @brief Initializes the none terrain scheme.
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
};
