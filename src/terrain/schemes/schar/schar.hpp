#pragma once
#include "../../base/topography.hpp"

/*This file contains the declaration of the Schär terrain scheme.
It manages the initialization of the Schär terrain scheme and the building of the topography.*/ 

class ScharScheme : public TerrainSchemeBase 
{
private:
    TerrainConfig config_;

public:
    ScharScheme();

    /*This function gets the name of the Schär terrain scheme.
    Takes in the name and returns the name of the Schär terrain scheme.*/
    std::string name() const override { return "schar"; }

    void initialize(const TerrainConfig& cfg) override;

    /*This function builds the topography.
    Takes in the configuration and the topography and builds the topography.*/
    void build_topography(const TerrainConfig& cfg, Topography2D& topo) override;

    void build_metrics(const TerrainConfig& cfg,
                      const Topography2D& topo,
                      TerrainMetrics3D& metrics,
                      TerrainDiagnostics* diag_opt = nullptr) override;

private:
    // Compute coordinates for grid point
    void get_coordinates(int i, int j, double& x, double& y) const;
};
