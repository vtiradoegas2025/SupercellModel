#pragma once
#include "../../base/surface_fluxes.hpp"
#include "boundary_layer_base.hpp"

/*This class implements the slab boundary layer scheme.
This scheme is a one-dimensional boundary layer scheme that 
is used to compute the boundary layer physics of the simulation.
It is based on the slab scheme developed by Betts (1973).*/
class SlabScheme : public BoundaryLayerSchemeBase 
{
private:
    // Slab model state (prognostic)
    double h_ = 500.0;        // PBL height [m]
    double theta_m_ = 300.0;  // mixed-layer potential temperature [K]
    double qv_m_ = 0.01;      // mixed-layer moisture [kg/kg]

    // Configuration parameters
    double entrainment_coeff_ = 0.2;  // entrainment efficiency
    double min_h_ = 100.0;            // minimum PBL height [m]
    double max_h_ = 2000.0;           // maximum PBL height [m]

public:
    SlabScheme();

    std::string name() const override { return "slab"; }
    int required_fields() const override;

    void initialize(const BoundaryLayerConfig& cfg) override;

    void compute_column(
        const BoundaryLayerConfig& cfg,
        const SurfaceConfig& sfc,
        const BoundaryLayerColumnStateView& col,
        BoundaryLayerTendencies& tend,
        BoundaryLayerDiagnostics* diag_opt = nullptr) override;

private:
    // Diagnose PBL height from Richardson number
    double diagnose_pbl_height(
        const BoundaryLayerColumnStateView& col,
        const BoundaryLayerConfig& cfg
    );

    // Update slab state with entrainment
    void update_slab_state(
        const BoundaryLayerColumnStateView& col,
        const SurfaceConfig& sfc,
        const surface_fluxes::BulkFluxes& fluxes,
        double dt
    );

    // Apply tendencies to grid
    void apply_slab_tendencies(
        const BoundaryLayerColumnStateView& col,
        BoundaryLayerTendencies& tend,
        const BoundaryLayerConfig& cfg
    );
};
