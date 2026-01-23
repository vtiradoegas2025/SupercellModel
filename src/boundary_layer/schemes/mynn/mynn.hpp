#pragma once
#include "../../base/surface_fluxes.hpp"
#include "boundary_layer_base.hpp"

/*This class implements the MYNN boundary layer scheme.
This scheme is a one-dimensional boundary layer scheme that is 
used to compute the boundary layer physics of the simulation.
It is based on the MYNN scheme developed by Hong et al. (2006).*/
class MYNNScheme : public BoundaryLayerSchemeBase 
{
private:
    // MYNN scheme parameters
    double a1_ = 1.18;        // stability function parameter
    double a2_ = 0.665;       // stability function parameter
    double b1_ = 24.0;        // stability function parameter
    double b2_ = 15.0;        // stability function parameter
    double c1_ = 0.137;       // stability function parameter
    double c2_ = 0.75;        // stability function parameter
    double c3_ = 0.352;       // stability function parameter
    double c4_ = 0.0;         // stability function parameter
    double c5_ = 0.2;         // stability function parameter

    // TKE equation coefficients
    double ce1_ = 1.0;        // dissipation coefficient
    double ce2_ = 1.33;       // dissipation coefficient

    // Mixing length parameters
    double lmax_ = 500.0;     // maximum mixing length [m]

    // Prognostic TKE field
    std::vector<double> tke_; // turbulent kinetic energy [m²/s²]

/*This constructor initializes the MYNN boundary layer scheme.*/
public:
    MYNNScheme();

    std::string name() const override { return "mynn"; }
    int required_fields() const override;

    void initialize(const BoundaryLayerConfig& cfg) override;

    void compute_column(
        const BoundaryLayerConfig& cfg,
        const SurfaceConfig& sfc,
        const BoundaryLayerColumnStateView& col,
        BoundaryLayerTendencies& tend,
        BoundaryLayerDiagnostics* diag_opt = nullptr) override;


/*This function returns the name of the MYNN boundary layer scheme.*/
private:
    // Diagnose PBL height
    double diagnose_pbl_height(
        const BoundaryLayerColumnStateView& col,
        const BoundaryLayerConfig& cfg
    );

    // Compute mixing length
    void compute_mixing_length(
        const BoundaryLayerColumnStateView& col,
        double h,
        std::vector<double>& l_mix
    );

    // Compute eddy diffusivities
    void compute_eddy_diffusivities(
        const BoundaryLayerColumnStateView& col,
        const std::vector<double>& l_mix,
        std::vector<double>& K_m,
        std::vector<double>& K_h,
        std::vector<double>& K_e
    );

    // Update TKE prognostic equation
    void update_tke(
        const BoundaryLayerColumnStateView& col,
        const surface_fluxes::BulkFluxes& fluxes,
        const std::vector<double>& K_m,
        const std::vector<double>& K_h,
        const std::vector<double>& K_e,
        double dt,
        std::vector<double>& tke_tend
    );

    // Apply diffusion tendencies
    void apply_diffusion_tendencies(
        const BoundaryLayerColumnStateView& col,
        const std::vector<double>& K_m,
        const std::vector<double>& K_h,
        BoundaryLayerTendencies& tend
    );
};
