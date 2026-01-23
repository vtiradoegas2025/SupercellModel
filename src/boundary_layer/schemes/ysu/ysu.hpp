#pragma once
#include "../../base/surface_fluxes.hpp"
#include "boundary_layer_base.hpp"

/*This class implements the YSU boundary layer scheme.
This scheme is a one-dimensional boundary layer scheme that 
is used to compute the boundary layer physics of the simulation.
It is based on the YSU scheme developed by Hong et al. (2006).*/
class YSUScheme : public BoundaryLayerSchemeBase {
private:
    // YSU scheme parameters (Hong et al. 2006)
    double pblfac_ = 1.0;     // PBL factor
    double cn_ = 0.75;        // coefficient for eddy diffusivity
    double ck_ = 0.1;         // coefficient for nonlocal transport
    double ce_ = 0.5;         // entrainment efficiency
    double c0_ = 0.15;        // surface exchange coefficient

    // Scheme state
    std::vector<double> K_m_;     // momentum diffusivity profile
    std::vector<double> K_h_;     // heat diffusivity profile

public:
    YSUScheme();

    std::string name() const override { return "ysu"; }
    int required_fields() const override;

    void initialize(const BoundaryLayerConfig& cfg) override;

    void compute_column(
        const BoundaryLayerConfig& cfg,
        const SurfaceConfig& sfc,
        const BoundaryLayerColumnStateView& col,
        BoundaryLayerTendencies& tend,
        BoundaryLayerDiagnostics* diag_opt = nullptr) override;

private:
    // Diagnose PBL height
    double diagnose_pbl_height(
        const BoundaryLayerColumnStateView& col,
        const BoundaryLayerConfig& cfg
    );

    // Compute K-profile diffusivities
    void compute_k_profile(
        const BoundaryLayerColumnStateView& col,
        double h, double ustar,
        std::vector<double>& K_m,
        std::vector<double>& K_h
    );

    // Compute nonlocal transport term
    void compute_nonlocal_transport(
        const BoundaryLayerColumnStateView& col,
        const surface_fluxes::BulkFluxes& fluxes,
        double h,
        std::vector<double>& nonlocal_theta,
        std::vector<double>& nonlocal_qv
    );

    // Apply diffusion tendencies
    void apply_diffusion_tendencies(
        const BoundaryLayerColumnStateView& col,
        const std::vector<double>& K_m,
        const std::vector<double>& K_h,
        BoundaryLayerTendencies& tend
    );
};
