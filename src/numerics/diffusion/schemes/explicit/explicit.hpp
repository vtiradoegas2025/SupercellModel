#pragma once
#include "diffusion_base.hpp"

/*This class implements the explicit diffusion scheme.
It is a subclass of the DiffusionSchemeBase class.
It implements the compute_diffusion_tendencies method.*/
class ExplicitDiffusionScheme : public DiffusionSchemeBase 
{
private:
    DiffusionConfig config_;

public:
    ExplicitDiffusionScheme();

    std::string name() const override { return "explicit"; }
    void initialize() override;
    void initialize(const DiffusionConfig& cfg) override;

    /*This function computes the diffusion tendencies.
    Takes in the configuration, state, tendencies, and diagnostics 
    and computes the diffusion tendencies.*/
    void compute_diffusion_tendencies(
        const DiffusionConfig& cfg,
        const DiffusionStateView& state,
        DiffusionTendencies& tendencies,
        DiffusionDiagnostics* diag_opt = nullptr) override;

    /*This function checks the stability of the diffusion scheme.
    Takes in the configuration and state and checks the stability of the diffusion scheme.*/
    double check_stability(
        const DiffusionConfig& cfg,
        const DiffusionStateView& state) override;

private:

    /*This function computes the scalar diffusion.
    Takes in the field, the diffusivity field, the grid, and the tendency 
    and computes the scalar diffusion.*/
    // Compute diffusion tendency for a single field
    void compute_scalar_diffusion(
        const std::vector<std::vector<std::vector<double>>>& field,
        const std::vector<std::vector<std::vector<double>>>& K_field,
        const GridMetrics& grid,
        std::vector<std::vector<std::vector<double>>>& tendency);

    /*This function computes the momentum diffusion.
    Takes in the u, v, w velocities, the diffusivity field, the grid, 
    and the tendency and computes the momentum diffusion.*/
    void compute_momentum_diffusion(
        const std::vector<std::vector<std::vector<double>>>& u,
        const std::vector<std::vector<std::vector<double>>>& v,
        const std::vector<std::vector<std::vector<double>>>& w,
        const std::vector<std::vector<std::vector<double>>>& nu_t,
        const GridMetrics& grid,
        std::vector<std::vector<std::vector<double>>>& du_dt,
        std::vector<std::vector<std::vector<double>>>& dv_dt,
        std::vector<std::vector<std::vector<double>>>& dw_dt);
};
