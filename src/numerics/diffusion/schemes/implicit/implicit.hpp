#pragma once
#include "diffusion_base.hpp"

class ImplicitDiffusionScheme : public DiffusionSchemeBase {
private:
    DiffusionConfig config_;

    /*This function solves the tridiagonal system.
    Takes in the a, b, c, rhs, and x and solves the tridiagonal system.*/
    void solve_tridiagonal(
        const std::vector<double>& a,
        const std::vector<double>& b,
        const std::vector<double>& c,
        const std::vector<double>& rhs,
        std::vector<double>& x);

public:
    ImplicitDiffusionScheme();

    std::string name() const override { return "implicit"; }
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
    // Implicit vertical diffusion for scalars

    /*This function computes the vertical diffusion.
    Takes in the field, the diffusivity field, the grid, the time step, 
    and the new field and computes the vertical diffusion.*/
    void implicit_vertical_diffusion(
        const std::vector<std::vector<std::vector<double>>>& field,
        const std::vector<std::vector<std::vector<double>>>& K_field,
        const GridMetrics& grid,
        double dt,
        std::vector<std::vector<std::vector<double>>>& field_new);

    /*This function computes the vertical momentum diffusion.
    Takes in the u, v, w velocities, the diffusivity field, the grid, 
    the time step, the new u, v, and w velocities, and computes the vertical momentum diffusion.*/
    void implicit_vertical_momentum_diffusion(
        const std::vector<std::vector<std::vector<double>>>& u,
        const std::vector<std::vector<std::vector<double>>>& v,
        const std::vector<std::vector<std::vector<double>>>& w,
        const std::vector<std::vector<std::vector<double>>>& nu_t,
        const GridMetrics& grid,
        double dt,
        std::vector<std::vector<std::vector<double>>>& u_new,
        std::vector<std::vector<std::vector<double>>>& v_new,
        std::vector<std::vector<std::vector<double>>>& w_new);
};
