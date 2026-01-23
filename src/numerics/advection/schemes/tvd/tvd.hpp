#pragma once
#include "advection_base.hpp"

/*This class implements the TVD advection scheme.
It is a subclass of the AdvectionSchemeBase class.
It implements the compute_flux_divergence method.*/
class TVDScheme : public AdvectionSchemeBase 
{
private:
    AdvectionConfig config_;
    double (*limiter_function_)(double) = nullptr;  // pointer to limiter function

    // Limiter functions
    static double minmod_limiter(double r);
    static double vanleer_limiter(double r);
    static double superbee_limiter(double r);
    static double mc_limiter(double r);
    static double universal_limiter(double r);

public:
    TVDScheme();

    std::string name() const override { return "tvd"; }
    void initialize() override;
    void initialize(const AdvectionConfig& cfg) override;

    /*This function computes the flux divergence.
    Takes in the configuration, state, tendencies, and diagnostics and computes the flux divergence.*/
    void compute_flux_divergence(
        const AdvectionConfig& cfg,
        const AdvectionStateView& state,
        AdvectionTendencies& tendencies,
        AdvectionDiagnostics* diag_opt = nullptr) override;

    /*This function suggests the time step.
    Takes in the configuration and state and suggests the time step.*/
    double suggest_dt(
        const AdvectionConfig& cfg,
        const AdvectionStateView& state) override;

private:
    // MUSCL reconstruction with limiter
    void muscl_reconstruct(
        const std::vector<double>& q,
        std::vector<double>& q_left,
        std::vector<double>& q_right,
        double (*limiter)(double));

    // Compute numerical flux (upwind)
    double numerical_flux(double q_left, double q_right, double velocity);

    // 1D advection in one direction
    void advect_1d(
        const std::vector<double>& q,
        const std::vector<double>& velocity,
        double dx, double dt,
        std::vector<double>& dqdt);

    // Apply positivity limiter if enabled
    void apply_positivity_limiter(std::vector<double>& q, double min_value = 0.0);
};
