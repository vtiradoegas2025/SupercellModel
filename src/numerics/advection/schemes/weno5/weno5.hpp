#pragma once
#include "advection_base.hpp"

/*This class implements the WENO5 advection scheme.
It is a subclass of the AdvectionSchemeBase class.
It implements the compute_flux_divergence method.*/
class WENO5Scheme : public AdvectionSchemeBase {
private:
    AdvectionConfig config_;
    double weno_epsilon_;  // small parameter for smoothness indicators

public:
    WENO5Scheme();

    std::string name() const override { return "weno5"; }
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
    // WENO5 reconstruction functions
    double weno5_reconstruct_left(const std::vector<double>& q, int i);
    double weno5_reconstruct_right(const std::vector<double>& q, int i);

    // Smoothness indicators
    double smoothness_indicator(const std::vector<double>& q, int start, int end);

    // Candidate stencils for WENO5
    double q_tilde_0(const std::vector<double>& q, int i);  // stencil 0
    double q_tilde_1(const std::vector<double>& q, int i);  // stencil 1
    double q_tilde_2(const std::vector<double>& q, int i);  // stencil 2

    // Linear weights
    static constexpr double d0 = 1.0/10.0;
    static constexpr double d1 = 6.0/10.0;
    static constexpr double d2 = 3.0/10.0;

    // Numerical flux (Rusanov)
    double rusanov_flux(double q_left, double q_right, double velocity);

    // 1D WENO5 advection
    void weno5_advect_1d(
        const std::vector<double>& q,
        const std::vector<double>& velocity,
        double dx, double dt,
        std::vector<double>& dqdt);

    // Apply positivity preservation if enabled
    void apply_positivity_preservation(std::vector<double>& q, double min_value = 0.0);
};
