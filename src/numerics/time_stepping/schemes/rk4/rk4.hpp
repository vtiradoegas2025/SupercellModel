#pragma once
#include "time_stepping_base.hpp"

class RK4Scheme : public TimeSteppingSchemeBase {
private:
    TimeSteppingConfig config_;
    RHSFunction rhs_function_;

    // Stage storage for RK4
    std::vector<NumericalState> stage1_;
    std::vector<NumericalState> stage2_;
    std::vector<NumericalState> stage3_;
    std::vector<NumericalTendencies> k1_;
    std::vector<NumericalTendencies> k2_;
    std::vector<NumericalTendencies> k3_;
    std::vector<NumericalTendencies> k4_;

public:
    RK4Scheme();

    /*This function returns the name of the RK4 time stepping scheme.
    Takes in the name of the RK4 time stepping scheme.*/
    std::string name() const override { return "rk4"; }
    void initialize() override;
    void initialize(const TimeSteppingConfig& cfg, RHSFunction rhs_func) override;

    /*This function steps the RK4 time stepping scheme.
    Takes in the configuration, state, and diagnostics and steps the RK4 time stepping scheme.*/
    void step(
        const TimeSteppingConfig& cfg,
        TimeSteppingState& state,
        TimeSteppingDiagnostics* diag_opt = nullptr) override;

    double suggest_dt(
        const TimeSteppingConfig& cfg,
        const TimeSteppingState& state) override;

private:
    // RK4 stages
    void compute_k1(const TimeSteppingState& state);
    void compute_k2(const TimeSteppingState& state);
    void compute_k3(const TimeSteppingState& state);
    void compute_k4(const TimeSteppingState& state);
    void update_final(TimeSteppingState& state);

    // Helper functions
    void copy_state(const std::vector<NumericalState>& src, std::vector<NumericalState>& dst);
    void add_states(const std::vector<NumericalState>& state1, double coef1,
                   const std::vector<NumericalTendencies>& tend, double coef2, double dt,
                   std::vector<NumericalState>& result);
};
