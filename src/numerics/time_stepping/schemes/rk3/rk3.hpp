#pragma once
#include "time_stepping_base.hpp"

class RK3Scheme : public TimeSteppingSchemeBase {
private:
    TimeSteppingConfig config_;
    RHSFunction rhs_function_;

    // Stage storage for RK3
    std::vector<NumericalState> stage1_;
    std::vector<NumericalState> stage2_;
    std::vector<NumericalTendencies> rhs_stage1_;
    std::vector<NumericalTendencies> rhs_stage2_;

public:
    RK3Scheme();

    /*This function returns the name of the RK3 time stepping scheme.
    Takes in the name of the RK3 time stepping scheme.*/
    std::string name() const override { return "rk3"; }
    void initialize() override;
    void initialize(const TimeSteppingConfig& cfg, RHSFunction rhs_func) override;

    /*This function steps the RK3 time stepping scheme.
    Takes in the configuration, state, and diagnostics and steps the RK3 time stepping scheme.*/
    void step(
        const TimeSteppingConfig& cfg,
        TimeSteppingState& state,
        TimeSteppingDiagnostics* diag_opt = nullptr) override;

    /*This function suggests the time step.
    Takes in the configuration and state and suggests the time step.*/
    double suggest_dt(
        const TimeSteppingConfig& cfg,
        const TimeSteppingState& state) override;

private:
    // RK3 stages
    /*This function computes the first stage.
    Takes in the state and the tendencies and computes the first stage.*/
    void compute_stage1(const TimeSteppingState& state, std::vector<NumericalTendencies>& tendencies);

    /*This function computes the second stage.
    Takes in the state and the tendencies and computes the second stage.*/
    void compute_stage2(const TimeSteppingState& state, std::vector<NumericalTendencies>& tendencies);

    void compute_final(TimeSteppingState& state, std::vector<NumericalTendencies>& tendencies);

    /*This function copies the state.
    Takes in the source and destination states and copies the state.*/
    void copy_state(const std::vector<NumericalState>& src, std::vector<NumericalState>& dst);

    /*This function adds the states.
    Takes in the state1, the coefficient 1, the tendencies, the coefficient 2, the time step, 
    and the result and adds the states.*/
    void add_states(const std::vector<NumericalState>& state1, double coef1,
                   const std::vector<NumericalTendencies>& tend, double coef2, double dt,
                   std::vector<NumericalState>& result);
    void update_state(std::vector<NumericalState>& state,
                     const std::vector<NumericalTendencies>& tendencies, double dt);
};
