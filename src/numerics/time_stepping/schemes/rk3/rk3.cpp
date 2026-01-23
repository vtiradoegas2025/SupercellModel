#include "rk3.hpp"
#include <iostream>
#include <algorithm>


/*This class implements the RK3 time stepping scheme.
It is a subclass of the TimeSteppingSchemeBase class.
It implements the initialize, step, suggest_dt, compute_stage1, compute_stage2, compute_final, copy_state, add_states, and update_state methods.*/

RK3Scheme::RK3Scheme() : rhs_function_(nullptr) {}

/*This function initializes the RK3 time stepping scheme.
Takes in the configuration and initializes the RK3 time stepping scheme.*/
void RK3Scheme::initialize() 
{
    initialize(TimeSteppingConfig{}, nullptr);
}

/*This function initializes the RK3 time stepping scheme.
Takes in the configuration and the right-hand side function and initializes the RK3 time stepping scheme.*/
void RK3Scheme::initialize(const TimeSteppingConfig& cfg, RHSFunction rhs_func) 
{
    config_ = cfg;
    rhs_function_ = rhs_func;

    std::cout << "Initialized RK3 (SSPRK3) time stepping scheme" << std::endl;
    if (cfg.split_acoustic) {
        std::cout << "  Split-explicit acoustic mode enabled" << std::endl;
    }
}

/*This function steps the RK3 time stepping scheme.
Takes in the configuration, state, and diagnostics and steps the RK3 time stepping scheme.*/

void RK3Scheme::step(
    const TimeSteppingConfig& cfg,
    TimeSteppingState& state,
    TimeSteppingDiagnostics* diag_opt
) 
{
    // If the right-hand side function is not set, throw an error.
    if (!rhs_function_) 
    {
        throw std::runtime_error("RHS function not set for RK3 scheme");
    }

    // If the stage storage is not set, resize the stage storage.
    if (stage1_.empty() || stage1_.size() != state.fields.size()) {
        stage1_.resize(state.fields.size());
        stage2_.resize(state.fields.size());
        rhs_stage1_.resize(state.fields.size());
        rhs_stage2_.resize(state.fields.size());
    }

    // Stage 1: q^{(1)} = q^n + Δt * L(q^n)
    compute_stage1(state, rhs_stage1_);

    // Stage 2: q^{(2)} = (3/4)q^n + (1/4)q^{(1)} + (1/4)Δt * L(q^{(1)})
    compute_stage2(state, rhs_stage2_);

    // Final: q^{n+1} = (1/3)q^n + (2/3)q^{(2)} + (2/3)Δt * L(q^{(2)})
    compute_final(state, rhs_stage2_);

    // Update time
    state.time += state.dt;

    // If the diagnostics are requested, increment the number of steps and store the time step.
    if (diag_opt) 
    {
        diag_opt->n_steps++;
        diag_opt->dt_history.push_back(state.dt);
    }
}

/*This function suggests the time step.
Takes in the configuration and state and suggests the time step.*/
double RK3Scheme::suggest_dt(
    const TimeSteppingConfig& cfg,
    const TimeSteppingState& state
) 
{
    // RK3 is stable for CFL ≈ 1.0 with appropriate spatial schemes
    return 1.0;  // Placeholder - would need CFL calculation
}

/*This function computes the first stage.
Takes in the state and the tendencies and computes the first stage.*/
void RK3Scheme::compute_stage1(const TimeSteppingState& state, std::vector<NumericalTendencies>& tendencies) 
{
    // Stage 1: RHS evaluation at current state
    rhs_function_(state.fields, state.time, tendencies);
}

/*This function computes the second stage.
Takes in the state and the tendencies and computes the second stage.*/
void RK3Scheme::compute_stage2(const TimeSteppingState& state, std::vector<NumericalTendencies>& tendencies) 
{
    // Create intermediate state: q^{(1)} = q^n + Δt * L(q^n)
    std::vector<NumericalState> intermediate_state;
    add_states(state.fields, 1.0, rhs_stage1_, 1.0, state.dt, intermediate_state);

    // Evaluate RHS at intermediate state
    rhs_function_(intermediate_state, state.time + state.dt, tendencies);
}

/*This function computes the final stage.
Takes in the state and the tendencies and computes the final stage.*/
void RK3Scheme::compute_final(TimeSteppingState& state, std::vector<NumericalTendencies>& tendencies) 
{
    // Create intermediate state: q^{(2)} = (3/4)q^n + (1/4)q^{(1)} + (1/4)Δt * L(q^{(1)})
    std::vector<NumericalState> intermediate_state;
    add_states(state.fields, 3.0/4.0, rhs_stage1_, 1.0/4.0, state.dt, intermediate_state);

    // Add the second stage contribution: q^{(2)} += (1/4)Δt * L(q^{(1)})
    // (This is already included in the add_states call above)

    // Evaluate RHS at second intermediate state
    std::vector<NumericalTendencies> rhs_stage2;
    rhs_stage2.resize(state.fields.size());
    rhs_function_(intermediate_state, state.time + state.dt, rhs_stage2);

    // Final update: q^{n+1} = (1/3)q^n + (2/3)q^{(2)} + (2/3)Δt * L(q^{(2)})
    // Iterate over the fields and compute the final update.
    for (size_t i = 0; i < state.fields.size(); ++i) {
        auto& field = state.fields[i];
        const auto& rhs_final = rhs_stage2[i];

        // q^{n+1} = (1/3)q^n + (2/3)q^{(2)} + (2/3)Δt * L(q^{(2)})
        // Since q^{(2)} = (3/4)q^n + (1/4)q^{(1)} + (1/4)Δt * L(q^{(1)})
        // We can substitute: q^{n+1} = (1/3)q^n + (2/3)[(3/4)q^n + (1/4)q^{(1)} + (1/4)Δt L(q^{(1)})] + (2/3)Δt L(q^{(2)})
        //                   = (1/3)q^n + (1/2)q^n + (1/6)q^{(1)} + (1/6)Δt L(q^{(1)}) + (2/3)Δt L(q^{(2)})
        //                   = (5/6)q^n + (1/6)q^{(1)} + (1/6)Δt L(q^{(1)}) + (2/3)Δt L(q^{(2)})

        // Iterate over the vertical levels and compute the final update.
        for (size_t k = 0; k < field.data.size(); ++k) 
        {
            // Iterate over the horizontal columns and vertical levels and compute the final update.
            for (size_t l = 0; l < field.data[k].size(); ++l) 
            {
                // Iterate over the vertical levels and compute the final update.
                for (size_t m = 0; m < field.data[k][l].size(); ++m) 
                {
                    // Get q^{(1)} contribution
                    double q1_contrib = (1.0/6.0) * (state.fields[i].data[k][l][m] +
                                                   state.dt * rhs_stage1_[i].tendencies[k][l][m]);

                    // Final update
                    field.data[k][l][m] = (5.0/6.0) * state.fields[i].data[k][l][m] +
                                         q1_contrib +
                                         (2.0/3.0) * state.dt * rhs_final.tendencies[k][l][m];
                }
            }
        }
    }
}

/*This function copies the state.
Takes in the source and destination states and copies the state.*/
void RK3Scheme::copy_state(const std::vector<NumericalState>& src, std::vector<NumericalState>& dst) 
{
    dst.resize(src.size());

    // Iterate over the fields and copy the state.
    for (size_t i = 0; i < src.size(); ++i) 
    {
        dst[i].data = src[i].data;  // Deep copy of 3D data
        dst[i].name = src[i].name;
    }
}

/*This function adds the states.
Takes in the state1, the coefficient 1, the tendencies, the coefficient 2, the time step, 
and the result and adds the states.*/
void RK3Scheme::add_states(const std::vector<NumericalState>& state1, double coef1,
                          const std::vector<NumericalTendencies>& tend, double coef2, double dt,
                          std::vector<NumericalState>& result) {
    result.resize(state1.size());

    // Iterate over the fields and add the states.
    for (size_t i = 0; i < state1.size(); ++i) 
    {
        result[i].name = state1[i].name;
        result[i].data.resize(state1[i].data.size());

        // Iterate over the vertical levels and add the states.
        for (size_t k = 0; k < state1[i].data.size(); ++k) 
        {
            result[i].data[k].resize(state1[i].data[k].size());

            // Iterate over the horizontal columns and vertical levels and add the states.
            for (size_t l = 0; l < state1[i].data[k].size(); ++l) 
            {
                result[i].data[k][l].resize(state1[i].data[k][l].size());

                // Iterate over the vertical levels and add the states.
                for (size_t m = 0; m < state1[i].data[k][l].size(); ++m) 
                {
                    result[i].data[k][l][m] = coef1 * state1[i].data[k][l][m] +
                                            coef2 * dt * tend[i].tendencies[k][l][m];
                }
            }
        }
    }
}

/*This function updates the state.
Takes in the state, the tendencies, and the time step and updates the state.*/
void RK3Scheme::update_state(std::vector<NumericalState>& state,
                           const std::vector<NumericalTendencies>& tendencies, double dt) 
{
    // Iterate over the fields and update the state.
    for (size_t i = 0; i < state.size(); ++i) 
    {
        // Iterate over the vertical levels and update the state.
        for (size_t k = 0; k < state[i].data.size(); ++k) 
        {
            // Iterate over the horizontal columns and vertical levels and update the state.
            for (size_t l = 0; l < state[i].data[k].size(); ++l) 
            {
                // Iterate over the vertical levels and update the state.
                for (size_t m = 0; m < state[i].data[k][l].size(); ++m) 
                {
                    // Update the state.
                    state[i].data[k][l][m] += dt * tendencies[i].tendencies[k][l][m];
                }
            }
        }
    }
}
