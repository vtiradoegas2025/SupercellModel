#include "rk4.hpp"
#include <iostream>

RK4Scheme::RK4Scheme() : rhs_function_(nullptr) {}

/*This function initializes the RK4 time stepping scheme.
Takes in the configuration and initializes the RK4 time stepping scheme.*/
void RK4Scheme::initialize() 
{
    initialize(TimeSteppingConfig{}, nullptr);
}

/*This function initializes the RK4 time stepping scheme.
Takes in the configuration and the right-hand side function and initializes the RK4 time stepping scheme.*/
void RK4Scheme::initialize(const TimeSteppingConfig& cfg, RHSFunction rhs_func) 
{
    config_ = cfg;
    rhs_function_ = rhs_func;

    std::cout << "Initialized RK4 time stepping scheme" << std::endl;
}

/*This function steps the RK4 time stepping scheme.
Takes in the configuration, state, and diagnostics and steps the RK4 time stepping scheme.*/
void RK4Scheme::step(
    const TimeSteppingConfig& cfg,
    TimeSteppingState& state,
    TimeSteppingDiagnostics* diag_opt
) 
{
    // If the right-hand side function is not set, throw an error.
    if (!rhs_function_) 
    {
        throw std::runtime_error("RHS function not set for RK4 scheme");
    }

    // If the stage storage is not set, resize the stage storage.
    if (stage1_.empty() || stage1_.size() != state.fields.size()) {
        stage1_.resize(state.fields.size());
        stage2_.resize(state.fields.size());
        stage3_.resize(state.fields.size());
        k1_.resize(state.fields.size());
        k2_.resize(state.fields.size());
        k3_.resize(state.fields.size());
        k4_.resize(state.fields.size());
    }

    // Compute RK4 stages
    compute_k1(state);
    compute_k2(state);
    compute_k3(state);
    compute_k4(state);
    update_final(state);

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
double RK4Scheme::suggest_dt(
    const TimeSteppingConfig& cfg,
    const TimeSteppingState& state
) 
{
    // RK4 typically allows CFL ≈ 0.5-0.8 depending on the spatial scheme
    return 0.5;  // Conservative estimate
}

/*This function computes the first stage.
Takes in the state and the tendencies and computes the first stage.*/
void RK4Scheme::compute_k1(const TimeSteppingState& state) 
{
    // k1 = L(q^n)
    rhs_function_(state.fields, state.time, k1_);
}

/*This function computes the second stage.
Takes in the state and the tendencies and computes the second stage.*/
void RK4Scheme::compute_k2(const TimeSteppingState& state) 
{
    // q_temp = q^n + (Δt/2) * k1
    add_states(state.fields, 1.0, k1_, 0.5, state.dt, stage1_);

    // k2 = L(q_temp)
    rhs_function_(stage1_, state.time + 0.5 * state.dt, k2_);
}

/*This function computes the third stage.
Takes in the state and the tendencies and computes the third stage.*/
void RK4Scheme::compute_k3(const TimeSteppingState& state) 
{
    // q_temp = q^n + (Δt/2) * k2
    add_states(state.fields, 1.0, k2_, 0.5, state.dt, stage2_);

    // k3 = L(q_temp)
    rhs_function_(stage2_, state.time + 0.5 * state.dt, k3_);
}

void RK4Scheme::compute_k4(const TimeSteppingState& state) {
    // q_temp = q^n + Δt * k3
    add_states(state.fields, 1.0, k3_, 1.0, state.dt, stage3_);

    // k4 = L(q_temp)
    rhs_function_(stage3_, state.time + state.dt, k4_);
}

/*This function updates the final stage.
Takes in the state and updates the final stage.*/
void RK4Scheme::update_final(TimeSteppingState& state) 
{
    // q^{n+1} = q^n + (Δt/6) * (k1 + 2*k2 + 2*k3 + k4)

    // Iterate over the fields and update the final stage.  
    for (size_t i = 0; i < state.fields.size(); ++i) 
    {

        // Iterate over the vertical levels and update the final stage.
        for (size_t k = 0; k < state.fields[i].data.size(); ++k) 
        {
            // Iterate over the horizontal columns and vertical levels and update the final stage.
            for (size_t l = 0; l < state.fields[i].data[k].size(); ++l) 
            {
                // Iterate over the vertical levels and update the final stage.
                for (size_t m = 0; m < state.fields[i].data[k][l].size(); ++m) 
                {
                    // Get the values of the first, second, third, and fourth stages.
                    double k1_val = k1_[i].tendencies[k][l][m];
                    double k2_val = k2_[i].tendencies[k][l][m];
                    double k3_val = k3_[i].tendencies[k][l][m];
                    double k4_val = k4_[i].tendencies[k][l][m];

                    state.fields[i].data[k][l][m] += (state.dt / 6.0) *
                        (k1_val + 2.0 * k2_val + 2.0 * k3_val + k4_val);
                }
            }
        }
    }
}

/*This function copies the state.
Takes in the source and destination states and copies the state.*/
void RK4Scheme::copy_state(const std::vector<NumericalState>& src, std::vector<NumericalState>& dst) 
{
    // Resize the destination state.
    dst.resize(src.size());

    // Iterate over the fields and copy the state.
    for (size_t i = 0; i < src.size(); ++i) 
    {
        // Copy the data.
        dst[i].data = src[i].data;  // Deep copy of 3D data
        // Copy the name.
        dst[i].name = src[i].name;
    }
}

/*This function adds the states.
Takes in the state1, the coefficient 1, the tendencies, the coefficient 2, the time step, 
and the result and adds the states.*/
void RK4Scheme::add_states(const std::vector<NumericalState>& state1, double coef1,
                          const std::vector<NumericalTendencies>& tend, double coef2, double dt,
                          std::vector<NumericalState>& result) 
{
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
