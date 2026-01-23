#include "simulation.hpp"
#include "factory.hpp"
#include <iostream>
#include <memory>

/*This file contains the implementation of the chaos perturbation module.
It manages stochastic perturbations for ensemble-style sensitivity experiments
and "chaos seeding" for storm predictability studies in convection-permitting setups.
*/

// Global chaos scheme instance
std::unique_ptr<chaos::ChaosScheme> chaos_scheme = nullptr;

// Global chaos configuration
chaos::ChaosConfig global_chaos_config;

// Chaos diagnostics
chaos::ChaosDiagnostics global_chaos_diagnostics;

/**
 * @brief Initialize the chaos perturbation scheme
 * @param cfg Chaos configuration from YAML or defaults
 */
void initialize_chaos(const chaos::ChaosConfig& cfg) 
{
    try {
        global_chaos_config = cfg;
        chaos_scheme = create_chaos_scheme(cfg.scheme_id);

        // Initialize with grid metrics (will be set after grid initialization)
        // For now, create a dummy grid - will be properly initialized later
        GridMetrics dummy_grid;
        chaos_scheme->initialize(cfg, dummy_grid);

        std::cout << "Initialized chaos scheme: " << cfg.scheme_id << std::endl;

        if (cfg.scheme_id != "none") {
            std::cout << "  Perturbation amplitudes:" << std::endl;
            if (!cfg.apply_to_ic.empty()) {
                std::cout << "    IC variables: ";
                for (const auto& var : cfg.apply_to_ic) {
                    auto it = cfg.sigma_ic.find(var);
                    if (it != cfg.sigma_ic.end()) {
                        std::cout << var << "(σ=" << it->second << ") ";
                    }
                }
                std::cout << std::endl;
            }
            if (!cfg.apply_to_tendencies.empty()) {
                std::cout << "    Tendency blocks: ";
                for (const auto& block : cfg.apply_to_tendencies) {
                    auto it = cfg.alpha_tend.find(block);
                    if (it != cfg.alpha_tend.end()) {
                        std::cout << block << "(α=" << it->second << ") ";
                    }
                }
                std::cout << std::endl;
            }
            std::cout << "  Correlation: Lx=" << cfg.Lx << "m, Ly=" << cfg.Ly << "m" << std::endl;
            std::cout << "  Temporal decorrelation: τ=" << cfg.tau_t << "s" << std::endl;
        }

    } catch (const std::exception& e) {
        std::cerr << "Error initializing chaos: " << e.what() << std::endl;
        throw;
    }
}

/**
 * @brief Apply initial condition perturbations
 * Called once at simulation startup (t=0)
 */
void apply_chaos_initial_conditions() 
{
    if (!chaos_scheme) 
    {
        std::cerr << "Warning: Chaos scheme not initialized" << std::endl;
        return;
    }

    try {
        // Create state view for chaos scheme
        chaos::ChaosStateView state_view;
        state_view.grid = &global_grid_metrics;
        state_view.u = &u;
        state_view.v_theta = &v_theta;
        state_view.w = &w;
        state_view.theta = &theta;
        state_view.qv = &qv;
        state_view.qc = &qc;
        state_view.qr = &qr;

        // Apply perturbations (state is modified in-place)
        // TODO: Proper SimulationState interface
        // This is a COMEBACK SECTION - using null pointer placeholder for initial integration
        //
        // Full implementation should:
        // 1. Create proper SimulationState object with references to global fields
        // 2. Or modify scheme interface to work with global variables directly
        // 3. Ensure thread safety and proper state management
        //
        // Current limitation: Scheme interface expects SimulationState but globals are used directly

        chaos_scheme->apply_initial_conditions(
            global_chaos_config,
            global_grid_metrics,
            *(SimulationState*)nullptr,  // TODO: Replace with proper state object
            &global_chaos_diagnostics
        );

        // Log diagnostics
        if (!global_chaos_diagnostics.warnings.empty()) 
        {
            std::cout << "Chaos IC warnings:" << std::endl;
            for (const auto& warning : global_chaos_diagnostics.warnings) 
            {
                std::cout << "  " << warning << std::endl;
            }
        }

    }
    catch (const std::exception& e) 
    {
        std::cerr << "Error applying chaos IC perturbations: " << e.what() << std::endl;
        throw;
    }
}

/**
 * @brief Apply tendency perturbations
 * Called each physics timestep before tendencies are applied
 */
void apply_chaos_tendencies() 
{
    if (!chaos_scheme) 
    {
        return;  // Silently skip if not initialized
    }

    try 
    {
        // Create state view
        chaos::ChaosStateView state_view;
        state_view.grid = &global_grid_metrics;
        state_view.u = &u;
        state_view.v_theta = &v_theta;
        state_view.w = &w;
        state_view.theta = &theta;
        state_view.qv = &qv;
        state_view.qc = &qc;
        state_view.qr = &qr;
        state_view.pbl_height = 1000.0;  // Default PBL height

        // Create tendencies view (this will need to be populated with actual tendencies)
        chaos::ChaosTendencies tendencies;
        // Note: In a full implementation, these would point to the actual physics tendencies
        // For now, this is a placeholder structure

        // Apply perturbations
        chaos_scheme->apply_tendencies(
            global_chaos_config,
            global_grid_metrics,
            state_view,
            tendencies,
            &global_chaos_diagnostics
        );

    } 
    catch (const std::exception& e) 
    {
        std::cerr << "Error applying chaos tendency perturbations: " << e.what() << std::endl;
        // Don't throw here - physics should continue
    }
}

/**
 * @brief Evolve chaos noise fields
 * Called each timestep for schemes with temporal evolution
 * @param dt Timestep size
 */
void step_chaos_noise(double dt) 
{
    if(!chaos_scheme) 
    {
        return;
    }

    try 
    {
        chaos_scheme->step_noise(global_chaos_config, global_grid_metrics, dt);
    } 
    catch (const std::exception& e) 
    {
        std::cerr << "Error stepping chaos noise: " << e.what() << std::endl;
        // Don't throw - simulation should continue
    }
}

/**
 * @brief Get current chaos diagnostics
 */
const chaos::ChaosDiagnostics& get_chaos_diagnostics() 
{
    return global_chaos_diagnostics;
}

/**
 * @brief Reset chaos diagnostics
 */
void reset_chaos_diagnostics() 
{
    global_chaos_diagnostics = chaos::ChaosDiagnostics{};
}
