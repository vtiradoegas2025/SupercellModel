#include "simulation.hpp"
#include "factory.hpp"
#include "field3d.hpp"
#include <iostream>
#include <memory>

/*This file contains the implementation of the chaos perturbation module.
It manages stochastic perturbations for ensemble-style sensitivity experiments
and "chaos seeding" for storm predictability studies in convection-permitting setups.
*/

// Simple SimulationState wrapper for chaos perturbations
// This is a temporary solution until proper SimulationState interface is implemented
struct SimulationState {
    // References to global simulation fields
    Field3D* u;
    Field3D* v_theta;
    Field3D* w;
    Field3D* theta;
    Field3D* qv;
    Field3D* qc;
    Field3D* qr;
    
    SimulationState() 
        : u(&::u), v_theta(&::v_theta), w(&::w), theta(&::theta), 
          qv(&::qv), qc(&::qc), qr(&::qr) {}
};

// Global chaos scheme instance
std::unique_ptr<chaos::ChaosScheme> chaos_scheme = nullptr;

// Global chaos configuration
chaos::ChaosConfig global_chaos_config;  // Defined here, declared extern in simulation.hpp

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

        // Initialize with actual grid metrics (should be set by initialize_numerics() first)
        // If grid metrics not yet initialized, use defaults
        if (global_grid_metrics.dx > 0.0) {
            chaos_scheme->initialize(cfg, global_grid_metrics);
        } else {
            // Fallback: create grid metrics from global grid variables
            GridMetrics grid;
            grid.dx = dr;
            grid.dy = dtheta * 1000.0;  // approximate arc length
            grid.dz.assign(NZ, dz);
            grid.z_int.resize(NZ + 1);
            grid.z_int[0] = 0.0;
            for (int k = 1; k <= NZ; ++k) {
                grid.z_int[k] = grid.z_int[k-1] + dz;
            }
            chaos_scheme->initialize(cfg, grid);
        }

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
        // Create SimulationState wrapper with references to global fields
        SimulationState sim_state;
        
        chaos_scheme->apply_initial_conditions(
            global_chaos_config,
            global_grid_metrics,
            sim_state,
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
        
        // Debug: Check values after chaos perturbations
        float theta_min = 1e10, theta_max = -1e10;
        int nan_count = 0, inf_count = 0;
        for (int i = 0; i < NR; ++i) {
            for (int j = 0; j < NTH; ++j) {
                for (int k = 0; k < NZ; ++k) {
                    float theta_val = theta[i][j][k];
                    if (std::isnan(theta_val)) nan_count++;
                    if (std::isinf(theta_val)) inf_count++;
                    if (theta_val < theta_min) theta_min = theta_val;
                    if (theta_val > theta_max) theta_max = theta_val;
                }
            }
        }
        std::cout << "\n[CHAOS DEBUG] After chaos perturbations:" << std::endl;
        std::cout << "  Theta: min=" << theta_min << "K, max=" << theta_max << "K" << std::endl;
        std::cout << "  NaN count: " << nan_count << ", Inf count: " << inf_count << std::endl;
        if (theta_min < 0 || theta_max > 500) {
            std::cerr << "  ⚠️  WARNING: Theta corrupted by chaos perturbations!" << std::endl;
        }
        std::cout << std::endl;

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
