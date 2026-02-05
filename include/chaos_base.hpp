#pragma once
#include <vector>
#include <memory>
#include <string>
#include <unordered_map>
#include <functional>
#include "field3d.hpp"

/**
 * @file chaos_base.hpp
 * @brief Base definitions for stochastic chaos perturbations in convection-permitting models
 *
 * This module implements controlled stochastic variability for ensemble-style
 * sensitivity experiments and "chaos seeding" for storm predictability studies.
 *
 * CURRENT STATUS (v0.1 - Initial Integration):
 * -  Core infrastructure: ChaosConfig, ChaosStateView, ChaosEffects, ChaosScheme
 * - Random number generation with reproducible seeding
 * - Perturbation field utilities (generation, scaling, evolution, bounds)
 * - Recursive Gaussian correlation filtering
 * - Factory pattern and scheme registration
 * - Basic schemes: "none" (complete), "initial_conditions" (simplified)
 * - Integration into tornado_sim.cpp
 *
 * COMEBACK SECTIONS (Intentionally incomplete for initial integration):
 * - TODO: SpectralGaussianFilter (temporarily removed due to compilation issues)
 * - TODO: Full initial_conditions scheme implementation
 * - TODO: boundary_layer scheme implementation
 * - TODO: full_stochastic scheme implementation
 * - TODO: Proper SimulationState interface
 */

// Forward declarations
struct GridMetrics;
class SimulationState;

namespace chaos 
{

//==============================================================================
// Configuration structures
//==============================================================================

/**
 * @brief Configuration for chaos perturbation schemes
 */
struct ChaosConfig 
{
    // Scheme selection
    std::string scheme_id = "none";  // "none", "initial_conditions", "boundary_layer", "full_stochastic"

    // Reproducibility
    uint64_t seed = 42;              // Base random seed
    int member_id = 0;               // Ensemble member ID (for stream separation)

    // Temporal parameters
    double dt_chaos = 60.0;          // Chaos timestep in seconds (cadence for noise updates)
    double tau_t = 21600.0;          // Temporal correlation time (seconds, ~6 hours)

    // Spatial correlation parameters
    double Lx = 50000.0;             // Horizontal correlation length X (meters)
    double Ly = 50000.0;             // Horizontal correlation length Y (meters)
    std::string filter_id = "recursive_gaussian";  // Correlation filter type
    // TODO: Change back to "spectral_gaussian" once SpectralGaussianFilter is restored

    // Bounding and tapering
    double xi_max = 2.0;             // Maximum perturbation amplitude (for clipping)
    std::string taper_id = "none";   // Vertical tapering: "none", "pbl_only", "cosine"
    double taper_z1 = 1000.0;        // Lower taper height (m)
    double taper_z2 = 3000.0;        // Upper taper height (m)

    // What to perturb
    std::vector<std::string> apply_to_ic;        // Variables for IC perturbations: ["u", "v", "w", "theta", "qv", "qc", "qr"]
    std::vector<std::string> apply_to_tendencies; // Physics blocks for tendency perturbations: ["microphysics", "pbl", "diffusion"]

    // Perturbation amplitudes
    std::unordered_map<std::string, double> sigma_ic;    // Additive IC perturbations: var -> sigma
    std::unordered_map<std::string, double> alpha_tend;  // Multiplicative tendency perturbations: block/var -> alpha

    // Default amplitudes (from AMS literature)
    ChaosConfig() 
    {
        // IC perturbations (additive, conservative values)
        sigma_ic["u"] = 0.5;      // m/s wind perturbations
        sigma_ic["v"] = 0.5;
        sigma_ic["w"] = 0.1;
        sigma_ic["theta"] = 0.5;  // K potential temperature
        sigma_ic["qv"] = 0.0005; // kg/kg water vapor (0.5 g/kg)
        sigma_ic["qc"] = 0.0;     // cloud water (usually not perturbed IC)
        sigma_ic["qr"] = 0.0;     // rain water

        // Tendency perturbations (multiplicative, SPPT-style)
        alpha_tend["microphysics"] = 0.3;  // 30% relative perturbation
        alpha_tend["pbl"] = 0.2;           // 20% relative perturbation
        alpha_tend["diffusion"] = 0.1;     // 10% relative perturbation
    }
};

/**
 * @brief Read-only view of simulation state for chaos perturbations
 * This structure is used to view the simulation state for the chaos perturbations.
 */
struct ChaosStateView 
{
    // Grid information
    const GridMetrics* grid = nullptr;

    // State fields (read-only references)
    const Field3D* u = nullptr;     // radial wind
    const Field3D* v_theta = nullptr; // azimuthal wind
    const Field3D* w = nullptr;     // vertical wind
    const Field3D* theta = nullptr; // potential temperature
    const Field3D* qv = nullptr;    // water vapor
    const Field3D* qc = nullptr;    // cloud water
    const Field3D* qr = nullptr;    // rain water

    // Diagnostic helpers
    double pbl_height = 1000.0;  // Planetary boundary layer height (m)
    const std::vector<std::vector<double>>* terrain_mask = nullptr;  // Optional terrain mask
};

/**
 * @brief Physics tendencies that can be perturbed (SPPT-style)
 * This structure is used to store the physics tendencies for the chaos perturbations.
 */
struct ChaosTendencies 
{
    // Microphysics tendencies
    std::vector<std::vector<std::vector<float>>>* dq_micro_dt;  // microphysics moisture tendencies

    // PBL tendencies
    std::vector<std::vector<std::vector<float>>>* du_pbl_dt;    // PBL momentum tendencies
    std::vector<std::vector<std::vector<float>>>* dv_theta_pbl_dt;
    std::vector<std::vector<std::vector<float>>>* dw_pbl_dt;
    std::vector<std::vector<std::vector<float>>>* dtheta_pbl_dt; // PBL theta tendency

    // Diffusion tendencies
    std::vector<std::vector<std::vector<float>>>* du_diff_dt;   // diffusion momentum tendencies
    std::vector<std::vector<std::vector<float>>>* dv_theta_diff_dt;
    std::vector<std::vector<std::vector<float>>>* dw_diff_dt;
    std::vector<std::vector<std::vector<float>>>* dtheta_diff_dt; // diffusion theta tendency
};

/**
 * @brief Effects/outputs from chaos perturbations
 * This structure is used to store the effects of the chaos perturbations.
 */
struct ChaosEffects 
{
    // Perturbed state deltas (for IC perturbations)
    std::unordered_map<std::string, std::vector<std::vector<std::vector<float>>>> delta_state;

    // Multipliers applied to tendencies (for SPPT-style perturbations)
    std::unordered_map<std::string, std::vector<std::vector<std::vector<float>>>> multipliers;

    // Optional diagnostics
    std::vector<std::vector<std::vector<float>>> xi_field;  // Current noise field (for debugging)
    double realized_variance = 0.0;  // Realized perturbation variance
    double min_multiplier = 1.0;     // Minimum applied multiplier
    double max_multiplier = 1.0;     // Maximum applied multiplier
};

/**
 * @brief Diagnostics for chaos perturbations
 * This structure is used to store the diagnostics for the chaos perturbations.
 */
struct ChaosDiagnostics 
{
    // Statistical properties
    double mean_perturbation = 0.0;
    double variance_perturbation = 0.0;
    double correlation_length_x = 0.0;  // Realized correlation lengths
    double correlation_length_y = 0.0;

    // Physical bounds checking
    bool negative_moisture_detected = false;
    int negative_moisture_count = 0;
    double min_qv_after_perturbation = 1.0;

    // Warnings/errors
    std::vector<std::string> warnings;
    std::vector<std::string> errors;

    // Timing
    double time_noise_update = 0.0;    // Time spent updating noise (seconds)
    double time_apply_ic = 0.0;        // Time spent applying IC perturbations
    double time_apply_tend = 0.0;      // Time spent applying tendency perturbations
};

//==============================================================================
// Base class for chaos schemes
//==============================================================================

/**
 * @brief Base class for all chaos perturbation schemes
 */
class ChaosScheme 
{
public:
    virtual ~ChaosScheme() = default;

    /**
     * @brief Get scheme name/identifier
     */
    virtual std::string name() const = 0;

    /**
     * @brief Initialize the scheme with configuration
     * @param cfg Chaos configuration
     * @param grid Grid metrics
     */
    virtual void initialize(const ChaosConfig& cfg, const GridMetrics& grid) = 0;

    /**
     * @brief Apply initial condition perturbations (called once at t=0)
     * @param cfg Configuration
     * @param grid Grid metrics
     * @param state Simulation state (modified in-place)
     * @param diag Optional diagnostics output
     */
    virtual void apply_initial_conditions(
        const ChaosConfig& cfg,
        const GridMetrics& grid,
        SimulationState& state,
        ChaosDiagnostics* diag = nullptr
    ) {}

    /**
     * @brief Apply tendency perturbations (called each physics timestep)
     * @param cfg Configuration
     * @param grid Grid metrics
     * @param state_view Read-only state view
     * @param tendencies Physics tendencies (modified in-place)
     * @param diag Optional diagnostics output
     */
    virtual void apply_tendencies(
        const ChaosConfig& cfg,
        const GridMetrics& grid,
        const ChaosStateView& state_view,
        ChaosTendencies& tendencies,
        ChaosDiagnostics* diag = nullptr
    ) {}

    /**
     * @brief Update evolving noise fields (for temporal correlation)
     * @param cfg Configuration
     * @param grid Grid metrics
     * @param dt Actual timestep size
     */
    virtual void step_noise(
        const ChaosConfig& cfg,
        const GridMetrics& grid,
        double dt
    ) {}
};

//==============================================================================
// Factory function declaration
//==============================================================================

/**
 * @brief Create a chaos scheme instance
 * @param scheme_name Name of the scheme ("none", "initial_conditions", etc.)
 * @return Unique pointer to chaos scheme
 */
std::unique_ptr<ChaosScheme> create_chaos_scheme(const std::string& scheme_name);

} // namespace chaos
