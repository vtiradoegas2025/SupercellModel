#pragma once
#include "chaos_base.hpp"
#include "../../base/random_generator.hpp"
#include "../../base/correlation_filter.hpp"

namespace chaos {

/**
 * @brief Initial conditions perturbation scheme
 *
 * Applies stochastic perturbations to initial conditions (IC) only.
 * Represents uncertainty in the initial storm environment without ongoing
 * model error forcing. Based on AMS ensemble practices.
 */
class InitialConditionsScheme : public ChaosScheme {
public:
    std::string name() const override { return "initial_conditions"; }

    void initialize(const ChaosConfig& cfg, const GridMetrics& grid) override;

    void apply_initial_conditions(
        const ChaosConfig& cfg,
        const GridMetrics& grid,
        SimulationState& state,
        ChaosDiagnostics* diag = nullptr
    ) override;

    // No tendency perturbations for IC-only scheme
    void apply_tendencies(
        const ChaosConfig& cfg,
        const GridMetrics& grid,
        const ChaosStateView& state_view,
        ChaosTendencies& tendencies,
        ChaosDiagnostics* diag = nullptr
    ) override;

    // No temporal evolution needed
    void step_noise(
        const ChaosConfig& cfg,
        const GridMetrics& grid,
        double dt
    ) override;

private:
    ChaosRNG rng_;
    std::unique_ptr<CorrelationFilter> correlation_filter_;
};

} // namespace chaos
