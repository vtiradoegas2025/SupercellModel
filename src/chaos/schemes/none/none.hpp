#pragma once
#include "chaos_base.hpp"

namespace chaos {

/**
 * @brief "None" chaos scheme - no perturbations applied
 *
 * This scheme serves as a baseline and does nothing. Useful for deterministic
 * runs or when chaos perturbations are disabled.
 */
class NoneScheme : public ChaosScheme {
public:
    std::string name() const override { return "none"; }

    void initialize(const ChaosConfig& cfg, const GridMetrics& grid) override;

    // IC perturbations - none applied
    void apply_initial_conditions(
        const ChaosConfig& cfg,
        const GridMetrics& grid,
        SimulationState& state,
        ChaosDiagnostics* diag = nullptr
    ) override;

    // Tendency perturbations - none applied
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
};

} // namespace chaos
