#include "none.hpp"
#include <iostream>

namespace chaos {

void NoneScheme::initialize(const ChaosConfig& cfg, const GridMetrics& grid) {
    std::cout << "Initialized chaos scheme: none (no perturbations)" << std::endl;
}

void NoneScheme::apply_initial_conditions(
    const ChaosConfig& cfg,
    const GridMetrics& grid,
    SimulationState& state,
    ChaosDiagnostics* diag
) {
    // No perturbations applied - do nothing
    if (diag) {
        diag->mean_perturbation = 0.0;
        diag->variance_perturbation = 0.0;
        diag->warnings.push_back("None scheme: no IC perturbations applied");
    }
}

void NoneScheme::apply_tendencies(
    const ChaosConfig& cfg,
    const GridMetrics& grid,
    const ChaosStateView& state_view,
    ChaosTendencies& tendencies,
    ChaosDiagnostics* diag
) {
    // No perturbations applied - do nothing
    if (diag) {
        diag->mean_perturbation = 0.0;
        diag->variance_perturbation = 0.0;
        diag->warnings.push_back("None scheme: no tendency perturbations applied");
    }
}

void NoneScheme::step_noise(
    const ChaosConfig& cfg,
    const GridMetrics& grid,
    double dt
) {
    // No noise to evolve - do nothing
}

} // namespace chaos
