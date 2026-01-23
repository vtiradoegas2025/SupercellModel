#pragma once
#include "chaos_base.hpp"
#include "../../base/random_generator.hpp"
#include "../../base/correlation_filter.hpp"

namespace chaos {

/**
 * @brief Full stochastic perturbation scheme
 *
 * Represents time-evolving model uncertainty (subgrid physics errors) in
 * convection-permitting ensembles. Uses SPPT-like multiplicative noise
 * applied to physics tendencies. Based on AMS research on stochastic
 * physics in CAM ensembles.
 */
class FullStochasticScheme : public ChaosScheme {
public:
    std::string name() const override { return "full_stochastic"; }

    void initialize(const ChaosConfig& cfg, const GridMetrics& grid) override;

    void apply_initial_conditions(
        const ChaosConfig& cfg,
        const GridMetrics& grid,
        SimulationState& state,
        ChaosDiagnostics* diag = nullptr
    ) override;

    void apply_tendencies(
        const ChaosConfig& cfg,
        const GridMetrics& grid,
        const ChaosStateView& state_view,
        ChaosTendencies& tendencies,
        ChaosDiagnostics* diag = nullptr
    ) override;

    void step_noise(
        const ChaosConfig& cfg,
        const GridMetrics& grid,
        double dt
    ) override;

private:
    ChaosRNG rng_;
    std::unique_ptr<CorrelationFilter> correlation_filter_;

    // Evolving noise fields for each physics block
    std::vector<std::vector<std::vector<double>>> xi_microphysics_;
    std::vector<std::vector<std::vector<double>>> xi_pbl_;
    std::vector<std::vector<std::vector<double>>> xi_diffusion_;

    // Previous timestep noise (for AR(1) evolution)
    std::vector<std::vector<std::vector<double>>> xi_micro_prev_;
    std::vector<std::vector<std::vector<double>>> xi_pbl_prev_;
    std::vector<std::vector<std::vector<double>>> xi_diff_prev_;

    uint64_t time_step_counter_ = 0;
};

} // namespace chaos
