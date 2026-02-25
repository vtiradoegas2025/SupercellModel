/**
 * @file perturbation_field.hpp
 * @brief Declarations for the chaos module.
 *
 * Defines interfaces, data structures, and contracts used by
 * the chaos runtime and scheme implementations.
 * This file is part of the src/chaos subsystem.
 */

#pragma once
#include "random_generator.hpp"
#include "field3d.hpp"
#include <vector>
#include <string>

namespace chaos 
{

/**
 * @brief Generates the white noise perturbation field for the 3D field.
 */
/**
 * @brief Generate white noise perturbation field
 * @param rng Random number generator
 * @param NR Radial grid points
 * @param NTH Azimuthal grid points
 * @param NZ Vertical grid points
 * @param stream_key Stream identifier for reproducibility
 * @param field_name Field identifier for stream separation
 * @return 3D perturbation field with N(0,1) values
 */
std::vector<std::vector<std::vector<double>>> generate_white_noise_3d(
    ChaosRNG& rng,
    size_t NR,
    size_t NTH,
    size_t NZ,
    uint64_t stream_key,
    const std::string& field_name = ""
);

/**
 * @brief Generate white noise perturbation field in contiguous Field3D layout
 */
Field3D generate_white_noise_field3d(
    ChaosRNG& rng,
    int NR,
    int NTH,
    int NZ,
    uint64_t stream_key,
    const std::string& field_name = ""
);

/**
 * @brief Generate white noise perturbation field (2D horizontal slice)
 * @param rng Random number generator
 * @param NR Radial grid points
 * @param NTH Azimuthal grid points
 * @param stream_key Stream identifier for reproducibility
 * @param field_name Field identifier for stream separation
 * @return 2D perturbation field with N(0,1) values
 */
std::vector<std::vector<double>> generate_white_noise_2d(
    ChaosRNG& rng,
    size_t NR,
    size_t NTH,
    uint64_t stream_key,
    const std::string& field_name = ""
);

/**
 * @brief Scale perturbation field by specified amplitude
 * @param field Input/output perturbation field
 * @param sigma Scaling factor (standard deviation)
 */
void scale_perturbation_field(
    std::vector<std::vector<std::vector<double>>>& field,
    double sigma
);

void scale_perturbation_field(
    Field3D& field,
    double sigma
);

/**
 * @brief Scale perturbation field by specified amplitude (2D version)
 * @param field Input/output perturbation field
 * @param sigma Scaling factor (standard deviation)
 */
void scale_perturbation_field(
    std::vector<std::vector<double>>& field,
    double sigma
);

/**
 * @brief Ensure perturbation field has unit variance (renormalization)
 * @param field Input/output perturbation field
 * @return Realized variance of the field
 */
double renormalize_to_unit_variance(
    std::vector<std::vector<std::vector<double>>>& field
);

double renormalize_to_unit_variance(
    Field3D& field
);

/**
 * @brief Compute statistical properties of perturbation field
 * @param field Perturbation field
 * @return Tuple of (mean, variance, min, max)
 */
std::tuple<double, double, double, double> compute_field_statistics(
    const std::vector<std::vector<std::vector<double>>>& field
);

/**
 * @brief Apply physical bounds to prevent unphysical states
 * @param perturbations Perturbation values to apply
 * @param base_values Original field values
 * @param variable_name Variable being perturbed ("qv", "qc", etc.)
 * @param min_value Minimum allowed value (e.g., 0 for moisture)
 */
void apply_physical_bounds(
    std::vector<std::vector<std::vector<double>>>& perturbations,
    const std::vector<std::vector<std::vector<double>>>& base_values,
    const std::string& variable_name,
    double min_value = 0.0
);


/**
 * @brief Evolve perturbation field using AR(1) process
 * @param xi Current perturbation field (modified in-place)
 * @param xi_prev Previous perturbation field
 * @param rho_t Temporal correlation coefficient (exp(-dt/tau_t))
 * @param rng Random number generator for innovation term
 * @param stream_key Stream identifier for reproducibility
 * @param time_step Current time step counter
 */
void evolve_ar1_3d(
    std::vector<std::vector<std::vector<double>>>& xi,
    const std::vector<std::vector<std::vector<double>>>& xi_prev,
    double rho_t,
    ChaosRNG& rng,
    uint64_t stream_key,
    uint64_t time_step
);

void evolve_ar1_3d(
    Field3D& xi,
    const Field3D& xi_prev,
    double rho_t,
    ChaosRNG& rng,
    uint64_t stream_key,
    uint64_t time_step
);

/**
 * @brief Evolve perturbation field using AR(1) process (2D version)
 * @param xi Current perturbation field (modified in-place)
 * @param xi_prev Previous perturbation field
 * @param rho_t Temporal correlation coefficient
 * @param rng Random number generator for innovation term
 * @param stream_key Stream identifier for reproducibility
 * @param time_step Current time step counter
 */
void evolve_ar1_2d(
    std::vector<std::vector<double>>& xi,
    const std::vector<std::vector<double>>& xi_prev,
    double rho_t,
    ChaosRNG& rng,
    uint64_t stream_key,
    uint64_t time_step
);

/**
 * @brief Compute temporal correlation coefficient from decorrelation time
 * @param dt Timestep size
 * @param tau_t Decorrelation time
 * @return rho_t = exp(-dt/tau_t)
 */
double compute_temporal_correlation(double dt, double tau_t);

}
