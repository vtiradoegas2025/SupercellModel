/**
 * @file perturbation_field.cpp
 * @brief Implementation for the chaos module.
 *
 * Provides executable logic for the chaos runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/chaos subsystem.
 */

#include "perturbation_field.hpp"
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream>

namespace chaos 
{

/**
 * @brief Generates the white noise for the 3D field.
 */
std::vector<std::vector<std::vector<double>>> generate_white_noise_3d(
    ChaosRNG& rng,
    size_t NR,
    size_t NTH,
    size_t NZ,
    uint64_t stream_key,
    const std::string& field_name
) 
{
    std::vector<std::vector<std::vector<double>>> field(
        NR, std::vector<std::vector<double>>(
            NTH, std::vector<double>(NZ, 0.0)
        )
    );

    rng.fill_field_normal(field, stream_key, field_name);
    return field;
}

Field3D generate_white_noise_field3d(
    ChaosRNG& rng,
    int NR,
    int NTH,
    int NZ,
    uint64_t stream_key,
    const std::string& field_name
)
{
    Field3D field(NR, NTH, NZ);
    rng.fill_field_normal(field, stream_key, field_name);
    return field;
}

/**
 * @brief Generates the white noise for the 2D field.
 */
std::vector<std::vector<double>> generate_white_noise_2d(
    ChaosRNG& rng,
    size_t NR,
    size_t NTH,
    uint64_t stream_key,
    const std::string& field_name
) 
{
    std::vector<std::vector<double>> field(
        NR, std::vector<double>(NTH, 0.0)
    );

    rng.fill_field_normal(field, stream_key, field_name);
    return field;
}

/**
 * @brief Scales the perturbation field.
 */
void scale_perturbation_field(std::vector<std::vector<std::vector<double>>>& field, double sigma) 
{

    if (field.empty() || field[0].empty() || field[0][0].empty()) {return;}

    const size_t NR = field.size();
    const size_t NTH = field[0].size();
    const size_t NZ = field[0][0].size();

    for (size_t i = 0; i < NR; ++i) 
    {
        for (size_t j = 0; j < NTH; ++j) 
        {
            for (size_t k = 0; k < NZ; ++k) 
            {
                field[i][j][k] *= sigma;
            }
        }
    }
}

void scale_perturbation_field(Field3D& field, double sigma)
{
    if (field.empty())
    {
        return;
    }

    float* values = field.data();
    const std::size_t count = field.size();
    for (std::size_t idx = 0; idx < count; ++idx)
    {
        values[idx] = static_cast<float>(static_cast<double>(values[idx]) * sigma);
    }
}

/**
 * @brief Scales the perturbation field for the 2D field.
 */
void scale_perturbation_field(std::vector<std::vector<double>>& field, double sigma) 
{
    if (field.empty() || field[0].empty()) {return;}

    const size_t NR = field.size();
    const size_t NTH = field[0].size();

    for (size_t i = 0; i < NR; ++i) 
    {
        for (size_t j = 0; j < NTH; ++j) 
        {
            field[i][j] *= sigma;
        }
    }
}

/**
 * @brief Renormalizes the perturbation field to unit variance.
 */
double renormalize_to_unit_variance(std::vector<std::vector<std::vector<double>>>& field) 
{
    if (field.empty() || field[0].empty() || field[0][0].empty()) {return 0.0;}

    const size_t NR = field.size();
    const size_t NTH = field[0].size();
    const size_t NZ = field[0][0].size();
    const size_t total_points = NR * NTH * NZ;

    double sum = 0.0;
    double sum_sq = 0.0;

    for (size_t i = 0; i < NR; ++i) 
    {
        for (size_t j = 0; j < NTH; ++j) 
        {
            for (size_t k = 0; k < NZ; ++k) 
            {
                double val = field[i][j][k];
                sum += val;
                sum_sq += val * val;
            }
        }
    }

    double mean = sum / total_points;
    double variance = (sum_sq / total_points) - (mean * mean);

    if (variance > 1e-12) 
    {
        double scale_factor = 1.0 / std::sqrt(variance);

        for (size_t i = 0; i < NR; ++i) 
        {
            for (size_t j = 0; j < NTH; ++j) 
            {
                for (size_t k = 0; k < NZ; ++k) 
                {
                    field[i][j][k] = (field[i][j][k] - mean) * scale_factor;
                }
            }
        }
    }

    return variance;
}

double renormalize_to_unit_variance(Field3D& field)
{
    if (field.empty())
    {
        return 0.0;
    }

    const size_t total_points = field.size();
    const float* values = field.data();
    double sum = 0.0;
    double sum_sq = 0.0;
    for (size_t idx = 0; idx < total_points; ++idx)
    {
        const double val = static_cast<double>(values[idx]);
        sum += val;
        sum_sq += val * val;
    }

    const double mean = sum / static_cast<double>(total_points);
    const double variance = (sum_sq / static_cast<double>(total_points)) - (mean * mean);
    if (variance <= 1e-12)
    {
        return variance;
    }

    const double scale_factor = 1.0 / std::sqrt(variance);
    float* writable = field.data();
    for (size_t idx = 0; idx < total_points; ++idx)
    {
        const double centered = static_cast<double>(writable[idx]) - mean;
        writable[idx] = static_cast<float>(centered * scale_factor);
    }

    return variance;
}

/**
 * @brief Computes the statistics of the perturbation field.
 */
std::tuple<double, double, double, double> compute_field_statistics(
const std::vector<std::vector<std::vector<double>>>& field
) 
{
    if (field.empty() || field[0].empty() || field[0][0].empty()) {return {0.0, 0.0, 0.0, 0.0};}

    const size_t NR = field.size();
    const size_t NTH = field[0].size();
    const size_t NZ = field[0][0].size();
    const size_t total_points = NR * NTH * NZ;

    double sum = 0.0;
    double sum_sq = 0.0;
    double min_val = field[0][0][0];
    double max_val = field[0][0][0];

    for (size_t i = 0; i < NR; ++i) 
    {
        for (size_t j = 0; j < NTH; ++j) 
        {
            for (size_t k = 0; k < NZ; ++k) 
            {
                double val = field[i][j][k];
                sum += val;
                sum_sq += val * val;
                min_val = std::min(min_val, val);
                max_val = std::max(max_val, val);
            }
        }
    }

    double mean = sum / total_points;
    double variance = (sum_sq / total_points) - (mean * mean);

    return {mean, variance, min_val, max_val};
}

/**
 * @brief Applies the physical bounds to the perturbation field.
 */
void apply_physical_bounds(
    std::vector<std::vector<std::vector<double>>>& perturbations,
    const std::vector<std::vector<std::vector<double>>>& base_values,
    const std::string& variable_name,
    double min_value
) 
{
    if (perturbations.empty() || base_values.empty() ||
        perturbations.size() != base_values.size() ||
        perturbations[0].size() != base_values[0].size() ||
        perturbations[0][0].size() != base_values[0][0].size()) 
        {return;}

    const size_t NR = perturbations.size();
    const size_t NTH = perturbations[0].size();
    const size_t NZ = perturbations[0][0].size();

    bool is_moisture = (variable_name == "qv" || variable_name == "qc" ||
                       variable_name == "qr" || variable_name == "qi" ||
                       variable_name == "qs" || variable_name == "qg");

    for (size_t i = 0; i < NR; ++i) 
    {
        for (size_t j = 0; j < NTH; ++j) 
        {
            for (size_t k = 0; k < NZ; ++k) 
            {
                double base_val = base_values[i][j][k];
                double& pert_val = perturbations[i][j][k];
                double final_val = base_val + pert_val;

                if (is_moisture && final_val < min_value) 
                {
                    pert_val = min_value - base_val;
                }
            }
        }
    }
}

/**
 * @brief Evolves the perturbation field for the 3D field.
 */
void evolve_ar1_3d(
    std::vector<std::vector<std::vector<double>>>& xi,
    const std::vector<std::vector<std::vector<double>>>& xi_prev,
    double rho_t,
    ChaosRNG& rng,
    uint64_t stream_key,
    uint64_t time_step
) 
{
    if (xi.empty() || xi_prev.empty() ||
        xi.size() != xi_prev.size() ||
        xi[0].size() != xi_prev[0].size() ||
        xi[0][0].size() != xi_prev[0][0].size()) {return;}

    const size_t NR = xi.size();
    const size_t NTH = xi[0].size();
    const size_t NZ = xi[0][0].size();

    auto innovation = generate_white_noise_3d(rng, NR, NTH, NZ,
                                            stream_key + time_step, "innovation");

    double innovation_scale = std::sqrt(std::max(0.0, 1.0 - rho_t * rho_t));

    for (size_t i = 0; i < NR; ++i) 
    {
        for (size_t j = 0; j < NTH; ++j) 
        {
            for (size_t k = 0; k < NZ; ++k) 
            {
                xi[i][j][k] = rho_t * xi_prev[i][j][k] +
                             innovation_scale * innovation[i][j][k];
            }
        }
    }
}

void evolve_ar1_3d(
    Field3D& xi,
    const Field3D& xi_prev,
    double rho_t,
    ChaosRNG& rng,
    uint64_t stream_key,
    uint64_t time_step
)
{
    if (xi.empty() || xi_prev.empty() ||
        xi.size_r() != xi_prev.size_r() ||
        xi.size_th() != xi_prev.size_th() ||
        xi.size_z() != xi_prev.size_z())
    {
        return;
    }

    const size_t total_points = xi.size();
    const auto innovation = rng.normal_stream(total_points, stream_key + time_step);
    const double innovation_scale = std::sqrt(std::max(0.0, 1.0 - rho_t * rho_t));

    float* current = xi.data();
    const float* previous = xi_prev.data();
    for (size_t idx = 0; idx < total_points; ++idx)
    {
        current[idx] = static_cast<float>(
            rho_t * static_cast<double>(previous[idx]) +
            innovation_scale * innovation[idx]
        );
    }
}

/**
 * @brief Evolves the perturbation field for the 2D field.
 */
void evolve_ar1_2d(
    std::vector<std::vector<double>>& xi,
    const std::vector<std::vector<double>>& xi_prev,
    double rho_t,
    ChaosRNG& rng,
    uint64_t stream_key,
    uint64_t time_step
) 

{
    if (xi.empty() || xi_prev.empty() ||
        xi.size() != xi_prev.size() ||
        xi[0].size() != xi_prev[0].size()) {return;}

    const size_t NR = xi.size();
    const size_t NTH = xi[0].size();

    auto innovation = generate_white_noise_2d(rng, NR, NTH,
                                            stream_key + time_step, "innovation");

    double innovation_scale = std::sqrt(std::max(0.0, 1.0 - rho_t * rho_t));

    for (size_t i = 0; i < NR; ++i) 
    {
        for (size_t j = 0; j < NTH; ++j) 
        {
            xi[i][j] = rho_t * xi_prev[i][j] + innovation_scale * innovation[i][j];
        }
    }
}

/**
 * @brief Computes the temporal correlation.
 */
double compute_temporal_correlation(double dt, double tau_t) 
{
    if (tau_t <= 0.0) 
    {
        return 0.0;
    }
    return std::exp(-dt / tau_t);
}

}
