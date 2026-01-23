#include "perturbation_field.hpp"
#include <cmath>
#include <numeric>
#include <algorithm>
#include <iostream>

/*This file contains the implementation of the perturbation field.
This file contains the implementation of the white noise generation,
scaling and renormalization, physical bounds enforcement,
and temporal evolution.*/
namespace chaos 
{

/*This function generates the white noise for the 3D field.
Takes in the random number generator, the number of rows, the number of columns,
the number of levels, the stream key, and the field name
and generates the white noise for the 3D field.*/
std::vector<std::vector<std::vector<double>>> generate_white_noise_3d(
    ChaosRNG& rng,
    size_t NR,
    size_t NTH,
    size_t NZ,
    uint64_t stream_key,
    const std::string& field_name
) 
{
    // Initialize the 3D field.
    std::vector<std::vector<std::vector<double>>> field(
        NR, std::vector<std::vector<double>>(
            NTH, std::vector<double>(NZ, 0.0)
        )
    );

    // Fill the 3D field with normal random numbers.
    rng.fill_field_normal(field, stream_key, field_name);
    return field;
}

/*This function generates the white noise for the 2D field.
Takes in the random number generator, the number of rows, the number of columns,
the number of levels, the stream key, and the field name
and generates the white noise for the 2D field.*/
std::vector<std::vector<double>> generate_white_noise_2d(
    ChaosRNG& rng,
    size_t NR,
    size_t NTH,
    uint64_t stream_key,
    const std::string& field_name
) 
{
    // Initialize the 2D field.
    std::vector<std::vector<double>> field(
        NR, std::vector<double>(NTH, 0.0)
    );

    rng.fill_field_normal(field, stream_key, field_name);
    return field;
}

/*This function scales the perturbation field.
Takes in the field, and the scaling factor
and scales the perturbation field.*/
void scale_perturbation_field(std::vector<std::vector<std::vector<double>>>& field, double sigma) 
{

    // If the field is empty, return.
    if (field.empty() || field[0].empty() || field[0][0].empty()) {return;}

    // Get the number of rows, columns, and levels.
    const size_t NR = field.size();
    const size_t NTH = field[0].size();
    const size_t NZ = field[0][0].size();

    // Iterate over the rows
    for (size_t i = 0; i < NR; ++i) 
    {
        // Iterate over the columns
        for (size_t j = 0; j < NTH; ++j) 
        {
            // Iterate over the levels
            for (size_t k = 0; k < NZ; ++k) 
            {
                field[i][j][k] *= sigma; // Scale the field.
            }
        }
    }
}

/*This function scales the perturbation field for the 2D field.
Takes in the field, and the scaling factor
and scales the perturbation field for the 2D field.*/
void scale_perturbation_field(std::vector<std::vector<double>>& field, double sigma) 
{
    // If the field is empty, return.
    if (field.empty() || field[0].empty()) {return;}

    const size_t NR = field.size();
    const size_t NTH = field[0].size();

    // Iterate over the rows
    for (size_t i = 0; i < NR; ++i) 
    {
        // Iterate over the columns
        for (size_t j = 0; j < NTH; ++j) 
        {
            field[i][j] *= sigma; // Scale the field.
        }
    }
}

/*This function renormalizes the perturbation field to unit variance.
Takes in the field and renormalizes the perturbation field to unit variance.*/
double renormalize_to_unit_variance(std::vector<std::vector<std::vector<double>>>& field) 
{
    // If the field is empty, return.
    if (field.empty() || field[0].empty() || field[0][0].empty()) {return 0.0;}

    // Get the number of rows, columns, and levels.
    const size_t NR = field.size();
    const size_t NTH = field[0].size();
    const size_t NZ = field[0][0].size();
    const size_t total_points = NR * NTH * NZ;

    // Compute current variance
    double sum = 0.0;
    double sum_sq = 0.0;

    // Iterate over the rows
    for (size_t i = 0; i < NR; ++i) 
    {
        // Iterate over the columns
        for (size_t j = 0; j < NTH; ++j) 
        {
            // Iterate over the levels
            for (size_t k = 0; k < NZ; ++k) 
            {
                double val = field[i][j][k]; // Get the value of the field.
                sum += val; // Add the value to the sum.
                sum_sq += val * val; // Add the value squared to the sum of squares.
            }
        }
    }

    double mean = sum / total_points;
    double variance = (sum_sq / total_points) - (mean * mean);

    // If the variance is greater than 1e-12, renormalize the field.
    if (variance > 1e-12) 
    {
        double scale_factor = 1.0 / std::sqrt(variance);

        // Iterate over the rows 
        for (size_t i = 0; i < NR; ++i) 
        {
            // Iterate over the columns
            for (size_t j = 0; j < NTH; ++j) 
            {
                // Iterate over the levels
                for (size_t k = 0; k < NZ; ++k) 
                {
                    field[i][j][k] = (field[i][j][k] - mean) * scale_factor; // Renormalize the field.
                }
            }
        }
    }

    return variance;
}

/*This function computes the statistics of the perturbation field.
Takes in the field and computes the statistics of the perturbation field.*/
std::tuple<double, double, double, double> compute_field_statistics(
const std::vector<std::vector<std::vector<double>>>& field
) 
{
    // If the field is empty, return.
    if (field.empty() || field[0].empty() || field[0][0].empty()) {return {0.0, 0.0, 0.0, 0.0};}

    // Get the number of rows, columns, and levels.
    const size_t NR = field.size();
    const size_t NTH = field[0].size();
    const size_t NZ = field[0][0].size();
    const size_t total_points = NR * NTH * NZ;

    double sum = 0.0;
    double sum_sq = 0.0;
    double min_val = field[0][0][0];
    double max_val = field[0][0][0];

    // Iterate over the rows
    for (size_t i = 0; i < NR; ++i) 
    {
        // Iterate over the columns
        for (size_t j = 0; j < NTH; ++j) 
        {
            // Iterate over the levels
            for (size_t k = 0; k < NZ; ++k) 
            {
                double val = field[i][j][k]; // Get the value of the field.
                sum += val; // Add the value to the sum.
                sum_sq += val * val; // Add the value squared to the sum of squares.
                min_val = std::min(min_val, val); // Update the minimum value.
                max_val = std::max(max_val, val); // Update the maximum value.
            }
        }
    }

    double mean = sum / total_points;
    double variance = (sum_sq / total_points) - (mean * mean);

    return {mean, variance, min_val, max_val};
}

/*This function applies the physical bounds to the perturbation field.
Takes in the perturbation field, the base values, the variable name, and the minimum value
and applies the physical bounds to the perturbation field.*/
void apply_physical_bounds(
    std::vector<std::vector<std::vector<double>>>& perturbations,
    const std::vector<std::vector<std::vector<double>>>& base_values,
    const std::string& variable_name,
    double min_value
) 
{
    // If the perturbation field, base values, or variable name is empty, return.
    if (perturbations.empty() || base_values.empty() ||
        perturbations.size() != base_values.size() ||
        perturbations[0].size() != base_values[0].size() ||
        perturbations[0][0].size() != base_values[0][0].size()) 
        {return;}

    const size_t NR = perturbations.size();
    const size_t NTH = perturbations[0].size();
    const size_t NZ = perturbations[0][0].size();

    // For moisture variables, ensure non-negative values
    bool is_moisture = (variable_name == "qv" || variable_name == "qc" ||
                       variable_name == "qr" || variable_name == "qi" ||
                       variable_name == "qs" || variable_name == "qg");

    // Iterate over the rows
    for (size_t i = 0; i < NR; ++i) 
    {
        // Iterate over the columns
        for (size_t j = 0; j < NTH; ++j) 
        {
            // Iterate over the levels
            for (size_t k = 0; k < NZ; ++k) 
            {
                double base_val = base_values[i][j][k];
                double& pert_val = perturbations[i][j][k];
                double final_val = base_val + pert_val;

                // If the final value is less than the minimum value, adjust the perturbation to maintain non-negativity.
                if (is_moisture && final_val < min_value) 
                {
                    // Adjust perturbation to maintain non-negativity
                    pert_val = min_value - base_val;
                }
            }
        }
    }
}


/
/*This function evolves the perturbation field for the 3D field.
Takes in the perturbation field, the previous perturbation field, the temporal correlation,
the random number generator, the stream key, and the time step
and evolves the perturbation field for the 3D field.*/
void evolve_ar1_3d(
    std::vector<std::vector<std::vector<double>>>& xi,
    const std::vector<std::vector<std::vector<double>>>& xi_prev,
    double rho_t,
    ChaosRNG& rng,
    uint64_t stream_key,
    uint64_t time_step
) 
{
    // If the perturbation field, previous perturbation field, or time step is empty, return.
    if (xi.empty() || xi_prev.empty() ||
        xi.size() != xi_prev.size() ||
        xi[0].size() != xi_prev[0].size() ||
        xi[0][0].size() != xi_prev[0][0].size()) {return;}

    // Get the number of rows, columns, and levels.
    const size_t NR = xi.size();
    const size_t NTH = xi[0].size();
    const size_t NZ = xi[0][0].size();

    // Generate innovation term (white noise)
    auto innovation = generate_white_noise_3d(rng, NR, NTH, NZ,
                                            stream_key + time_step, "innovation");

    // AR(1) evolution: xi = rho_t * xi_prev + sqrt(1-rho_t^2) * innovation
    double innovation_scale = std::sqrt(1.0 - rho_t * rho_t);

    // Iterate over the rows
    for (size_t i = 0; i < NR; ++i) 
    {
        // Iterate over the columns
        for (size_t j = 0; j < NTH; ++j) 
        {
            // Iterate over the levels
            for (size_t k = 0; k < NZ; ++k) 
            {
                // Evolve the perturbation field.
                xi[i][j][k] = rho_t * xi_prev[i][j][k] +
                             innovation_scale * innovation[i][j][k];
            }
        }
    }
}

/*This function evolves the perturbation field for the 2D field.
Takes in the perturbation field, the previous perturbation field, the temporal correlation,
the random number generator, the stream key, and the time step
and evolves the perturbation field for the 2D field.*/
void evolve_ar1_2d(
    std::vector<std::vector<double>>& xi,
    const std::vector<std::vector<double>>& xi_prev,
    double rho_t,
    ChaosRNG& rng,
    uint64_t stream_key,
    uint64_t time_step
) 

{
    // If the perturbation field, previous perturbation field, or time step is empty, return.
    if (xi.empty() || xi_prev.empty() ||
        xi.size() != xi_prev.size() ||
        xi[0].size() != xi_prev[0].size()) {return;}

    const size_t NR = xi.size();
    const size_t NTH = xi[0].size();

    // Generate innovation term
    auto innovation = generate_white_noise_2d(rng, NR, NTH,
                                            stream_key + time_step, "innovation");

    // Compute the innovation scale.
    double innovation_scale = std::sqrt(1.0 - rho_t * rho_t);

    // Iterate over the rows
    for (size_t i = 0; i < NR; ++i) 
    {
        // Iterate over the columns
        for (size_t j = 0; j < NTH; ++j) 
        {
            // Evolve the perturbation field.
            xi[i][j] = rho_t * xi_prev[i][j] + innovation_scale * innovation[i][j];
        }
    }
}

/*This function computes the temporal correlation.
Takes in the time step and the temporal correlation time
and computes the temporal correlation.*/
double compute_temporal_correlation(double dt, double tau_t) 
{
    // If the temporal correlation time is less than or equal to 0, return 0.
    if (tau_t <= 0.0) 
    {
        return 0.0;  // No correlation
    }
    return std::exp(-dt / tau_t); // Compute the temporal correlation.
}

} // namespace chaos
