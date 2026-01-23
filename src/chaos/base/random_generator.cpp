#include "random_generator.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

/*This file contains the implementation of the random number generator.
This file contains the implementation of the normal distribution,
the uniform distribution, and the stream generator.*/
namespace chaos 
{


ChaosRNG::ChaosRNG(uint64_t base_seed, int member_id)
    : base_seed_(base_seed), member_id_(member_id),
      normal_dist_(0.0, 1.0), uniform_dist_(0.0, 1.0)
{
    base_generator_.seed(base_seed_ + static_cast<uint64_t>(member_id_));
}

/*This function generates the normal distribution for the stream.
Takes in the number of points and the stream key
and generates the normal distribution for the stream.*/
std::vector<double> ChaosRNG::normal_stream(size_t count, uint64_t stream_key) 
{
    std::vector<double> result(count);

    // Create deterministic generator for this stream
    uint64_t stream_seed = make_stream_seed(stream_key, 0);
    auto generator = create_stream_generator(stream_seed);

    // Iterate over the number of points
    for (size_t i = 0; i < count; ++i) 
    {
        result[i] = normal_dist_(generator); // Generate the normal distribution for the stream.
    }

    return result;
}

/*This function generates the normal distribution for the stream.
Takes in the stream key and generates the normal distribution for the stream.*/
double ChaosRNG::normal(uint64_t stream_key) 
{
    uint64_t stream_seed = make_stream_seed(stream_key, 0);
    auto generator = create_stream_generator(stream_seed);
    return normal_dist_(generator);
}

/*This function generates the uniform distribution for the stream.
Takes in the number of points and the stream key
and generates the uniform distribution for the stream.*/
std::vector<double> ChaosRNG::uniform_stream(size_t count, uint64_t stream_key) 
{
    std::vector<double> result(count);

    uint64_t stream_seed = make_stream_seed(stream_key, 0);
    auto generator = create_stream_generator(stream_seed);

    // Iterate over the number of points
    for (size_t i = 0; i < count; ++i) 
    {
        result[i] = uniform_dist_(generator); // Generate the uniform distribution for the stream.
    }

    return result;
}

/*This function generates the uniform distribution for the stream.
Takes in the stream key and generates the uniform distribution for the stream.*/
double ChaosRNG::uniform(uint64_t stream_key) 
{
    uint64_t stream_seed = make_stream_seed(stream_key, 0);
    auto generator = create_stream_generator(stream_seed);
    return uniform_dist_(generator);
}

/*This function fills the field with the normal distribution.
Takes in the field, the stream key, and the field name
and fills the field with the normal distribution.*/
void ChaosRNG::fill_field_normal(
    std::vector<std::vector<std::vector<double>>>& field,
    uint64_t stream_key,
    const std::string& field_name
) 
{
    if (field.empty() || field[0].empty() || field[0][0].empty()) {return;}

    // Get the number of rows, columns, and levels.
    const size_t NR = field.size();
    const size_t NTH = field[0].size();
    const size_t NZ = field[0][0].size();
    const size_t total_points = NR * NTH * NZ;

    // Generate all random numbers for this field at once
    uint64_t field_key = stream_key;
    if (!field_name.empty()) 
    {
        // Hash field name into key for additional separation
        for (char c : field_name) 
        {
            field_key = field_key * 31 + static_cast<uint64_t>(c); // Hash the field name into key for additional separation.
        }
    }

    auto random_values = normal_stream(total_points, field_key);

    // Iterate over the rows
    size_t idx = 0;
    for (size_t i = 0; i < NR; ++i) 
    {
        // Iterate over the columns
        for (size_t j = 0; j < NTH; ++j) 
        {
            // Iterate over the levels
            for (size_t k = 0; k < NZ; ++k) 
            {
                field[i][j][k] = random_values[idx++]; // Fill the field with the normal distribution.
            }
        }
    }
}

/*This function fills the field with the normal distribution.
Takes in the field, the stream key, and the field name
and fills the field with the normal distribution.*/
void ChaosRNG::fill_field_normal(
    std::vector<std::vector<double>>& field,
    uint64_t stream_key,
    const std::string& field_name
) 
{
    if (field.empty() || field[0].empty()) {return;}

    const size_t NR = field.size();
    const size_t NTH = field[0].size();
    const size_t total_points = NR * NTH;

    uint64_t field_key = stream_key;

    // If the field name is not empty, hash the field name into key for additional separation.
    if (!field_name.empty()) 
    {
        for (char c : field_name) 
        {
            field_key = field_key * 31 + static_cast<uint64_t>(c);
        }
    }

    // Generate all random numbers for this field at once
    auto random_values = normal_stream(total_points, field_key);

    size_t idx = 0;

    // Iterate over the rows
    for (size_t i = 0; i < NR; ++i) 
    {
        // Iterate over the columns
        for (size_t j = 0; j < NTH; ++j) 
        {
            field[i][j] = random_values[idx++]; // Fill the field with the normal distribution.
        }
    }
}

/*This function resets the random number generator.
Takes in the base seed and the member id
and resets the random number generator.*/
void ChaosRNG::reset() 
{
    base_generator_.seed(base_seed_ + static_cast<uint64_t>(member_id_));
}

/*This function makes the stream seed.
Takes in the stream key and the sub key
and makes the stream seed.*/
uint64_t ChaosRNG::make_stream_seed(uint64_t stream_key, uint64_t sub_key) const 
{
    // Combine base seed, member ID, and stream key deterministically
    // Use a simple hash-like combination to avoid collisions
    uint64_t combined = base_seed_;
    combined = combined * 6364136223846793005ULL + static_cast<uint64_t>(member_id_);
    combined = combined * 6364136223846793005ULL + stream_key;
    combined = combined * 6364136223846793005ULL + sub_key;
    return combined;
}

/*This function creates the stream generator.
Takes in the stream seed and creates the stream generator.*/
std::mt19937_64 ChaosRNG::create_stream_generator(uint64_t stream_seed) const 
{
    std::mt19937_64 generator;
    generator.seed(stream_seed);
    return generator;
}

/*This function bounds the perturbation field.
Takes in the perturbation field, the maximum value, and the use of tanh
and bounds the perturbation field.*/

void bound_perturbation_field(
    std::vector<std::vector<std::vector<double>>>& xi,
    double xi_max,
    bool use_tanh
) 
{
    // If the perturbation field is empty, return.
    if (xi.empty() || xi[0].empty() || xi[0][0].empty()) {return;}

    // Get the number of rows, columns, and levels.
    const size_t NR = xi.size();
    const size_t NTH = xi[0].size();
    const size_t NZ = xi[0][0].size();

    // Iterate over the rows
    for (size_t i = 0; i < NR; ++i) 
    {
        // Iterate over the columns
        for (size_t j = 0; j < NTH; ++j) 
        {
            // Iterate over the levels
            for (size_t k = 0; k < NZ; ++k) 
            {
                double& val = xi[i][j][k];

                // If the use of tanh is true, smooth saturation using tanh.
                if (use_tanh) 
                {
                    // Smooth saturation using tanh
                    val = xi_max * std::tanh(val / xi_max);
                } 
                else 
                {
                    // Hard clipping
                    val = std::max(-xi_max, std::min(xi_max, val));
                }
            }
        }
    }
}

/*This function applies the vertical taper to the perturbation field.
Takes in the perturbation field, the z levels, the taper id, the z1, and the z2
and applies the vertical taper to the perturbation field.*/
void apply_vertical_taper(
    std::vector<std::vector<std::vector<double>>>& xi,
    const std::vector<double>& z_levels,
    const std::string& taper_id,
    double z1,
    double z2
) 
{
    // If the perturbation field, z levels, or taper id is empty, return.
    if (xi.empty() || xi[0].empty() || xi[0][0].empty() || z_levels.empty()) {return;}

    // Get the number of rows, columns, and levels.
    const size_t NR = xi.size();
    const size_t NTH = xi[0].size();
    const size_t NZ = xi[0][0].size();

    // If the taper id is none, return.
    if (taper_id == "none") 
    {
        return;  // No tapering
    }

    // Iterate over the rows
    for (size_t i = 0; i < NR; ++i) 
    {
        // Iterate over the columns
        for (size_t j = 0; j < NTH; ++j) 
        {
            // Iterate over the levels
            for (size_t k = 0; k < NZ; ++k) 
            {
                double z = z_levels[k];
                double taper_factor = 1.0;

                // If the taper id is pbl_only, apply the perturbation below PBL height.
                if (taper_id == "pbl_only") 
                {
                    // Only apply perturbations below PBL height
                    taper_factor = (z <= z1) ? 1.0 : 0.0;
                } 
                else if (taper_id == "cosine") 
                {
                    // if the z is less than or equal to z1, set the taper factor to 1.
                    if (z <= z1) 
                    {
                        taper_factor = 1.0;
                    }
                     else if (z >= z2) 
                    {
                        taper_factor = 0.0;
                    }
                     else 
                    {
                        // Cosine taper: 1 at z1, 0 at z2
                        double phase = M_PI * (z - z1) / (z2 - z1);
                        taper_factor = 0.5 * (1.0 + std::cos(phase));
                    }
                }
                xi[i][j][k] *= taper_factor; // Apply the vertical taper to the perturbation field.
            }
        }
    }
}

/*This function applies the horizontal mask to the perturbation field.
Takes in the perturbation field and the mask
and applies the horizontal mask to the perturbation field.*/
void apply_horizontal_mask(
    std::vector<std::vector<std::vector<double>>>& xi,
    const std::vector<std::vector<double>>& mask
) 
{
    // If the perturbation field, mask, or mask is empty, return.
    if (xi.empty() || xi[0].empty() || xi[0][0].empty() ||
        mask.empty() || mask[0].empty()) {return;}

    const size_t NR = xi.size();
    const size_t NTH = xi[0].size();
    const size_t NZ = xi[0][0].size();

    // If the mask dimensions do not match the perturbation field dimensions, return.
    if (mask.size() != NR || mask[0].size() != NTH) 
    {
        std::cerr << "Warning: Mask dimensions don't match perturbation field" << std::endl;
        return;
    }

    // Iterate over the rows
    for (size_t i = 0; i < NR; ++i) 
    {
        // Iterate over the columns
        for (size_t j = 0; j < NTH; ++j) 
        {
            double mask_val = mask[i][j];

            // Iterate over the levels
            for (size_t k = 0; k < NZ; ++k) 
            {
                xi[i][j][k] *= mask_val; // Apply the horizontal mask to the perturbation field.
            }
        }
    }
}

} // namespace chaos
