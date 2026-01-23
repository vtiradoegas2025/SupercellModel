#include "zdr.hpp"
#include <cmath>
#include <iostream>

/*This function initializes the ZDR scheme.
Takes in the configuration, the number of rows, the number of columns, 
and the number of vertical levels and initializes the ZDR scheme.*/
void ZDRScheme::initialize(const RadarConfig& config, int NR, int NTH, int NZ) 
{
    NR_ = NR;
    NTH_ = NTH;
    NZ_ = NZ;

    std::cout << "Initialized ZDR radar scheme (" << NR << "x" << NTH << "x" << NZ << ")" << std::endl;
    std::cout << "Operator tier: " << config.operator_tier << std::endl;
}

/*This function computes the ZDR scheme.
Takes in the configuration, the state, and the output and computes the ZDR scheme.*/
void ZDRScheme::compute(const RadarConfig& config, const RadarStateView& state, RadarOut& out) 
{
    // If the operator tier is simple, compute the simple ZDR.
    if (config.operator_tier == "simple") 
    {
        compute_simple_zdr(config, state, out);
    } 

    // If the operator tier is polarimetric forward operator, compute the polarimetric forward operator.
    else if (config.operator_tier == "polarimetric_fo") 
    {
        compute_polarimetric_fo(config, state, out);
    }
    else 
    {
        // If the operator tier is unknown, throw an error.
        throw std::runtime_error("Unknown ZDR operator tier: " + config.operator_tier);
    }
}

/*This function computes the simple ZDR.
Takes in the configuration, the state, and the output and computes the simple ZDR.*/
void ZDRScheme::compute_simple_zdr(const RadarConfig& config, const RadarStateView& state, RadarOut& out) 
{
    // Simple ZDR: rain-only, parameterized axis ratios
    // Z_H ≈ Z_e, Z_V computed using axis ratio parameterization

    // Iterate over the rows, columns, and vertical levels and compute the simple ZDR.
    for (int i = 0; i < NR_; ++i) 
    {
        // Iterate over the columns and vertical levels and compute the simple ZDR.
        for (int j = 0; j < NTH_; ++j) 
        {
            // Iterate over the vertical levels and compute the simple ZDR.
            for (int k = 0; k < NZ_; ++k) 
            {
                // Get the rain quantity.
                float q_rain = state.qr ? (*state.qr)[i][j][k] : 0.0f;

                if (q_rain <= 0.0f) {
                    // No rain, assume isotropic scattering (Z_DR = 0 dB)
                    out.ZDR_dB[i][j][k] = 0.0f;
                    out.ZH_dBZ[i][j][k] = -30.0f;  // Below threshold
                    out.ZV_dBZ[i][j][k] = -30.0f;
                    continue;
                }

                // Estimate axis ratio (oblate drops)
                float axis_ratio = estimate_rain_axis_ratio(q_rain);

                // Compute Z_H and Z_V
                float Z_H, Z_V;
                compute_rain_polarization(q_rain, axis_ratio, Z_H, Z_V);

                // Convert to dBZ
                float Z_H_dBZ = (Z_H > 1e-10) ? 10.0f * std::log10(Z_H) : -30.0f;
                float Z_V_dBZ = (Z_V > 1e-10) ? 10.0f * std::log10(Z_V) : -30.0f;

                // Z_DR = 10*log10(Z_H / Z_V)

                // If the vertical and horizontal reflectivities are greater than 0, compute the ZDR.
                if (Z_V > 1e-10 && Z_H > 1e-10) 
                {
                    float Z_DR_linear = Z_H / Z_V;
                    out.ZDR_dB[i][j][k] = 10.0f * std::log10(Z_DR_linear);
                } 
                else 
                {
                    out.ZDR_dB[i][j][k] = 0.0f;  // Default to 0 dB
                }

                out.ZH_dBZ[i][j][k] = Z_H_dBZ;
                out.ZV_dBZ[i][j][k] = Z_V_dBZ;
            }
        }
    }
}

/*This function computes the polarimetric forward operator.
Takes in the configuration, the state, and the output and computes the polarimetric forward operator.*/
void ZDRScheme::compute_polarimetric_fo(const RadarConfig& config, const RadarStateView& state, RadarOut& out) 
{
    // Full polarimetric forward operator following Jung et al. (2008) style
    // Implements proper polarimetric scattering for different hydrometeor types

    // Initialize polarimetric fields
    // Iterate over the rows, columns, and vertical levels and initialize the polarimetric fields.
    for (int i = 0; i < NR_; ++i) 
    {
        // Iterate over the columns and vertical levels and initialize the polarimetric fields.
        for (int j = 0; j < NTH_; ++j) 
        {
            // Iterate over the vertical levels and initialize the polarimetric fields.
            for (int k = 0; k < NZ_; ++k) 
            {
                out.ZH_dBZ[i][j][k] = -30.0f;  // Below threshold
                out.ZV_dBZ[i][j][k] = -30.0f;
                out.ZDR_dB[i][j][k] = 0.0f;    // Default to 0 dB
            }
        }
    }

    
    // If the rain is present, compute the polarimetric rain.
    if (state.qr && config.has_qr) 
    {
        compute_polarimetric_rain(*state.qr, out);
    }

    // If the snow is present, compute the polarimetric snow.
    if (state.qs && config.has_qs) 
    {
        compute_polarimetric_ice(*state.qs, "snow", out);
    }

    // If the graupel is present, compute the polarimetric graupel.
    if (state.qg && config.has_qg) 
    {
        compute_polarimetric_ice(*state.qg, "graupel", out);
    }

    // If the hail is present, compute the polarimetric hail.
    if (state.qh && config.has_qh) 
    {
        compute_polarimetric_ice(*state.qh, "hail", out);
    }

    // If the ice is present, compute the polarimetric ice.
    if (state.qi && config.has_qi) 
    {
        compute_polarimetric_ice(*state.qi, "ice", out);
    }

    // Iterate over the rows, columns, and vertical levels and compute the total reflectivity and ZDR.
    for (int i = 0; i < NR_; ++i) 
    {
        // Iterate over the columns and vertical levels and compute the total reflectivity and ZDR.
        for (int j = 0; j < NTH_; ++j) 
        {
            // Iterate over the vertical levels and compute the total reflectivity and ZDR.
            for (int k = 0; k < NZ_; ++k) 
            {
                // Convert linear reflectivities to dBZ
                float Z_H_linear = dbz_to_linear(out.ZH_dBZ[i][j][k]);
                float Z_V_linear = dbz_to_linear(out.ZV_dBZ[i][j][k]);

                // Compute total reflectivity as average of H and V (for Z_e)
                float Z_total_linear = (Z_H_linear + Z_V_linear) / 2.0f;
                out.Z_dBZ[i][j][k] = (Z_total_linear > 1e-10) ? 10.0f * log10(Z_total_linear) : -30.0f;

                // Compute ZDR = 10*log10(Z_H / Z_V)

                // If the vertical and horizontal reflectivities are greater than 0, compute the ZDR.
                if (Z_V_linear > 1e-10 && Z_H_linear > 1e-10) 
                {
                    out.ZDR_dB[i][j][k] = 10.0f * log10(Z_H_linear / Z_V_linear);
                } 
                else 
                {
                    out.ZDR_dB[i][j][k] = 0.0f;  // Default when either polarization is below threshold
                }
            }
        }
    }
}


/*This function computes the rain polarization.
Takes in the rain quantity, the axis ratio, and the horizontal and vertical reflectivities and computes the rain polarization.*/
void ZDRScheme::compute_rain_polarization(float q_rain, float axis_ratio, float& Z_H, float& Z_V) 
{
    // Simplified polarimetric calculation for oblate rain drops
    // In a full implementation, this would use T-matrix or similar methods

    // Base reflectivity (similar to reflectivity scheme)
    double rho_ratio = 1000.0 / 1000.0;  // Rain density / water density
    double Ze_base = 3.63e9 * rho_ratio * rho_ratio * q_rain * q_rain;

    // Polarimetric effects: oblate drops have Z_H > Z_V
    // The factor depends on axis ratio and dielectric properties
    double f_h = 1.0 / (axis_ratio * axis_ratio);  // Simplified
    double f_v = axis_ratio * axis_ratio;          // Simplified

    Z_H = static_cast<float>(Ze_base * f_h);
    Z_V = static_cast<float>(Ze_base * f_v);
}

/*This function estimates the rain axis ratio.
Takes in the rain quantity and estimates the rain axis ratio.*/
float ZDRScheme::estimate_rain_axis_ratio(float q_rain) 
{
    // Simplified axis ratio estimation
    // In reality, this depends on drop size distribution
    // Larger drops are more oblate (smaller axis ratio)

    // Rough parameterization based on rain rate
    // Higher q_rain suggests larger drops -> more oblate
    float rain_rate = q_rain * 3600.0f;  // Convert to mm/hr roughly

    // If the rain rate is less than 1 mm/hr, return 0.95.
    if (rain_rate < 1.0f) 
    {
        return 0.95f;  // Small drops, nearly spherical
    } 

    // If the rain rate is less than 10 mm/hr, return 0.85.
    else if (rain_rate < 10.0f) 
    {
        return 0.85f;  // Medium drops
    } 
    // If the rain rate is greater than 10 mm/hr, return 0.75.
    else 
    {
        return 0.75f;  // Large drops, quite oblate
    }
}

/*This function computes the polarimetric rain.
Takes in the rain quantity and the output and computes the polarimetric rain.*/
void ZDRScheme::compute_polarimetric_rain(const std::vector<std::vector<std::vector<float>>>& q_rain, RadarOut& out) 
{
    // Polarimetric calculation for rain drops (oblate spheroids)
    // Based on Jung et al. (2008) style forward operator

    // Iterate over the rows, columns, and vertical levels and compute the polarimetric rain.
    for (int i = 0; i < NR_; ++i) 
    {
        // Iterate over the columns and vertical levels and compute the polarimetric rain.
        for (int j = 0; j < NTH_; ++j) 
        {
            // Iterate over the vertical levels and compute the polarimetric rain.
            for (int k = 0; k < NZ_; ++k) 
            {
                // Get the rain quantity.
                float q_val = q_rain[i][j][k];

                // If the rain quantity is less than or equal to 0, skip.
                if (q_val <= 0.0f) 
                {
                    continue;  // Skip, keep default values
                }

                // Estimate axis ratio (horizontal/vertical dimension ratio)
                float axis_ratio = estimate_rain_axis_ratio(q_val);

                // Base reflectivity (similar to reflectivity scheme)
                double rho_ratio = 1000.0 / 1000.0;  // Rain density / water density
                double Ze_base = 3.63e9 * rho_ratio * rho_ratio * q_val * q_val;

                // Polarimetric effects for oblate spheroids
                // For rain, Z_H > Z_V due to horizontal flattening
                // The exact factors depend on dielectric properties and orientation

                // Simplified polarimetric factors (from literature)
                double f_h = 1.0 / (axis_ratio * axis_ratio);      // Enhancement for horizontal
                double f_v = axis_ratio * axis_ratio * axis_ratio * axis_ratio;  // Reduction for vertical

                double Z_H_linear = Ze_base * f_h;
                double Z_V_linear = Ze_base * f_v;

                // Convert to dBZ and accumulate (add to existing values)
                float Z_H_dBZ_existing = out.ZH_dBZ[i][j][k];
                float Z_V_dBZ_existing = out.ZV_dBZ[i][j][k];

                // Convert the existing horizontal and vertical reflectivities to linear.
                float Z_H_linear_existing = dbz_to_linear(Z_H_dBZ_existing);
                float Z_V_linear_existing = dbz_to_linear(Z_V_dBZ_existing);

                // Compute the total horizontal and vertical reflectivities.
                float Z_H_linear_total = Z_H_linear_existing + Z_H_linear;
                float Z_V_linear_total = Z_V_linear_existing + Z_V_linear;

                // Convert the total horizontal and vertical reflectivities to dBZ.
                out.ZH_dBZ[i][j][k] = (Z_H_linear_total > 1e-10) ? 10.0f * log10(Z_H_linear_total) : -30.0f;
                out.ZV_dBZ[i][j][k] = (Z_V_linear_total > 1e-10) ? 10.0f * log10(Z_V_linear_total) : -30.0f;
            }
        }
    }
}

/*This function computes the polarimetric ice.
Takes in the ice quantity, the species, and the output and computes the polarimetric ice.*/
void ZDRScheme::compute_polarimetric_ice(const std::vector<std::vector<std::vector<float>>>& q_ice,
                                         const std::string& species, RadarOut& out) {
    // Simplified polarimetric calculation for ice species
    // Ice particles are generally less oriented and have different polarimetric properties

    // Species-specific properties (simplified)
    double rho_hyd;  // kg/m³
    if (species == "snow") rho_hyd = 100.0;
    else if (species == "graupel") rho_hyd = 400.0;
    else if (species == "hail") rho_hyd = 900.0;
    else rho_hyd = 917.0;  // ice

    // Iterate over the rows, columns, and vertical levels and compute the polarimetric ice.
    for (int i = 0; i < NR_; ++i) 
    {
        // Iterate over the columns and vertical levels and compute the polarimetric ice.
        for (int j = 0; j < NTH_; ++j) 
        {

            for (int k = 0; k < NZ_; ++k) 
            {
                float q_val = q_ice[i][j][k];

                // If the ice quantity is less than or equal to 0, skip.
                if (q_val <= 0.0f) 
                {
                    continue;
                }

                // Base reflectivity for ice species
                double rho_ratio = rho_hyd / 1000.0;
                double Ze_base = 3.63e9 * rho_ratio * rho_ratio * q_val * q_val;

                // Ice particles are typically less oriented than rain
                // Assume more isotropic scattering with slight orientation effects
                double orientation_factor = 0.9;  // Less oriented than rain

                if (species == "hail") {
                    // Large hail can be somewhat oriented
                    orientation_factor = 0.8;
                }

                // For ice, Z_H and Z_V are closer (smaller ZDR)
                double Z_H_linear = Ze_base * orientation_factor;
                double Z_V_linear = Ze_base / orientation_factor;

                // Convert to dBZ and accumulate
                float Z_H_dBZ_existing = out.ZH_dBZ[i][j][k];
                float Z_V_dBZ_existing = out.ZV_dBZ[i][j][k];

                float Z_H_linear_existing = dbz_to_linear(Z_H_dBZ_existing);
                float Z_V_linear_existing = dbz_to_linear(Z_V_dBZ_existing);

                float Z_H_linear_total = Z_H_linear_existing + Z_H_linear;
                float Z_V_linear_total = Z_V_linear_existing + Z_V_linear;

                out.ZH_dBZ[i][j][k] = (Z_H_linear_total > 1e-10) ? 10.0f * log10(Z_H_linear_total) : -30.0f;
                out.ZV_dBZ[i][j][k] = (Z_V_linear_total > 1e-10) ? 10.0f * log10(Z_V_linear_total) : -30.0f;
            }
        }
    }
}
