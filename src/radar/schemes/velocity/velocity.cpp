#include "velocity.hpp"
#include <iostream>

/*This function initializes the velocity scheme.
Takes in the configuration, the number of rows, the number of columns, 
and the number of vertical levels and initializes the velocity scheme.*/
void VelocityScheme::initialize(const RadarConfig& config, int NR, int NTH, int NZ)
{
    NR_ = NR;
    NTH_ = NTH;
    NZ_ = NZ;

    radar_x_ = config.radar_x;
    radar_y_ = config.radar_y;
    radar_z_ = config.radar_z;

    std::cout << "Initialized velocity radar scheme (" << NR << "x" << NTH << "x" << NZ << ")" << std::endl;
    std::cout << "Radar location: (" << radar_x_ << ", " << radar_y_ << ", " << radar_z_ << ")" << std::endl;
}

/*This function computes the velocity scheme.
Takes in the configuration, the state, and the output and computes the velocity scheme.*/
void VelocityScheme::compute(const RadarConfig& config, const RadarStateView& state, RadarOut& out) 
{
    // Compute basic radial velocity
    compute_radial_velocity(config, state, out);

    // Optional: apply scatterer motion correction
    // This would be controlled by a config flag in a full implementation
    // apply_scatterer_correction(config, state, out);
}

/*This function computes the radial velocity.
Takes in the configuration, the state, and the output and computes the radial velocity.*/
void VelocityScheme::compute_radial_velocity(const RadarConfig& config, const RadarStateView& state, RadarOut& out) 
{
    // If the wind fields are not present, throw an error.
    if (!state.u || !state.v || !state.w) 
    {
        std::cerr << "Warning: Velocity scheme requires wind fields (u, v, w)" << std::endl;
        return;
    }

    // Import grid resolution parameters from equations.cpp
    extern double dr, dz, dtheta;

    
    // Iterate over the rows, columns, and vertical levels and compute the radial velocity.
    for (int i = 0; i < NR_; ++i) 
    {
        for (int j = 0; j < NTH_; ++j) 
        {
            for (int k = 0; k < NZ_; ++k) 
            {
                // Get wind components (cylindrical coordinates)
                float u_cyl = (*state.u)[i][j][k];      // radial wind component
                float v_theta = (*state.v)[i][j][k];    // azimuthal wind component
                float w = (*state.w)[i][j][k];          // vertical wind component

                // Convert grid indices to cylindrical physical coordinates
                double r = static_cast<double>(i) * dr;       // radial distance from center
                double theta = static_cast<double>(j) * dtheta; // azimuthal angle
                double z = static_cast<double>(k) * dz;       // vertical coordinate

                // Convert to Cartesian coordinates for radar geometry
                double x = r * cos(theta);
                double y = r * sin(theta);

                // Convert cylindrical winds to Cartesian winds
                // u_cyl = u_cart * cos(theta) + v_cart * sin(theta)
                // v_theta = -u_cart * sin(theta) + v_cart * cos(theta)
                // Solving for Cartesian components:
                double u_cart = u_cyl * cos(theta) - v_theta * sin(theta);
                double v_cart = u_cyl * sin(theta) + v_theta * cos(theta);

                // Compute line-of-sight unit vector
                double e_r_x, e_r_y, e_r_z, R;
                RadarGeometry::compute_line_of_sight(radar_x_, radar_y_, radar_z_,
                                                   x, y, z, e_r_x, e_r_y, e_r_z, R);

                // Radial velocity: V_r = [u_cart, v_cart, w] · ê_r
                float Vr = static_cast<float>(u_cart * e_r_x + v_cart * e_r_y + w * e_r_z);

                out.Vr[i][j][k] = Vr;
            }
        }
    }
}

/*This function applies the scatterer correction.
Takes in the configuration, the state, and the output and applies the scatterer correction.*/
void VelocityScheme::apply_scatterer_correction(const RadarConfig& config, const RadarStateView& state, RadarOut& out) 
{
    // Optional scatterer motion correction
    // V_r_scatterer = V_r_air - V_terminal · ê_r

    // Import grid resolution parameters from equations.cpp
    extern double dr, dz, dtheta;

    // Iterate over the rows, columns, and vertical levels and apply the scatterer correction.
    for (int i = 0; i < NR_; ++i) 
    {
        // Iterate over the columns and vertical levels and apply the scatterer correction.
        for (int j = 0; j < NTH_; ++j) 
        {
            // Iterate over the vertical levels and apply the scatterer correction.
            for (int k = 0; k < NZ_; ++k) 
            {
                // Estimate fall speed
                float q_rain = state.qr ? (*state.qr)[i][j][k] : 0.0f;
                float q_snow = state.qs ? (*state.qs)[i][j][k] : 0.0f;
                float q_graupel = state.qg ? (*state.qg)[i][j][k] : 0.0f;
                float q_hail = state.qh ? (*state.qh)[i][j][k] : 0.0f;
                float q_ice = state.qi ? (*state.qi)[i][j][k] : 0.0f;

                float v_terminal = estimate_fall_speed(q_rain, q_snow, q_graupel, q_hail, q_ice);

                // Convert grid indices to cylindrical physical coordinates
                double r = static_cast<double>(i) * dr;
                double theta = static_cast<double>(j) * dtheta;
                double z = static_cast<double>(k) * dz;

                // Convert to Cartesian coordinates for radar geometry
                double x = r * cos(theta);
                double y = r * sin(theta);

                // Get unit vector for line-of-sight
                double e_r_x, e_r_y, e_r_z, R;
                RadarGeometry::compute_line_of_sight(radar_x_, radar_y_, radar_z_,
                                                   x, y, z, e_r_x, e_r_y, e_r_z, R);

                // Apply correction: V_r_scatterer = V_r_air - V_terminal * ê_r_z
                // (assuming fall speed is in vertical direction)
                float correction = -v_terminal * static_cast<float>(e_r_z);
                out.Vr[i][j][k] += correction;
            }
        }
    }
}

/*This function estimates the fall speed.
Takes in the rain, snow, graupel, hail, and ice and estimates the fall speed.*/
float VelocityScheme::estimate_fall_speed(float q_rain, float q_snow, float q_graupel, float q_hail, float q_ice) 
{
    // Simplified fall speed estimation
    // Real implementation would use proper relations based on particle size, density, etc.

    float v_fall = 0.0f;

    // If the rain is present, estimate the fall speed.
    // Rain fall speed (simplified)
    if (q_rain > 0.0f) 
    {
        v_fall += 5.0f * std::sqrt(q_rain * 1000.0f);  // Rough approximation
    }

    // If the snow is present, estimate the fall speed.
    if (q_snow > 0.0f) 
    {
        v_fall += 1.0f;  // Constant for snow
    }

    // If the graupel or hail is present, estimate the fall speed.
    if (q_graupel > 0.0f || q_hail > 0.0f) 
    {
        float q_gh = q_graupel + q_hail;
        v_fall += 8.0f * std::sqrt(q_gh * 1000.0f);  // Higher fall speed
    }

    // Ice fall speed
    if (q_ice > 0.0f) 
    {
        v_fall += 0.5f;  // Slow for small ice
    }

    return v_fall; // Return the fall speed.
}
