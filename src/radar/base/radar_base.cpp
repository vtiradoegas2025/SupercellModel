#include "radar_base.hpp"
#include <cmath>
#include <algorithm>

/**
 * @brief Shared radar geometry and sampling utilities implementation
 */

void RadarGeometry::compute_line_of_sight(double radar_x, double radar_y, double radar_z,
                                         double x, double y, double z,
                                         double& e_r_x, double& e_r_y, double& e_r_z,
                                         double& R) 
                                         
{
    // Vector from radar to sample point
    double dx = x - radar_x;
    double dy = y - radar_y;
    double dz = z - radar_z;

    // Distance
    R = std::sqrt(dx*dx + dy*dy + dz*dz);

    // Avoid division by zero
    if (R < 1e-6) {
        e_r_x = 0.0;
        e_r_y = 0.0;
        e_r_z = 1.0;  // Pointing up if at radar location
        return;
    }

    // Unit vector
    e_r_x = dx / R;
    e_r_y = dy / R;
    e_r_z = dz / R;
}

void RadarGeometry::sample_state_point(const RadarStateView& state, int i, int j, int k,
                                      RadarStateView& point_state) 
{
    // For v1: simple point sampling at grid indices
    // Future: implement interpolation or beam-volume averaging

    point_state.NR = 1;
    point_state.NTH = 1;
    point_state.NZ = 1;

    // Bounds checking
    i = std::max(0, std::min(i, state.NR - 1));
    j = std::max(0, std::min(j, state.NTH - 1));
    k = std::max(0, std::min(k, state.NZ - 1));

    // Point sampling is handled directly in the calling code for now
    // This function serves as a placeholder for future interpolation/beam-volume sampling
}