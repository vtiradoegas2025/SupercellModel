#include "thermodynamics.hpp"
#include <algorithm>


/* This is the utility to convert 3D field of theta to temperature
// It converts the potential temperature field to the temperature field
// by iterating over all grid points
// and converting the potential temperature to the temperature using the theta_to_temperature function
// The temperature field is then stored in the temperature variable. This is simplified for the purpose of this project
// In a full implementation, this would be a more complex function
// For example, the potential temperature to temperature conversion would be a more complex function*/


/*This function converts the potential temperature field to the temperature field.
Takes in the theta field, the pressure field, and the temperature field
and converts the potential temperature field to the temperature field.*/
void thermodynamics::convert_theta_to_temperature_field(
    const std::vector<std::vector<std::vector<float>>>& theta,
    const std::vector<std::vector<std::vector<float>>>& p,
    std::vector<std::vector<std::vector<float>>>& temperature
) 
{
    // Get the number of rows, columns, and levels.
    int NR = theta.size();
    if (NR == 0) return;
    int NTH = theta[0].size();
    if (NTH == 0) return;
    int NZ = theta[0][0].size();

    temperature.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));

    // Iterate over all rows.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal points
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical points
            for (int k = 0; k < NZ; ++k) 
            {
                temperature[i][j][k] = static_cast<float>(
                    theta_to_temperature(theta[i][j][k], p[i][j][k])
                );
            }
        }
    }
}

/*This function converts the temperature field to the potential temperature field.
Takes in the temperature field, the pressure field, and the potential temperature field
and converts the temperature field to the potential temperature field.*/
void thermodynamics::convert_temperature_to_theta_field(
    const std::vector<std::vector<std::vector<float>>>& temperature,
    const std::vector<std::vector<std::vector<float>>>& p,
    std::vector<std::vector<std::vector<float>>>& theta
) 
{
    // Get the number of rows, columns, and levels.
    int NR = temperature.size();
    if (NR == 0) return;
    int NTH = temperature[0].size();
    if (NTH == 0) return;
    int NZ = temperature[0][0].size();

    theta.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));

    // Iterate over all radial points
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal points
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical points
            for (int k = 0; k < NZ; ++k) 
            {
                theta[i][j][k] = static_cast<float>(
                    temperature_to_theta(temperature[i][j][k], p[i][j][k])
                );
            }
        }
    }
}
