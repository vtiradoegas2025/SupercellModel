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
    const Field3D& theta,
    const Field3D& p,
    Field3D& temperature
) 
{
    // Get the number of rows, columns, and levels.
    int NR = theta.size_r();
    if (NR == 0) return;
    int NTH = theta.size_th();
    if (NTH == 0) return;
    int NZ = theta.size_z();

    temperature.resize(NR, NTH, NZ, 0.0f);

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
                    theta_to_temperature(static_cast<double>(theta[i][j][k]), static_cast<double>(p[i][j][k]))
                );
            }
        }
    }
}

/*This function converts the temperature field to the potential temperature field.
Takes in the temperature field, the pressure field, and the potential temperature field
and converts the temperature field to the potential temperature field.*/
void thermodynamics::convert_temperature_to_theta_field(
    const Field3D& temperature,
    const Field3D& p,
    Field3D& theta
) 
{
    // Get the number of rows, columns, and levels.
    int NR = temperature.size_r();
    if (NR == 0) return;
    int NTH = temperature.size_th();
    if (NTH == 0) return;
    int NZ = temperature.size_z();

    theta.resize(NR, NTH, NZ, 0.0f);

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
                    temperature_to_theta(static_cast<double>(temperature[i][j][k]), static_cast<double>(p[i][j][k]))
                );
            }
        }
    }
}
