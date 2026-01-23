#include "simulation.hpp"
#include "radiation/factory.hpp"
#include <iostream>
#include <vector>

/*This file contains the implementation of the radiation system.
It manages the initialization of the radiation system and the computation of the radiation system.*/

// Global radiation scheme instance and configuration
std::unique_ptr<RadiationSchemeBase> radiation_scheme = nullptr;
RadiationConfig global_radiation_config;

// Radiation tendency field
std::vector<std::vector<std::vector<float>>> dtheta_dt_rad;

/*This function initializes the radiation scheme.
Takes in the scheme name and configuration and initializes the radiation scheme.*/
void initialize_radiation(const std::string& scheme_name, const RadiationConfig& cfg) 
{
    try 
    {
        global_radiation_config = cfg;
        radiation_scheme = create_radiation_scheme(scheme_name);
        radiation_scheme->initialize(global_radiation_config);
        std::cout << "Initialized radiation scheme: " << scheme_name << std::endl;

        // Resize radiation tendency field
        dtheta_dt_rad.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));

    } 
    catch (const std::exception& e) 
    {
        std::cerr << "Error initializing radiation: " << e.what() << std::endl;
        throw;
    }
}

/*This function steps the radiation forward in time.
Takes in the current time and steps the radiation forward in time.*/
void step_radiation(double current_time) 
{
    // If the radiation scheme is not set, return.
    if (!radiation_scheme) return;

    // Iterate over the rows, columns, and levels and step the radiation forward in time.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over the columns and step the radiation forward in time.
        for (int j = 0; j < NTH; ++j) 
        {
            // Extract vertical column data
            RadiationColumnStateView col;

            // Convert pressure to double vector
            std::vector<double> p_col(NZ);

            // Iterate over the levels and convert the pressure to a double vector.
            for (int k = 0; k < NZ; ++k) 
            {
                p_col[k] = static_cast<double>(p[i][j][k]);
            }
            col.p = &p_col;

            // Convert theta to temperature
            std::vector<double> T_col(NZ);

            // Iterate over the levels and convert the theta to temperature.
            for (int k = 0; k < NZ; ++k) 
            {
                // theta to T conversion: T = theta * (p/p0)^(R/cp)
                double kappa = R_d / cp;
                T_col[k] = theta[i][j][k] * pow(p[i][j][k] / p0, kappa);
            }
            col.T = &T_col;

            // Convert density and set vertical coordinates
            std::vector<double> rho_col(NZ);
            std::vector<double> dz_col(NZ, dz);

            // Iterate over the levels and convert the density to a double vector.
            for (int k = 0; k < NZ; ++k) 
            {
                rho_col[k] = static_cast<double>(rho[i][j][k]);
            }
            col.rho = &rho_col;
            col.dz = &dz_col;

            // Solar geometry (placeholder - would need proper solar position calculation)
            col.mu0 = 0.5;  // cos(solar zenith angle) - placeholder
            col.S0 = S0_solar;

            // Surface properties (could be made configurable)
            col.Tsfc = T_col[0];  // assume surface = lowest level temperature
            col.albedo_sw = 0.2;  // grassland albedo
            col.emissivity_lw = 0.95; // typical surface emissivity

            // Optional: water vapor (if available)
            // If the water vapor field is not set, return.
            if (!qv.empty()) 
            {
                std::vector<double> qv_col(NZ);

                // Iterate over the levels and convert the water vapor to a double vector.
                for (int k = 0; k < NZ; ++k) 
                {
                    qv_col[k] = static_cast<double>(qv[i][j][k]);
                }
                col.qv = &qv_col;
            }

            // Compute radiation
            RadiationColumnTendencies tend;
            radiation_scheme->compute_column(global_radiation_config, col, tend);

            // Iterate over the levels and convert the temperature tendency to theta tendency.
            for (int k = 0; k < NZ; ++k) 
            {
                double kappa = R_d / cp;
                dtheta_dt_rad[i][j][k] = tend.dTdt_rad[k] * pow(p0 / p[i][j][k], kappa);
            }
        }
    }
}
