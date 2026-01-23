#include "simulation.hpp"
#include "boundary_layer/factory.hpp"
#include <iostream>
#include <vector>

/*This file contains the implementation of the boundary layer scheme.
It manages the initialization of the boundary layer scheme and the computation of the boundary layer scheme.*/

// Global boundary layer scheme instance and configuration
std::unique_ptr<BoundaryLayerSchemeBase> boundary_layer_scheme = nullptr;
BoundaryLayerConfig global_boundary_layer_config;
SurfaceConfig global_surface_config;

// Boundary layer tendency fields
std::vector<std::vector<std::vector<float>>> dtheta_dt_pbl;  // potential temperature tendency
std::vector<std::vector<std::vector<float>>> dqv_dt_pbl;     // moisture tendency
std::vector<std::vector<std::vector<float>>> du_dt_pbl;      // u-wind tendency
std::vector<std::vector<std::vector<float>>> dv_dt_pbl;      // v-wind tendency
std::vector<std::vector<std::vector<float>>> dtke_dt_pbl;    // TKE tendency (MYNN only)

/*This function initializes the boundary layer scheme.
Takes in the scheme name, configuration, and surface configuration and initializes the boundary layer scheme.*/
void initialize_boundary_layer(const std::string& scheme_name,
                              const BoundaryLayerConfig& cfg,
                              const SurfaceConfig& sfc) {
    try 
    {
        global_boundary_layer_config = cfg;
        global_surface_config = sfc;
        boundary_layer_scheme = create_boundary_layer_scheme(scheme_name);
        boundary_layer_scheme->initialize(global_boundary_layer_config);

        // Resize boundary layer tendency fields
        dtheta_dt_pbl.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
        dqv_dt_pbl.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
        du_dt_pbl.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
        dv_dt_pbl.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
        dtke_dt_pbl.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));

        std::cout << "Initialized boundary layer scheme: " << scheme_name << std::endl;
        std::cout << "  PBL cadence: " << global_boundary_layer_config.dt_pbl << " s" << std::endl;
        std::cout << "  Surface fluxes: " << (global_boundary_layer_config.apply_surface_fluxes ? "enabled" : "disabled") << std::endl;

    } 
    catch (const std::exception& e) 
    {
        std::cerr << "Error initializing boundary layer: " << e.what() << std::endl;
        throw;
    }
}

/*This function steps the boundary layer forward in time.
Takes in the current time and steps the boundary layer forward in time.*/
void step_boundary_layer(double current_time) 
{
    if (!boundary_layer_scheme) return;

    // Only compute boundary layer at specified cadence
    static double last_boundary_layer_time = -global_boundary_layer_config.dt_pbl;

    // If the current time is less than the last boundary layer time, return.
    if (current_time - last_boundary_layer_time < global_boundary_layer_config.dt_pbl) 
    {
        // Hold previous tendencies
        return;
    }

    last_boundary_layer_time = current_time;

    // For each horizontal column (r, θ)
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over the columns and step the boundary layer forward in time.
        for (int j = 0; j < NTH; ++j) 
        {
            // Extract vertical column data
            BoundaryLayerColumnStateView col;

            // Convert pressure to double vector
            std::vector<double> p_col(NZ);

            // Iterate over the levels and convert the pressure to a double vector.
            for (int k = 0; k < NZ; ++k) 
            {
                p_col[k] = static_cast<double>(p[i][j][k]);
            }
            col.p = &p_col;

            // Convert theta to double vector
            std::vector<double> theta_col(NZ);

            // Iterate over the levels and convert the theta to a double vector.
            for (int k = 0; k < NZ; ++k) 
            {
                theta_col[k] = static_cast<double>(theta[i][j][k]);
            }
            col.theta = &theta_col;

            // Convert moisture
            std::vector<double> qv_col(NZ);

            // Iterate over the levels and convert the moisture to a double vector.
            for (int k = 0; k < NZ; ++k) 
            {
                qv_col[k] = static_cast<double>(qv[i][j][k]);
            }
            col.qv = &qv_col;

            // Convert wind components
            std::vector<double> u_col(NZ);
            std::vector<double> v_col(NZ);

            // Iterate over the levels and convert the wind components to a double vector.
            for (int k = 0; k < NZ; ++k) 
            {
                u_col[k] = static_cast<double>(u[i][j][k]);
                v_col[k] = static_cast<double>(v_theta[i][j][k]);
            }
            col.u = &u_col;
            col.v = &v_col;

            // Convert density
            std::vector<double> rho_col(NZ);

            // Iterate over the levels and convert the density to a double vector.
            for (int k = 0; k < NZ; ++k) 
            {
                rho_col[k] = static_cast<double>(rho[i][j][k]);
            }
            col.rho = &rho_col;

            // Set up vertical coordinates
            std::vector<double> z_int_col(NZ + 1);
            z_int_col[0] = 0.0;  // surface

            // Iterate over the levels and convert the vertical coordinates to a double vector.
            for (int k = 1; k <= NZ; ++k) 
            {
                z_int_col[k] = z_int_col[k-1] + dz;
            }
            col.z_int = &z_int_col;

            // Surface layer inputs (lowest model level)
            col.u_sfc = static_cast<double>(u[i][j][0]);
            col.v_sfc = static_cast<double>(v_theta[i][j][0]);
            col.theta_sfc = static_cast<double>(theta[i][j][0]);
            col.qv_sfc = static_cast<double>(qv[i][j][0]);
            col.z_sfc = dz * 0.5;  // height of lowest level

            // Optional TKE field (for MYNN)

            // If the TKE field is not set, return.
            if (!tke.empty()) 
            {
                std::vector<double> tke_col(NZ);
                for (int k = 0; k < NZ; ++k) {
                    tke_col[k] = static_cast<double>(tke[i][j][k]);
                }
                col.tke = &tke_col;
            }

            // Compute boundary layer
            BoundaryLayerTendencies tend;
            boundary_layer_scheme->compute_column(global_boundary_layer_config, global_surface_config, col, tend);

            // Store tendencies (convert back to float)
            // Iterate over the levels and store the tendencies.
            for (int k = 0; k < NZ; ++k) 
            {
                dtheta_dt_pbl[i][j][k] = static_cast<float>(tend.dthetadt_pbl[k]);
                dqv_dt_pbl[i][j][k] = static_cast<float>(tend.dqvdt_pbl[k]);
                du_dt_pbl[i][j][k] = static_cast<float>(tend.dudt_pbl[k]);
                dv_dt_pbl[i][j][k] = static_cast<float>(tend.dvdt_pbl[k]);
                dtke_dt_pbl[i][j][k] = static_cast<float>(tend.dtkedt_pbl[k]);
            }
        }
    }
}
