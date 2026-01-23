#include "simulation.hpp"
#include "dynamics/factory.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>


/*This file contains the implementation of the dynamics scheme.
It manages the initialization of the dynamics scheme and the computation of the dynamics scheme.*/

// Global dynamics scheme instance
std::unique_ptr<DynamicsScheme> dynamics_scheme = nullptr;

// Vorticity diagnostic fields
std::vector<std::vector<std::vector<float>>> vorticity_r;
std::vector<std::vector<std::vector<float>>> vorticity_theta;
std::vector<std::vector<std::vector<float>>> vorticity_z;
std::vector<std::vector<std::vector<float>>> stretching_term;
std::vector<std::vector<std::vector<float>>> tilting_term;
std::vector<std::vector<std::vector<float>>> baroclinic_term;

// Angular momentum diagnostics
std::vector<std::vector<std::vector<float>>> angular_momentum;
std::vector<std::vector<std::vector<float>>> angular_momentum_tendency;

// Pressure diagnostics
std::vector<std::vector<std::vector<float>>> p_prime;
std::vector<std::vector<std::vector<float>>> dynamic_pressure;
std::vector<std::vector<std::vector<float>>> buoyancy_pressure;

/*This function initializes the dynamics scheme.
Takes in the scheme name and initializes the dynamics scheme.*/
void initialize_dynamics(const std::string& scheme_name) 
{
    try 
    {
        dynamics_scheme = create_dynamics_scheme(scheme_name);
        std::cout << "Initialized dynamics scheme: " << scheme_name << std::endl;

        // Resize diagnostic fields
        vorticity_r.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
        vorticity_theta.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
        vorticity_z.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
        stretching_term.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
        tilting_term.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
        baroclinic_term.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
        angular_momentum.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
        angular_momentum_tendency.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
        p_prime.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
        dynamic_pressure.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
        buoyancy_pressure.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));

    } 
    catch (const std::runtime_error& e) 
    {
        std::cerr << "Error initializing dynamics scheme: " << e.what() << std::endl;
        throw;
    }
}

/*This function steps the dynamics forward in time.
Takes in the timestep and the current time and steps the dynamics forward in time.*/
void step_dynamics_new(double dt_dynamics, double current_time) 
{
    // If the dynamics scheme is not set, return.
    if (!dynamics_scheme) 
    {
        std::cerr << "Warning: No dynamics scheme initialized, using old dynamics" << std::endl;
        step_dynamics_old(current_time); // Fall back to old dynamics
        return;
    }

    // Create temporary arrays for tendencies
    std::vector<std::vector<std::vector<float>>> du_r_dt(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    std::vector<std::vector<std::vector<float>>> du_theta_dt(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    std::vector<std::vector<std::vector<float>>> du_z_dt(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    std::vector<std::vector<std::vector<float>>> drho_dt(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    std::vector<std::vector<std::vector<float>>> dp_dt(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));

    // Compute tendencies using the dynamics scheme
    dynamics_scheme->compute_momentum_tendencies(
        u, v_theta, w, rho, p, theta, dt_dynamics,
        du_r_dt, du_theta_dt, du_z_dt, drho_dt, dp_dt
    );

    //Iterate over all grid points and apply tendencies
    for (int i = 0; i < NR; ++i) 
    {
        //Iterate over all azimuthal angles
        for (int j = 0; j < NTH; ++j) 
        {
            //Iterate over all vertical levels
            for (int k = 0; k < NZ; ++k) 
            {
                // Apply tendencies (dynamics + boundary layer)
                u[i][j][k] += (du_r_dt[i][j][k] + du_dt_pbl[i][j][k]) * dt_dynamics;
                v_theta[i][j][k] += (du_theta_dt[i][j][k] + dv_dt_pbl[i][j][k]) * dt_dynamics;
                w[i][j][k] += du_z_dt[i][j][k] * dt_dynamics;
                rho[i][j][k] += drho_dt[i][j][k] * dt_dynamics;
                p[i][j][k] += dp_dt[i][j][k] * dt_dynamics;
            }
        }
    }

    // Apply radiation (updates theta tendencies)
    step_radiation(current_time);

    // Apply microphysics at the same timestep as dynamics
    step_microphysics(dt_dynamics);

    // Calculate radar observables from updated microphysics state
    calculate_radar_reflectivity();

    // Apply subgrid turbulence
    TurbulenceTendencies turb_tend;
    step_turbulence(current_time, turb_tend);
    // Add turbulence tendencies to momentum and scalars

    // Iterate over the rows, columns, and levels and add the turbulence tendencies to the momentum and scalars.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over the columns and add the turbulence tendencies to the momentum and scalars.
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over the levels and add the turbulence tendencies to the momentum and scalars.
            for (int k = 0; k < NZ; ++k) 
            {
                u[i][j][k] += turb_tend.dudt_sgs[i][j][k] * dt_dynamics;
                v_theta[i][j][k] += turb_tend.dvdt_sgs[i][j][k] * dt_dynamics;
                w[i][j][k] += turb_tend.dwdt_sgs[i][j][k] * dt_dynamics;
                theta[i][j][k] += turb_tend.dthetadt_sgs[i][j][k] * dt_dynamics;

                // If the moisture field is not set, return.
                if (!qv.empty()) 
                {
                    qv[i][j][k] += turb_tend.dqvdt_sgs[i][j][k] * dt_dynamics;
                }

                // If the TKE field is not set, return.
                if (!tke.empty()) {
                    tke[i][j][k] += turb_tend.dtkedt_sgs[i][j][k] * dt_dynamics;
                    tke[i][j][k] = std::max(0.001f, tke[i][j][k]);  // ensure positive TKE
                }
            }
        }
    }

    // Compute dynamics diagnostics
    compute_dynamics_diagnostics();

    // Apply boundary conditions
    apply_boundary_conditions();
}

/*This function computes the dynamics diagnostics.
Takes in the dynamics scheme and computes the dynamics diagnostics.*/
void compute_dynamics_diagnostics() 
{
    if (!dynamics_scheme) 
    {
        return;
    }

    // Compute vorticity diagnostics
    dynamics_scheme->compute_vorticity_diagnostics(
        u, v_theta, w, rho, p,
        vorticity_r, vorticity_theta, vorticity_z,
        stretching_term, tilting_term, baroclinic_term
    );

    // Compute angular momentum diagnostics (if supported)
    dynamics_scheme->compute_angular_momentum(
        u, v_theta, angular_momentum, angular_momentum_tendency
    );

    // Compute pressure diagnostics (if supported)
    dynamics_scheme->compute_pressure_diagnostics(
        u, v_theta, w, rho, theta,
        p_prime, dynamic_pressure, buoyancy_pressure
    );
}



// === Apply boundary conditions ===
void apply_boundary_conditions()
{
    // ================== Boundary Conditions ==================
    // Radial boundaries (i=0 and i=NR-1): reflective for u, zero-gradient for others

    //Iterate over all azimuthal angles
    for (int j = 0; j < NTH; ++j) 
    {
        for (int k = 0; k < NZ; ++k) 
        {
            //Iterate over all vertical levels
            u[0][j][k] = -u[1][j][k]; // reflective boundary condition
            u[NR-1][j][k] = -u[NR-2][j][k];
            w[0][j][k] = w[1][j][k]; // no-penetration boundary condition
            w[NR-1][j][k] = w[NR-2][j][k];
            rho[0][j][k] = rho[1][j][k]; // zero-gradient boundary condition
            rho[NR-1][j][k] = rho[NR-2][j][k];
            p[0][j][k] = p[1][j][k]; // zero-gradient boundary condition
            p[NR-1][j][k] = p[NR-2][j][k];
        }
    }

    // Vertical boundaries (k=0 bottom, k=NZ-1 top): no-penetration for w, zero-gradient for others
    for (int i = 0; i < NR; ++i) 
    {
        //Iterate over all azimuthal angles
        for (int j = 0; j < NTH; ++j) 
        {
            //Iterate over all vertical levels
            w[i][j][0] = 0.0f; w[i][j][NZ-1] = 0.0f; // no-penetration boundary condition
            u[i][j][0] = u[i][j][1]; u[i][j][NZ-1] = u[i][j][NZ-2];
            rho[i][j][0] = rho[i][j][1]; rho[i][j][NZ-1] = rho[i][j][NZ-2]; // zero-gradient boundary condition for rho
            p[i][j][0] = p[i][j][1]; p[i][j][NZ-1] = p[i][j][NZ-2]; // zero-gradient boundary condition for p
            theta[i][j][0] = theta[i][j][1]; theta[i][j][NZ-1] = theta[i][j][NZ-2]; // zero-gradient for theta
            qv[i][j][0] = qv[i][j][1]; qv[i][j][NZ-1] = qv[i][j][NZ-2]; // zero-gradient for qv
            qc[i][j][0] = qc[i][j][1]; qc[i][j][NZ-1] = qc[i][j][NZ-2]; // zero-gradient for qc
            qr[i][j][0] = qr[i][j][1]; qr[i][j][NZ-1] = qr[i][j][NZ-2]; // zero-gradient for qr
            qi[i][j][0] = qi[i][j][1]; qi[i][j][NZ-1] = qi[i][j][NZ-2]; // zero-gradient for qi
            qs[i][j][0] = qs[i][j][1]; qs[i][j][NZ-1] = qs[i][j][NZ-2]; // zero-gradient for qs
        }
    }

    // Radial boundaries for thermodynamic fields (zero-gradient)
    for (int j = 0; j < NTH; ++j) 
    {
        //Iterate over all vertical levels
        for (int k = 0; k < NZ; ++k) 
        {
            //Zero-gradient boundary condition for theta
            theta[0][j][k] = theta[1][j][k]; theta[NR-1][j][k] = theta[NR-2][j][k];
            qv[0][j][k] = qv[1][j][k]; qv[NR-1][j][k] = qv[NR-2][j][k];
            qc[0][j][k] = qc[1][j][k]; qc[NR-1][j][k] = qc[NR-2][j][k];
            qr[0][j][k] = qr[1][j][k]; qr[NR-1][j][k] = qr[NR-2][j][k];
            qi[0][j][k] = qi[1][j][k]; qi[NR-1][j][k] = qi[NR-2][j][k];
            qs[0][j][k] = qs[1][j][k]; qs[NR-1][j][k] = qs[NR-2][j][k];
        }
    }
}


// === Compute one timestep using modular dynamics system ===

/*This function steps the dynamics forward in time.
Takes in the current time and steps the dynamics forward in time.*/
void step_dynamics(double current_time)
{
    // If the dynamics scheme is set, step the dynamics forward in time using the new modular dynamics system.
    if (dynamics_scheme)
    {
        step_dynamics_new(dt, current_time);
        return;
    }

    // Fallback: use old dynamics
    step_dynamics_old(current_time);
}


/*This function steps the dynamics forward in time using the old dynamics system.
Takes in the current time and steps the dynamics forward in time using the old dynamics system.*/
void step_dynamics_old(double current_time)

    // If the dynamics scheme is set, step the dynamics forward in time using the new modular dynamics system.
    if (dynamics_scheme)
    {
        step_dynamics_new(dt, current_time);
        return;
    }

    // Fallback: simplified approach for stability when no dynamics scheme is loaded
    // This removes the complex momentum equations that were causing NaN issues

    // CFL-limited timestep
    double max_speed = 1e-6;

    // Iterate over the rows, columns, and levels and compute the maximum speed.
    for (int i = 1; i < NR - 1; ++i) 
    {
        // Iterate over the columns and compute the maximum speed.
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over the levels and compute the maximum speed.
            for (int k = 1; k < NZ - 1; ++k) 
            {
                double sp = std::max({std::abs((double)u[i][j][k]), std::abs((double)w[i][j][k]), std::abs((double)v_theta[i][j][k])});
                if (sp > max_speed) max_speed = sp; 
        }
    }

    double cfl_num = 0.5;
    double dt_eff = dt;

    // If the maximum speed is greater than 1e-6, compute the CFL-limited timestep
    if (max_speed > 1e-6) 
    {
        double dl = std::min(dr, dz);
        double cfl_dt = cfl_num * dl / max_speed;
        dt_eff = std::min(dt, cfl_dt);
    }

    // Apply numerical stability and reset any bad values
    for (int i = 0; i < NR; ++i) 
    {
        //Iterate over all azimuthal angles
        for (int j = 0; j < NTH; ++j) 
        {
            //Iterate over all vertical levels
            for (int k = 0; k < NZ; ++k) 
            {
                //If the vertical level is within the domain, set the density to the base state
                if (k >= 0 && k < NZ) 
                {
                    float base_rho = rho0_base[k];

                    //If the base density is not NaN and is greater than 0 and is finite, set the density to the base state
                    if (!std::isnan(base_rho) && base_rho > 0 && std::isfinite(base_rho)) 
                    {
                        rho[i][j][k] = base_rho;
                    } 

                    //If the base density is NaN or is less than 0 or is not finite, set the density to 1.0
                    else 
                    {
                        rho[i][j][k] = 1.0f;  // Fallback
                    }
                } 

                //If the vertical level is not within the domain, set the density to 1.0
                else 
                {
                    rho[i][j][k] = 1.0f;  // Bounds check
                }

                // Reset any NaN/inf values to reasonable defaults
                if (std::isnan(p[i][j][k]) || std::isinf(p[i][j][k])) 
                {
                    p[i][j][k] = p0;
                }
                //If the temperature is not NaN and is finite, set the temperature to the base state
                if (std::isnan(theta[i][j][k]) || std::isinf(theta[i][j][k])) 
                {
                    theta[i][j][k] = theta0;
                }

                //If the u velocity is not NaN and is finite, set the u velocity to 0
                if (std::isnan(u[i][j][k]) || std::isinf(u[i][j][k])) 
                {
                    u[i][j][k] = 0.0f;
                }

                //If the w velocity is not NaN and is finite, set the w velocity to 0
                if (std::isnan(w[i][j][k]) || std::isinf(w[i][j][k])) 
                {
                    w[i][j][k] = 0.0f;
                }

                //If the v_theta velocity is not NaN and is finite, set the v_theta velocity to 0
                if (std::isnan(v_theta[i][j][k]) || std::isinf(v_theta[i][j][k])) 
                {
                    v_theta[i][j][k] = 0.0f;
                }

                // Apply bounds
                p[i][j][k] = std::max(1000.0f, std::min(200000.0f, p[i][j][k]));
                theta[i][j][k] = std::max(200.0f, std::min(500.0f, theta[i][j][k]));
                u[i][j][k] = std::max(-100.0f, std::min(100.0f, u[i][j][k]));
                w[i][j][k] = std::max(-50.0f, std::min(50.0f, w[i][j][k]));
                v_theta[i][j][k] = std::max(-100.0f, std::min(100.0f, v_theta[i][j][k]));
            }
        }
    }

    // Re-enable advection and microphysics
    advect_tracer(dt_eff);
    advect_thermodynamics(dt_eff);
    step_microphysics(dt_eff);

    // Apply boundary conditions (from old implementation)
    // Radial boundaries: reflective for u, zero-gradient for others

    // Iterate over all azimuthal angles to apply boundary conditions
    for (int j = 0; j < NTH; ++j) 
    {
        //Iterate over all vertical levels
        for (int k = 0; k < NZ; ++k) 
        {
            // Only apply boundary conditions if interior values are valid
            if (!std::isnan(u[1][j][k]) && std::isfinite(u[1][j][k])) 
            {
                u[0][j][k] = -u[1][j][k];
            }

            //If the u velocity is not NaN and is finite, set the u velocity to the negative of the u velocity at the next radial level
            if (!std::isnan(u[NR-2][j][k]) && std::isfinite(u[NR-2][j][k])) 
            {
                u[NR-1][j][k] = -u[NR-2][j][k];
            }

            //If the w velocity is not NaN and is finite, set the w velocity to the w velocity at the next radial level
            if (!std::isnan(w[1][j][k]) && std::isfinite(w[1][j][k])) 
            {
                w[0][j][k] = w[1][j][k];
            }

            //If the w velocity is not NaN and is finite, set the w velocity to the negative of the w velocity at the next radial level
            if (!std::isnan(w[NR-2][j][k]) && std::isfinite(w[NR-2][j][k])) 
            {
                w[NR-1][j][k] = w[NR-2][j][k];
            }

            // For scalars, use base values if interior is invalid
            rho[0][j][k] = (!std::isnan(rho[1][j][k]) && std::isfinite(rho[1][j][k])) ? rho[1][j][k] : rho0_base[k];
            rho[NR-1][j][k] = (!std::isnan(rho[NR-2][j][k]) && std::isfinite(rho[NR-2][j][k])) ? rho[NR-2][j][k] : rho0_base[k];

            p[0][j][k] = (!std::isnan(p[1][j][k]) && std::isfinite(p[1][j][k])) ? p[1][j][k] : p0;
            p[NR-1][j][k] = (!std::isnan(p[NR-2][j][k]) && std::isfinite(p[NR-2][j][k])) ? p[NR-2][j][k] : p0;
        }
    }

    // Vertical boundaries: no-penetration for w, zero-gradient for others
    // Iterate over all azimuthal angles to apply boundary conditions
    for (int i = 0; i < NR; ++i) 
    {
        //Iterate over all azimuthal angles
        for (int j = 0; j < NTH; ++j) 
        {
            w[i][j][0] = 0.0f;  // Always zero at bottom
            w[i][j][NZ-1] = 0.0f;  // Always zero at top

            // Use defensive boundary conditions
            u[i][j][0] = (!std::isnan(u[i][j][1]) && std::isfinite(u[i][j][1])) ? u[i][j][1] : 0.0f;
            u[i][j][NZ-1] = (!std::isnan(u[i][j][NZ-2]) && std::isfinite(u[i][j][NZ-2])) ? u[i][j][NZ-2] : 0.0f;

            rho[i][j][0] = (!std::isnan(rho[i][j][1]) && std::isfinite(rho[i][j][1])) ? rho[i][j][1] : rho0_base[0];
            rho[i][j][NZ-1] = (!std::isnan(rho[i][j][NZ-2]) && std::isfinite(rho[i][j][NZ-2])) ? rho[i][j][NZ-2] : rho0_base[NZ-1];

            p[i][j][0] = (!std::isnan(p[i][j][1]) && std::isfinite(p[i][j][1])) ? p[i][j][1] : p0;
            p[i][j][NZ-1] = (!std::isnan(p[i][j][NZ-2]) && std::isfinite(p[i][j][NZ-2])) ? p[i][j][NZ-2] : p0;

            theta[i][j][0] = (!std::isnan(theta[i][j][1]) && std::isfinite(theta[i][j][1])) ? theta[i][j][1] : theta0;
            theta[i][j][NZ-1] = (!std::isnan(theta[i][j][NZ-2]) && std::isfinite(theta[i][j][NZ-2])) ? theta[i][j][NZ-2] : theta0;

            qv[i][j][0] = (!std::isnan(qv[i][j][1]) && std::isfinite(qv[i][j][1])) ? qv[i][j][1] : 0.0f;
            qv[i][j][NZ-1] = (!std::isnan(qv[i][j][NZ-2]) && std::isfinite(qv[i][j][NZ-2])) ? qv[i][j][NZ-2] : 0.0f;

            qc[i][j][0] = (!std::isnan(qc[i][j][1]) && std::isfinite(qc[i][j][1])) ? qc[i][j][1] : 0.0f;
            qc[i][j][NZ-1] = (!std::isnan(qc[i][j][NZ-2]) && std::isfinite(qc[i][j][NZ-2])) ? qc[i][j][NZ-2] : 0.0f;

            qr[i][j][0] = (!std::isnan(qr[i][j][1]) && std::isfinite(qr[i][j][1])) ? qr[i][j][1] : 0.0f;
            qr[i][j][NZ-1] = (!std::isnan(qr[i][j][NZ-2]) && std::isfinite(qr[i][j][NZ-2])) ? qr[i][j][NZ-2] : 0.0f;
        }
    }

    // Radial boundaries for thermodynamic fields

    // Iterate over all azimuthal angles to apply boundary conditions
    for (int j = 0; j < NTH; ++j) 
    {
        for (int k = 0; k < NZ; ++k) 
        {
            theta[0][j][k] = (!std::isnan(theta[1][j][k]) && std::isfinite(theta[1][j][k])) ? theta[1][j][k] : theta0;
            theta[NR-1][j][k] = (!std::isnan(theta[NR-2][j][k]) && std::isfinite(theta[NR-2][j][k])) ? theta[NR-2][j][k] : theta0;

            qv[0][j][k] = (!std::isnan(qv[1][j][k]) && std::isfinite(qv[1][j][k])) ? qv[1][j][k] : 0.0f;
            qv[NR-1][j][k] = (!std::isnan(qv[NR-2][j][k]) && std::isfinite(qv[NR-2][j][k])) ? qv[NR-2][j][k] : 0.0f;

            qc[0][j][k] = (!std::isnan(qc[1][j][k]) && std::isfinite(qc[1][j][k])) ? qc[1][j][k] : 0.0f;
            qc[NR-1][j][k] = (!std::isnan(qc[NR-2][j][k]) && std::isfinite(qc[NR-2][j][k])) ? qc[NR-2][j][k] : 0.0f;

            qr[0][j][k] = (!std::isnan(qr[1][j][k]) && std::isfinite(qr[1][j][k])) ? qr[1][j][k] : 0.0f;
            qr[NR-1][j][k] = (!std::isnan(qr[NR-2][j][k]) && std::isfinite(qr[NR-2][j][k])) ? qr[NR-2][j][k] : 0.0f;
        }
    }
}
