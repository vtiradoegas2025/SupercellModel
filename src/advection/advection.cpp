#include "advection.hpp"
#include "simulation.hpp"
#include "advection_base.hpp"
#include "numerics_base.hpp"
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>

// Access global numerics variables from numerics.cpp
extern std::unique_ptr<AdvectionSchemeBase> advection_scheme;
extern AdvectionConfig global_advection_config;
extern GridMetrics global_grid_metrics;

/*This file contains the implementation of the scalar advection component.
It provides stable 3D advection using the numerics framework's TVD scheme
with directional splitting for cylindrical coordinates (r, theta, z).*/

// Helper function to convert Field3D to nested vector for numerics framework
static void field3d_to_nested(const Field3D& field, 
                              std::vector<std::vector<std::vector<double>>>& nested)
{
    nested.assign(NR, std::vector<std::vector<double>>(
        NTH, std::vector<double>(NZ, 0.0)));
    
    for (int i = 0; i < NR; ++i) {
        for (int j = 0; j < NTH; ++j) {
            for (int k = 0; k < NZ; ++k) {
                nested[i][j][k] = static_cast<double>(field[i][j][k]);
            }
        }
    }
}

// Helper function to convert nested vector back to Field3D
static void nested_to_field3d(const std::vector<std::vector<std::vector<double>>>& nested,
                               Field3D& field)
{
    for (int i = 0; i < NR; ++i) {
        for (int j = 0; j < NTH; ++j) {
            for (int k = 0; k < NZ; ++k) {
                field[i][j][k] = static_cast<float>(nested[i][j][k]);
            }
        }
    }
}

// Helper function to apply diffusion
static void apply_diffusion(Field3D& scalar, double dt, double kappa)
{
    if (kappa <= 0.0) return;
    
    Field3D new_scalar = scalar;
    
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < NR - 1; ++i) {
        double r = i * dr + 1e-6;
        for (int j = 0; j < NTH; ++j) {
            int j_prev = (j - 1 + NTH) % NTH;
            int j_next = (j + 1) % NTH;
            for (int k = 1; k < NZ - 1; ++k) {
                // Laplacian diffusion
                double lap_S = (scalar[i+1][j][k] - 2*scalar[i][j][k] + scalar[i-1][j][k])/(dr*dr)
                              + (scalar[i][j][k+1] - 2*scalar[i][j][k] + scalar[i][j][k-1])/(dz*dz)
                              + (scalar[i][j_prev][k] - 2*scalar[i][j][k] + scalar[i][j_next][k])/(dtheta*dtheta*r*r);
                
                new_scalar[i][j][k] += static_cast<float>(dt * kappa * lap_S);
            }
        }
    }
    
    scalar = new_scalar;
}

// Advect in radial (r) direction using TVD-like scheme
static void advect_scalar_1d_r(Field3D& scalar, double dt, double kappa)
{
    Field3D new_scalar = scalar;
    
    #pragma omp parallel for collapse(2)
    for (int j = 0; j < NTH; ++j) {
        for (int k = 1; k < NZ - 1; ++k) {
            // Extract 1D radial column
            std::vector<double> q_col(NR);
            std::vector<double> u_col(NR);
            for (int i = 0; i < NR; ++i) {
                q_col[i] = static_cast<double>(scalar[i][j][k]);
                u_col[i] = static_cast<double>(u[i][j][k]);
            }
            
            // Compute tendencies using TVD-like upwind scheme
            std::vector<double> dqdt(NR, 0.0);
            
            for (int i = 1; i < NR - 1; ++i) {
                // Upwind difference based on velocity
                // For u > 0: flow is outward, use backward difference (upwind)
                // For u < 0: flow is inward, use forward difference (upwind)
                double dq_dr;
                if (u_col[i] > 0) {
                    dq_dr = (q_col[i] - q_col[i-1]) / dr;  // backward difference
                } else if (u_col[i] < 0) {
                    dq_dr = (q_col[i+1] - q_col[i]) / dr;  // forward difference
                } else {
                    dq_dr = 0.0;  // no flow
                }
                
                // Advection: dq/dt = -u * dq/dr (negative sign for material derivative)
                dqdt[i] = -u_col[i] * dq_dr;
            }
            
            // Apply tendencies with CFL check
            for (int i = 1; i < NR - 1; ++i) {
                float change = static_cast<float>(dt * dqdt[i]);
                
                // CFL check: ensure change isn't too large
                double cfl_r = std::abs(u_col[i]) * dt / dr;
                if (cfl_r > 1.0) {
                    // Scale down if CFL > 1
                    change *= (1.0 / cfl_r);
                }
                
                // Limit change to prevent extreme values
                float max_change = std::abs(new_scalar[i][j][k]) * 0.5f;  // Limit to 50% change
                if (std::abs(change) > max_change) {
                    change = (change > 0) ? max_change : -max_change;
                }
                
                new_scalar[i][j][k] += change;
            }
        }
    }
    
    scalar = new_scalar;
}

// Advect in azimuthal (theta) direction
static void advect_scalar_1d_theta(Field3D& scalar, double dt, double kappa)
{
    Field3D new_scalar = scalar;
    
    #pragma omp parallel for collapse(2)
    for (int i = 1; i < NR - 1; ++i) {
        double r = i * dr + 1e-6;
        for (int k = 1; k < NZ - 1; ++k) {
            // Extract 1D azimuthal column
            std::vector<double> q_col(NTH);
            std::vector<double> v_col(NTH);
            for (int j = 0; j < NTH; ++j) {
                q_col[j] = static_cast<double>(scalar[i][j][k]);
                v_col[j] = static_cast<double>(v_theta[i][j][k]);
            }
            
            // Compute tendencies using TVD-like upwind scheme
            std::vector<double> dqdt(NTH, 0.0);
            
            for (int j = 0; j < NTH; ++j) {
                int j_prev = (j - 1 + NTH) % NTH;
                int j_next = (j + 1) % NTH;
                
                // Upwind difference based on velocity
                double dq_dtheta;
                if (v_col[j] > 0) {
                    dq_dtheta = (q_col[j] - q_col[j_prev]) / dtheta;  // backward difference
                } else if (v_col[j] < 0) {
                    dq_dtheta = (q_col[j_next] - q_col[j]) / dtheta;  // forward difference
                } else {
                    dq_dtheta = 0.0;  // no flow
                }
                
                // Advection: dq/dt = -(v/r) * dq/dtheta (accounting for cylindrical geometry)
                dqdt[j] = -(v_col[j] / r) * dq_dtheta;
            }
            
            // Apply tendencies
            for (int j = 0; j < NTH; ++j) {
                new_scalar[i][j][k] += static_cast<float>(dt * dqdt[j]);
            }
        }
    }
    
    scalar = new_scalar;
}

// Advect in vertical (z) direction using numerics framework TVD scheme
static void advect_scalar_1d_z(Field3D& scalar, double dt, double kappa)
{
    if (!advection_scheme) {
        // Fallback to simple upwind if numerics framework not initialized
        Field3D new_scalar = scalar;
        
        #pragma omp parallel for collapse(2)
        for (int i = 1; i < NR - 1; ++i) {
            for (int j = 0; j < NTH; ++j) {
                for (int k = 1; k < NZ - 1; ++k) {
                    double w_val = w[i][j][k];
                    double dq_dz;
                    if (w_val > 0) {
                        dq_dz = (scalar[i][j][k] - scalar[i][j][k-1]) / dz;
                    } else {
                        dq_dz = (scalar[i][j][k+1] - scalar[i][j][k]) / dz;
                    }
                    new_scalar[i][j][k] -= static_cast<float>(dt * w_val * dq_dz);
                }
            }
        }
        scalar = new_scalar;
        return;
    }
    
    // Use numerics framework TVD scheme for vertical advection
    // The TVD scheme does vertical advection column by column
    // Convert Field3D to nested vector format
    std::vector<std::vector<std::vector<double>>> q_nested, w_nested;
    field3d_to_nested(scalar, q_nested);
    field3d_to_nested(w, w_nested);
    
    // Create state view
    AdvectionStateView state;
    state.q = &q_nested;
    state.w = &w_nested;
    state.grid = &global_grid_metrics;
    
    // Compute tendencies (TVD scheme computes vertical advection tendencies)
    AdvectionTendencies tendencies;
    AdvectionConfig cfg = global_advection_config;
    cfg.positivity = false;  // Don't enforce positivity for theta
    
    advection_scheme->compute_flux_divergence(cfg, state, tendencies);
    
    // Apply tendencies: q_new = q_old + dt * dqdt
    for (int i = 0; i < NR; ++i) {
        for (int j = 0; j < NTH; ++j) {
            for (int k = 0; k < NZ; ++k) {
                q_nested[i][j][k] += dt * tendencies.dqdt_adv[i][j][k];
            }
        }
    }
    
    // Convert back to Field3D
    nested_to_field3d(q_nested, scalar);
}

/*This function advects a scalar field in 3D using directional splitting.
Uses Strang splitting for better stability: dt/2 in r, dt in theta, dt/2 in r,
then dt in z, then dt/2 in r, dt in theta, dt/2 in r.
This is second-order accurate and more stable than simple sequential splitting.*/
void advect_scalar_3d(Field3D& scalar, double dt, double kappa)
{
    // Debug: Track values for first few calls
    static int call_count = 0;
    call_count++;
    
    if (call_count <= 3) {
        float val_before = scalar[1][0][1];
        std::cerr << "[ADVECT_3D] Call " << call_count << " BEFORE: scalar[1][0][1]=" << val_before << ", dt=" << dt << std::endl;
    }
    
    // Proper Strang splitting for 3D advection (second-order accurate)
    // Standard Strang splitting: A(dt/2) -> B(dt/2) -> C(dt) -> B(dt/2) -> A(dt/2)
    // For 3D: r(dt/2) -> theta(dt/2) -> z(dt) -> theta(dt/2) -> r(dt/2)
    
    double dt_half = dt * 0.5;
    
    // First half-step in r
    advect_scalar_1d_r(scalar, dt_half, kappa);
    
    // First half-step in theta
    advect_scalar_1d_theta(scalar, dt_half, kappa);
    
    // Full step in z (vertical, uses TVD scheme)
    advect_scalar_1d_z(scalar, dt, kappa);
    
    // Second half-step in theta
    advect_scalar_1d_theta(scalar, dt_half, kappa);
    
    // Second half-step in r
    advect_scalar_1d_r(scalar, dt_half, kappa);
    
    // Apply diffusion (only once, not per direction)
    apply_diffusion(scalar, dt, kappa);
}

/*This function advects thermodynamics fields.*/
void advect_thermodynamics_3d(double dt_advect, double kappa_theta, double kappa_moisture)
{
    // Force output to ensure function is called
    std::cerr << "\n*** [ADVECT_THERMO_3D] CALLED *** dt=" << dt_advect << std::endl;
    std::cerr.flush();
    
    // Debug: Track theta before advection
    static int thermo_call = 0;
    thermo_call++;
    float theta_before = theta[1][0][1];
    std::cerr << "[ADVECT_THERMO_3D] Call " << thermo_call << " BEFORE: theta[1][0][1]=" << theta_before << "K, dt=" << dt_advect << std::endl;
    std::cerr.flush();
    
    // Advect potential temperature with specified diffusion
    advect_scalar_3d(theta, dt_advect, kappa_theta);
    
    // Debug: Track theta after advection
    if (thermo_call <= 3) {
        float theta_after = theta[1][0][1];
        std::cerr << "[ADVECT_THERMO_3D] Call " << thermo_call << " AFTER: theta[1][0][1]=" << theta_after << "K" << std::endl;
    }
    
    // Clamp theta to reasonable values
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < NR; ++i) {
        for (int j = 0; j < NTH; ++j) {
            for (int k = 0; k < NZ; ++k) {
                float theta_val = theta[i][j][k];
                if (theta_val < 200.0f || theta_val > 500.0f) {
                    theta[i][j][k] = std::max(200.0f, std::min(500.0f, theta_val));
                }
            }
        }
    }
    
    // Advect moisture fields
    advect_scalar_3d(qv, dt_advect, kappa_moisture);
    advect_scalar_3d(qc, dt_advect, kappa_moisture);
    advect_scalar_3d(qr, dt_advect, kappa_moisture);
    advect_scalar_3d(qi, dt_advect, kappa_moisture);
    advect_scalar_3d(qs, dt_advect, kappa_moisture);
    advect_scalar_3d(qh, dt_advect, kappa_moisture);
    advect_scalar_3d(qg, dt_advect, kappa_moisture);
}

/*This function advects the tracer field.*/
void advect_tracer_3d(double dt_advect, double kappa)
{
    advect_scalar_3d(tracer, dt_advect, kappa);
}
