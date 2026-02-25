#pragma once

#include "field3d.hpp"

/**
 * @file advection.hpp
 * @brief Public advection entry points used by the simulation runtime.
 *
 * Declares high-level helpers that advect scalar, tracer, and
 * thermodynamic fields on the model grid.
 * These functions wrap the configured advection implementation.
 */

/**
 * @brief Advects a scalar field over one advection step.
 * @param scalar Scalar field updated in place.
 * @param dt Advection step length in seconds.
 * @param kappa Diffusion-like stabilization coefficient.
 */
void advect_scalar_3d(Field3D& scalar, double dt, double kappa = 0.01);

/**
 * @brief Advects thermodynamic fields over one advection step.
 * @param dt_advect Advection step length in seconds.
 * @param kappa_theta Stabilization coefficient for potential temperature.
 * @param kappa_moisture Stabilization coefficient for moisture scalars.
 */
void advect_thermodynamics_3d(double dt_advect,
                              double kappa_theta = 0.01,
                              double kappa_moisture = 0.01);

/**
 * @brief Advects the passive tracer field over one advection step.
 * @param dt_advect Advection step length in seconds.
 * @param kappa Diffusion-like stabilization coefficient.
 */
void advect_tracer_3d(double dt_advect, double kappa = 0.01);

/**
 * @brief Resets advection performance counters.
 */
void reset_advection_perf_stats();

/**
 * @brief Logs a summary of advection performance counters.
 */
void log_advection_perf_summary();
