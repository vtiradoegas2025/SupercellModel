#pragma once
#include "field3d.hpp"

/*This header file contains the public interface for scalar advection.
The advection component provides stable 3D advection for Field3D scalars
using the numerics framework's TVD scheme with directional splitting.*/

/*This function advects a scalar field in 3D using directional splitting.
Takes in the scalar field, timestep, and diffusion coefficient.
Uses the numerics framework's TVD scheme for stability.
Performs advection in r, theta, and z directions sequentially.*/
void advect_scalar_3d(Field3D& scalar, double dt, double kappa = 0.01);

/*This function advects multiple scalar fields (for thermodynamics).
Takes in timestep and diffusion coefficient.
Advects theta, qv, qc, qr, qi, qs, qh, qg fields.*/
void advect_thermodynamics_3d(double dt_advect, double kappa_theta = 0.01, double kappa_moisture = 0.01);

/*This function advects the tracer field.
Takes in timestep and diffusion coefficient.*/
void advect_tracer_3d(double dt_advect, double kappa = 0.01);
