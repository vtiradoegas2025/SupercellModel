#include <iostream>
#include <cmath>
#include <memory>
#include "simulation.hpp"
#include "microphysics/factory.hpp"
#include "radar/factory.hpp"


/*This file contains the implementation of the equations.
It manages the initialization of the equations and the computation of the equations.*/

// Forward declaration for wind profile computation
void compute_wind_profile(const WindProfile& profile, double z, double& u, double& v);

// Local base-state parameters (not exported)
// p0 is now defined in simulation.hpp

// Grid parameters (initialized with defaults)
int NR = 200;         // radial points
int NTH = 128;        // azimuthal points (theta)
int NZ = 150;         // vertical points

// Grid resolution parameters
double dr = 100.0;    // radial resolution (m)
double dz = 100.0;    // vertical resolution (m)
double dt = 0.1;      // time step (s)
double dtheta = 2.0 * 3.14159265358979323846 / NTH; // azimuthal resolution (rad)

// Simulation time tracking
double simulation_time = 0.0;  // current simulation time [s]

/*This function updates the grid resolution parameters.
Takes in the grid resolution parameters and updates the grid resolution parameters.*/
void update_grid_resolution() 
{
    const double pi = 3.14159265358979323846;
    // dtheta depends only on NTH (azimuthal points)
    dtheta = 2.0 * pi / NTH;
    // dr, dz, dt are now loaded from config in tornado_sim.cpp
}

// Field variables [r][th][z] - will be resized at runtime
std::vector<std::vector<std::vector<float>>> rho;
std::vector<std::vector<std::vector<float>>> p;
std::vector<std::vector<std::vector<float>>> u;
std::vector<std::vector<std::vector<float>>> w;
std::vector<std::vector<std::vector<float>>> v_theta;
std::vector<std::vector<std::vector<float>>> tracer;

// Thermodynamic fields
std::vector<std::vector<std::vector<float>>> theta;
std::vector<std::vector<std::vector<float>>> qv;
std::vector<std::vector<std::vector<float>>> qc;
std::vector<std::vector<std::vector<float>>> qr;
std::vector<std::vector<std::vector<float>>> qi;  // cloud ice mixing ratio
std::vector<std::vector<std::vector<float>>> qs;  // snow mixing ratio
std::vector<std::vector<std::vector<float>>> qh;
std::vector<std::vector<std::vector<float>>> qg;
std::vector<std::vector<std::vector<float>>> tke; // turbulent kinetic energy (MYNN)

// Radar reflectivity field
std::vector<std::vector<std::vector<float>>> radar_reflectivity;

// Anelastic base state
std::vector<double> rho0_base;  // base-state density profile

// Global physics schemes
std::unique_ptr<MicrophysicsScheme> microphysics_scheme;
std::unique_ptr<RadarSchemeBase> radar_scheme;

// Nested grid configuration and fields
NestedGridConfig nested_config;
std::vector<std::vector<std::vector<float>>> nest_rho;
std::vector<std::vector<std::vector<float>>> nest_p;
std::vector<std::vector<std::vector<float>>> nest_u;
std::vector<std::vector<std::vector<float>>> nest_w;
std::vector<std::vector<std::vector<float>>> nest_v_theta;
std::vector<std::vector<std::vector<float>>> nest_theta;
std::vector<std::vector<std::vector<float>>> nest_qv;
std::vector<std::vector<std::vector<float>>> nest_qc;
std::vector<std::vector<std::vector<float>>> nest_qr;
std::vector<std::vector<std::vector<float>>> nest_qh;
std::vector<std::vector<std::vector<float>>> nest_qg;
std::vector<std::vector<std::vector<float>>> nest_tracer;

/*This function resizes the fields.
Takes in the fields and resizes the fields.*/
void resize_fields() 
{
    update_grid_resolution(); // Update resolution parameters first

    rho.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    p.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    u.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    w.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    v_theta.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    tracer.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    theta.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    qv.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    qc.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    qr.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    qi.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    qs.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    qh.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    qg.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    tke.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    radar_reflectivity.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));

    // Boundary layer tendency fields
    dtheta_dt_pbl.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dqv_dt_pbl.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    du_dt_pbl.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dv_dt_pbl.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dtke_dt_pbl.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
}

    /*This function initializes the base state.
    Takes in the base state and initializes the base state.*/
void initialize()
{
    // Pre-compute CAPE scaling factor
    double cape_scaling = global_cape_target / 2500.0; // normalize to classic case

    // Initialize anelastic base state density profile
    rho0_base.resize(NZ);  // Resize vector to match NZ

    // Iterate over the levels and initialize the base state.
    for (int k = 0; k < NZ; ++k)
    {
        double z = k * dz;
        double T = theta0 - 0.0065 * z;  // base temperature profile
        double p_local = p0 * pow(1 - (0.0065 * z / theta0), 5.255);
        rho0_base[k] = std::max(p_local / (R_d * T), 0.1);
    }
    std::cout << "Base state density initialized: rho0_base[0]=" << rho0_base[0]
              << ", rho0_base[" << NZ-1 << "]=" << rho0_base[NZ-1] << std::endl;

    // Iterate over the rows, columns, and levels and initialize the base state.
    for (int i = 0; i < NR; ++i)
    {
        // double r = i * dr;  // radial distance from center (unused)

        // Iterate over the columns and initialize the base state.
        for (int j = 0; j < NTH; ++j)
        {
            // double th = j * dtheta;  // azimuthal angle (unused)

            // Iterate over the levels and initialize the base state.
            for (int k = 0; k < NZ; ++k)
            {
                double z = k * dz;
                // Use hydrostatic balance for pressure and density
                // dP/dz = -ρg, with P = ρRT for ideal gas
                double T = theta0 - 0.0065 * z;  // actual temperature (dry adiabatic lapse rate)
                double p_local = p0 * pow(1 - (0.0065 * z / theta0), 5.255);  // hydrostatic pressure
                double rho_local = p_local / (R_d * T);  // ideal gas law

                p[i][j][k] = static_cast<float>(p_local);
                rho[i][j][k] = static_cast<float>(std::max(rho_local, 0.1));  // ensure positive density

                // Initialize potential temperature profile
                // First calculate actual temperature profile with desired stability
                double T_actual;
                
                // If the level is in the near-surface stable layer, initialize the base state.
                if (z < 1000.0) 
                {
                    // Near-surface stable layer (slight inversion)
                    T_actual = theta0 + 1.0;  // +1K inversion
                } 

                // If the level is in the unstable layer, initialize the base state.
                else if (z < 6000.0) 
                {
                    // Unstable layer with lapse rate scaled by CAPE
                    double lapse_rate = 0.004 + 0.002 * cape_scaling; // steeper lapse for higher CAPE
                    T_actual = theta0 + 1.0 - lapse_rate * (z - 1000.0);
                } 
                else 
                {
                    // Stable upper troposphere
                    T_actual = theta0 + 1.0 - 0.006 * 5000.0 - 0.003 * (z - 6000.0);
                }

                // Convert to potential temperature: θ = T * (p0/p)^(R/cp)
                double kappa = R_d / cp;  // R/cp
                double theta_potential = T_actual * pow(p0 / p_local, kappa);
                theta[i][j][k] = static_cast<float>(theta_potential);

                // Moisture profile: higher near surface, decreasing with height
                // Scale base moisture with CAPE target (higher CAPE = more moisture)
                double base_moisture = 0.016 * cape_scaling;
                double qv_base;

                // If the level is in the boundary layer, initialize the base state.
                if (z < 2000.0) 
                {
                    // Well-mixed boundary layer
                    qv_base = base_moisture;
                } 
                else 
                {
                    // Exponential decrease above boundary layer
                    qv_base = base_moisture * exp(-(z - 2000.0) / 3000.0);
                }
                qv[i][j][k] = static_cast<float>(qv_base);

                qc[i][j][k] = 0.0f;    // no cloud water initially
                qr[i][j][k] = 0.0f;    // no rain water initially
                qi[i][j][k] = 0.0f;    // no cloud ice initially
                qs[i][j][k] = 0.0f;    // no snow initially
                qh[i][j][k] = 0.0f;    // no hail initially
                qg[i][j][k] = 0.0f;    // no graupel initially
                tke[i][j][k] = 0.1f;   // background TKE [m²/s²]

                // Initialize winds from hodograph profile (Cartesian coordinates)
                double wind_u_cart, wind_v_cart;
                compute_wind_profile(global_wind_profile, z, wind_u_cart, wind_v_cart);

                // Convert Cartesian winds to cylindrical coordinates
                // In cylindrical coordinates: u_r = u_cart * cos(theta) + v_cart * sin(theta)
                //                          v_theta = -u_cart * sin(theta) + v_cart * cos(theta)
                double th = j * dtheta;  // azimuthal angle
                double u_r = wind_u_cart * cos(th) + wind_v_cart * sin(th);
                double v_th = -wind_u_cart * sin(th) + wind_v_cart * cos(th);

                u[i][j][k] = static_cast<float>(u_r);
                v_theta[i][j][k] = static_cast<float>(v_th);
                w[i][j][k] = 0.0f;
                tracer[i][j][k] = 0.0f;
            }
        }
    }

    // Add thermal bubble trigger for convection (Weisman-Klemp style)
    // Bubble centered at r=50km, z=1.5km, radius=10km, dtheta=2K
    // int bubble_center_r = static_cast<int>(50000.0 / dr);  // 50km in grid points (unused)
    // int bubble_center_z = static_cast<int>(1500.0 / dz);   // 1.5km in grid points (unused)
    double bubble_radius = 10000.0;  // 10km radius
    double bubble_dtheta = 2.0;      // 2K temperature perturbation

    // Iterate over the rows, columns, and levels and add the thermal bubble trigger for convection.
    for (int i = 0; i < NR; ++i)
    {
        double r_dist = i * dr;

        // Iterate over the columns and add the thermal bubble trigger for convection.
        for (int j = 0; j < NTH; ++j)
        {
            // Iterate over the levels and add the thermal bubble trigger for convection.
            for (int k = 0; k < NZ; ++k)
            {
                double z_dist = k * dz;
                double dist_from_center = sqrt(pow(r_dist - 50000.0, 2) + pow(z_dist - 1500.0, 2));

                // If the distance from the center is less than the bubble radius, add the thermal bubble trigger for convection.
                if (dist_from_center <= bubble_radius)
                {
                    // Gaussian bubble profile
                    double bubble_factor = exp(-pow(dist_from_center / (bubble_radius / 3.0), 2));
                    theta[i][j][k] += bubble_dtheta * bubble_factor;
                }
            }
        }
    }

    // Add initial vortex tracer
    int ic = NR / 4;
    int kc = NZ / 4;
    for (int j = 0; j < NTH; ++j)
    {
        tracer[ic][j][kc] = 1.0f;
    }

    // Initialize nested grid if enabled
    initialize_nested_grid();
}

/*This function initializes the microphysics scheme.
Takes in the scheme name and initializes the microphysics scheme.*/
void initialize_microphysics(const std::string& scheme_name) 
{
    try 
    {
        microphysics_scheme = create_microphysics_scheme(scheme_name);
        std::cout << "Initialized microphysics scheme: " << scheme_name << std::endl;
    } 
    catch (const std::runtime_error& e)
    {
        std::cerr << "Error initializing microphysics: " << e.what() << std::endl;
        // Fallback to Kessler
        microphysics_scheme = create_microphysics_scheme("kessler");
        std::cout << "Falling back to Kessler microphysics scheme" << std::endl;
    }
}


/*This function initializes the radar scheme.
Takes in the scheme name and initializes the radar scheme.*/    
void initialize_radar(const std::string& scheme_name) 
{
    try 
    {
        radar_scheme = RadarFactory::create(scheme_name);

        // Configure radar scheme
        RadarConfig config;
        config.scheme_id = scheme_name;
        config.operator_tier = "fast_da";  // Default tier
        config.has_qr = true;  // Enable rain reflectivity

        radar_scheme->initialize(config, NR, NTH, NZ);
        std::cout << "Initialized radar scheme: " << scheme_name << std::endl;
    } 
    catch (const std::runtime_error& e) {
        std::cerr << "Error initializing radar: " << e.what() << std::endl;
        // For now, radar is optional - don't fail if it can't initialize
        std::cout << "Radar scheme initialization failed, radar calculations disabled" << std::endl;
    }
}

/*This function advects the scalar field.
Takes in the scalar field, the timestep, and the diffusion coefficient and advects the scalar field.*/
void advect_scalar(std::vector<std::vector<std::vector<float>>>& scalar, double dt_advect, double kappa = 0.01)
{
    std::vector<std::vector<std::vector<float>>> new_scalar = scalar;

    // Iterate over the rows, columns, and levels and advect the scalar field.
    for (int i = 1; i < NR - 1; ++i)
    {
        double r = i * dr + 1e-6;

        // Iterate over the columns and advect the scalar field.
        for (int j = 0; j < NTH; ++j)
        {
            int j_prev = (j - 1 + NTH) % NTH;
            for (int k = 1; k < NZ - 1; ++k)
            {
                double dS_dr = (scalar[i][j][k] - scalar[i - 1][j][k]) / dr;
                double dS_dz = (scalar[i][j][k] - scalar[i][j][k - 1]) / dz;
                double dS_dth = (scalar[i][j][k] - scalar[i][j_prev][k]) / dtheta;

                // Laplacian diffusion
                double lap_S = (scalar[i+1][j][k] - 2*scalar[i][j][k] + scalar[i-1][j][k])/(dr*dr)
                              + (scalar[i][j][k+1] - 2*scalar[i][j][k] + scalar[i][j][k-1])/(dz*dz)
                              + (scalar[i][j_prev][k] - 2*scalar[i][j][k] + scalar[i][(j+1)%NTH][k])/(dtheta*dtheta*r*r);

                new_scalar[i][j][k] -= static_cast<float>(dt_advect *
                (
                    u[i][j][k] * dS_dr +
                    (v_theta[i][j][k] / r) * dS_dth +
                    w[i][j][k] * dS_dz
                ) - dt_advect * kappa * lap_S);
            }
        }
    }

    scalar = new_scalar;
}

/*This function advects the tracer field.
Takes in the timestep and advects the tracer field.*/
void advect_tracer(double dt_advect)
{
    advect_scalar(tracer, dt_advect, 0.01);  // tracer with diffusion
}

/*This function advects the thermodynamics.
Takes in the timestep and advects the thermodynamics.*/
void advect_thermodynamics(double dt_advect)
{
    advect_scalar(theta, dt_advect, 0.0);  // potential temperature (no diffusion)
    advect_scalar(qv, dt_advect, 0.0);     // water vapor (no diffusion)
    advect_scalar(qc, dt_advect, 0.0);     // cloud water (no diffusion)
    advect_scalar(qr, dt_advect, 0.0);     // rain water (no diffusion)
    advect_scalar(qi, dt_advect, 0.0);     // cloud ice (no diffusion)
    advect_scalar(qs, dt_advect, 0.0);     // snow (no diffusion)
    advect_scalar(qh, dt_advect, 0.0);     // hail (no diffusion)
    advect_scalar(qg, dt_advect, 0.0);     // graupel (no diffusion)
}

/*This function steps the microphysics forward in time.
Takes in the timestep and steps the microphysics forward in time.*/
void step_microphysics(double dt_micro)
{
    // If the microphysics scheme is not set, initialize the microphysics scheme.
    if (!microphysics_scheme) 
    {
        std::cerr << "Warning: Microphysics scheme not initialized, using default Kessler" << std::endl;
        initialize_microphysics("kessler");
    }

    // Temporary arrays for tendencies
    std::vector<std::vector<std::vector<float>>> dtheta_dt;
    std::vector<std::vector<std::vector<float>>> dqv_dt;
    std::vector<std::vector<std::vector<float>>> dqc_dt;
    std::vector<std::vector<std::vector<float>>> dqr_dt;
    std::vector<std::vector<std::vector<float>>> dqi_dt;
    std::vector<std::vector<std::vector<float>>> dqs_dt;
    std::vector<std::vector<std::vector<float>>> dqg_dt;
    std::vector<std::vector<std::vector<float>>> dqh_dt;

    // Compute microphysics tendencies using the selected scheme
    microphysics_scheme->compute_tendencies(
        p, theta, qv, qc, qr, qi, qs, qg, qh,
        dt_micro,
        dtheta_dt, dqv_dt, dqc_dt, dqr_dt, dqi_dt, dqs_dt, dqg_dt, dqh_dt
    );

    // Iterate over the rows, columns, and levels and apply the microphysics tendencies.
    for (int i = 0; i < NR; ++i)
    {
        // Iterate over the columns and apply the microphysics tendencies.
        for (int j = 0; j < NTH; ++j)
        {
            // Iterate over the levels and apply the microphysics tendencies.
    
            for (int k = 0; k < NZ; ++k)
            {
                // Apply microphysics + radiation + boundary layer tendencies
                theta[i][j][k] += (dtheta_dt[i][j][k] + dtheta_dt_rad[i][j][k] + dtheta_dt_pbl[i][j][k]) * dt_micro;
                qv[i][j][k] += (dqv_dt[i][j][k] + dqv_dt_pbl[i][j][k]) * dt_micro;
                qc[i][j][k] += dqc_dt[i][j][k] * dt_micro;
                qr[i][j][k] += dqr_dt[i][j][k] * dt_micro;
                qi[i][j][k] += dqi_dt[i][j][k] * dt_micro;
                qs[i][j][k] += dqs_dt[i][j][k] * dt_micro;
                qg[i][j][k] += dqg_dt[i][j][k] * dt_micro;
                qh[i][j][k] += dqh_dt[i][j][k] * dt_micro;

                // Ensure non-negative values
                qv[i][j][k] = std::max(0.0f, qv[i][j][k]);
                qc[i][j][k] = std::max(0.0f, qc[i][j][k]);
                qr[i][j][k] = std::max(0.0f, qr[i][j][k]);
                qi[i][j][k] = std::max(0.0f, qi[i][j][k]);
                qs[i][j][k] = std::max(0.0f, qs[i][j][k]);
                qg[i][j][k] = std::max(0.0f, qg[i][j][k]);
                qh[i][j][k] = std::max(0.0f, qh[i][j][k]);
            }
        }
    }
}

// === Nested Grid Functions ===

/*This function initializes the nested grid.
Takes in the nested grid and initializes the nested grid.*/
void initialize_nested_grid()
{
    // If the nested grid is not enabled, return.
    if (!nested_config.enabled) return;

    int nr_nest = nested_config.nest_size_r;
    int nth_nest = nested_config.nest_size_th;
    int nz_nest = nested_config.nest_size_z;

    // Resize nested grid fields
    nest_rho.assign(nr_nest, std::vector<std::vector<float>>(nth_nest, std::vector<float>(nz_nest, 0.0f)));
    nest_p.assign(nr_nest, std::vector<std::vector<float>>(nth_nest, std::vector<float>(nz_nest, 0.0f)));
    nest_u.assign(nr_nest, std::vector<std::vector<float>>(nth_nest, std::vector<float>(nz_nest, 0.0f)));
    nest_w.assign(nr_nest, std::vector<std::vector<float>>(nth_nest, std::vector<float>(nz_nest, 0.0f)));
    nest_v_theta.assign(nr_nest, std::vector<std::vector<float>>(nth_nest, std::vector<float>(nz_nest, 0.0f)));
    nest_theta.assign(nr_nest, std::vector<std::vector<float>>(nth_nest, std::vector<float>(nz_nest, 0.0f)));
    nest_qv.assign(nr_nest, std::vector<std::vector<float>>(nth_nest, std::vector<float>(nz_nest, 0.0f)));
    nest_qc.assign(nr_nest, std::vector<std::vector<float>>(nth_nest, std::vector<float>(nz_nest, 0.0f)));
    nest_qr.assign(nr_nest, std::vector<std::vector<float>>(nth_nest, std::vector<float>(nz_nest, 0.0f)));
    nest_qh.assign(nr_nest, std::vector<std::vector<float>>(nth_nest, std::vector<float>(nz_nest, 0.0f)));
    nest_qg.assign(nr_nest, std::vector<std::vector<float>>(nth_nest, std::vector<float>(nz_nest, 0.0f)));
    nest_tracer.assign(nr_nest, std::vector<std::vector<float>>(nth_nest, std::vector<float>(nz_nest, 0.0f)));

    // Interpolate from parent grid to initialize nested grid
    // This is a simplified bilinear/trilinear interpolation
    double ref = nested_config.refinement;
    int ci = nested_config.center_i;
    int cj = nested_config.center_j;
    int ck = nested_config.center_k;

    // Iterate over the nested grid and interpolate from the parent grid to initialize the nested grid.
    for (int i_nest = 0; i_nest < nr_nest; ++i_nest)
    {
        // Iterate over the columns and interpolate from the parent grid to initialize the nested grid.
        for (int j_nest = 0; j_nest < nth_nest; ++j_nest)
        {
            // Iterate over the levels and interpolate from the parent grid to initialize the nested grid.
            for (int k_nest = 0; k_nest < nz_nest; ++k_nest)
            {
                // Map nested grid coordinates to parent grid coordinates
                double i_parent = ci + (i_nest - nr_nest/2.0) / ref;
                double j_parent = cj + (j_nest - nth_nest/2.0) / ref;
                double k_parent = ck + (k_nest - nz_nest/2.0) / ref;

                // Bilinear interpolation (simplified - should use proper interpolation)
                int i0 = std::max(0, std::min(NR-2, (int)std::floor(i_parent)));
                int j0 = std::max(0, std::min(NTH-2, (int)std::floor(j_parent)));
                int k0 = std::max(0, std::min(NZ-2, (int)std::floor(k_parent)));

                double fi = i_parent - i0;
                double fj = j_parent - j0;
                double fk = k_parent - k0;

                // Interpolate each field
                auto interpolate = [&](const auto& field, int i0, int j0, int k0, double fi, double fj, double fk) {
                    double v000 = field[i0][j0][k0];
                    double v001 = field[i0][j0][k0+1];
                    double v010 = field[i0][j0+1][k0];
                    double v011 = field[i0][j0+1][k0+1];
                    double v100 = field[i0+1][j0][k0];
                    double v101 = field[i0+1][j0][k0+1];
                    double v110 = field[i0+1][j0+1][k0];
                    double v111 = field[i0+1][j0+1][k0+1];

                    return v000 * (1-fi)*(1-fj)*(1-fk) +
                           v001 * (1-fi)*(1-fj)*fk +
                           v010 * (1-fi)*fj*(1-fk) +
                           v011 * (1-fi)*fj*fk +
                           v100 * fi*(1-fj)*(1-fk) +
                           v101 * fi*(1-fj)*fk +
                           v110 * fi*fj*(1-fk) +
                           v111 * fi*fj*fk;
                };

                nest_rho[i_nest][j_nest][k_nest] = interpolate(rho, i0, j0, k0, fi, fj, fk);
                nest_p[i_nest][j_nest][k_nest] = interpolate(p, i0, j0, k0, fi, fj, fk);
                nest_u[i_nest][j_nest][k_nest] = interpolate(u, i0, j0, k0, fi, fj, fk);
                nest_w[i_nest][j_nest][k_nest] = interpolate(w, i0, j0, k0, fi, fj, fk);
                nest_v_theta[i_nest][j_nest][k_nest] = interpolate(v_theta, i0, j0, k0, fi, fj, fk);
                nest_theta[i_nest][j_nest][k_nest] = interpolate(theta, i0, j0, k0, fi, fj, fk);
                nest_qv[i_nest][j_nest][k_nest] = interpolate(qv, i0, j0, k0, fi, fj, fk);
                nest_qc[i_nest][j_nest][k_nest] = interpolate(qc, i0, j0, k0, fi, fj, fk);
                nest_qr[i_nest][j_nest][k_nest] = interpolate(qr, i0, j0, k0, fi, fj, fk);
                nest_qh[i_nest][j_nest][k_nest] = interpolate(qh, i0, j0, k0, fi, fj, fk);
                nest_qg[i_nest][j_nest][k_nest] = interpolate(qg, i0, j0, k0, fi, fj, fk);
                nest_tracer[i_nest][j_nest][k_nest] = interpolate(tracer, i0, j0, k0, fi, fj, fk);
            }
        }
    }
}

/*This function feedbacks the nested grid to the parent grid.
Takes in the nested grid and feedbacks the nested grid to the parent grid.*/
void feedback_to_parent()
{
    if (!nested_config.enabled) return;

    double ref = nested_config.refinement;
    int ci = nested_config.center_i;
    int cj = nested_config.center_j;
    int ck = nested_config.center_k;
    int nr_nest = nested_config.nest_size_r;
    int nth_nest = nested_config.nest_size_th;
    int nz_nest = nested_config.nest_size_z;

    // Simple replacement of parent grid values in the nested region
    // (More sophisticated methods would use weighted averaging)

    // Iterate over the nested grid and feedback the nested grid to the parent grid.
    for (int i_nest = 0; i_nest < nr_nest; ++i_nest)
    {
        // Iterate over the columns and feedback the nested grid to the parent grid.
        for (int j_nest = 0; j_nest < nth_nest; ++j_nest)   
        {
            // Iterate over the levels and feedback the nested grid to the parent grid.
            for (int k_nest = 0; k_nest < nz_nest; ++k_nest)
            {
                // Map nested coordinates to parent coordinates
                int i_parent = ci + (int)std::round((i_nest - nr_nest/2.0) / ref);
                int j_parent = cj + (int)std::round((j_nest - nth_nest/2.0) / ref);
                int k_parent = ck + (int)std::round((k_nest - nz_nest/2.0) / ref);

                // Check bounds
                if (i_parent >= 0 && i_parent < NR &&
                    j_parent >= 0 && j_parent < NTH &&
                    k_parent >= 0 && k_parent < NZ)
                {
                    rho[i_parent][j_parent][k_parent] = nest_rho[i_nest][j_nest][k_nest];
                    p[i_parent][j_parent][k_parent] = nest_p[i_nest][j_nest][k_nest];
                    u[i_parent][j_parent][k_parent] = nest_u[i_nest][j_nest][k_nest];
                    w[i_parent][j_parent][k_parent] = nest_w[i_nest][j_nest][k_nest];
                    v_theta[i_parent][j_parent][k_parent] = nest_v_theta[i_nest][j_nest][k_nest];
                    theta[i_parent][j_parent][k_parent] = nest_theta[i_nest][j_nest][k_nest];
                    qv[i_parent][j_parent][k_parent] = nest_qv[i_nest][j_nest][k_nest];
                    qc[i_parent][j_parent][k_parent] = nest_qc[i_nest][j_nest][k_nest];
                    qr[i_parent][j_parent][k_parent] = nest_qr[i_nest][j_nest][k_nest];
                    qh[i_parent][j_parent][k_parent] = nest_qh[i_nest][j_nest][k_nest];
                    qg[i_parent][j_parent][k_parent] = nest_qg[i_nest][j_nest][k_nest];
                    tracer[i_parent][j_parent][k_parent] = nest_tracer[i_nest][j_nest][k_nest];
                }
            }
        }
    }
}


/*This function calculates the radar reflectivity.
Takes in the radar reflectivity and calculates the radar reflectivity.*/
void calculate_radar_reflectivity()
{
    // If the radar scheme is not set, return.
    if (!radar_scheme) 
    {
        std::cerr << "Warning: Radar scheme not initialized, using microphysics fallback" << std::endl;
        // Fallback to microphysics scheme for backward compatibility
        if (microphysics_scheme) {
            microphysics_scheme->compute_radar_reflectivity(
                qc, qr, qi, qs, qg, qh, radar_reflectivity
            );
        }
        return;
    }

    // Create radar state view with current model state
    RadarStateView state_view;
    state_view.NR = NR;
    state_view.NTH = NTH;
    state_view.NZ = NZ;

    // Winds
    state_view.u = &u;
    state_view.v = &v_theta;  // azimuthal wind
    state_view.w = &w;

    // Hydrometeor mixing ratios
    state_view.qr = &qr;
    state_view.qs = &qs;
    state_view.qg = &qg;
    state_view.qh = &qh;
    state_view.qi = &qi;

    // Configure radar
    RadarConfig config;
    config.scheme_id = "reflectivity";
    config.operator_tier = "fast_da";
    config.has_qr = true;
    config.has_qs = true;
    config.has_qg = true;
    config.has_qh = true;
    config.has_qi = true;

    // Output structure
    RadarOut radar_out;
    radar_out.initialize(NR, NTH, NZ);

    // Compute radar observables
    radar_scheme->compute(config, state_view, radar_out);

    // Copy reflectivity to the global radar_reflectivity field
    // (radar_out.Z_dBZ contains dBZ values, but we need linear for now)

    // Iterate over the rows, columns, and levels and copy the radar reflectivity to the global radar_reflectivity field.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over the columns and copy the radar reflectivity to the global radar_reflectivity field.
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over the levels and copy the radar reflectivity to the global radar_reflectivity field.
            for (int k = 0; k < NZ; ++k) 
            {   
                // Convert back from dBZ to linear for compatibility
                float Z_dbz = radar_out.Z_dBZ[i][j][k];
                float Z_linear = (Z_dbz > -30.0f) ? std::pow(10.0f, Z_dbz / 10.0f) : 0.0f;
                radar_reflectivity[i][j][k] = Z_linear;
            }
        }
    }
}
