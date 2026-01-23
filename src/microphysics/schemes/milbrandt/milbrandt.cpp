#include "milbrandt.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

/*This constructor initializes the Milbrandt scheme with default parameters.
Takes in the autoconversion threshold, autoconversion rate, alpha parameters for the rain, ice, snow, graupel, and hail,
c parameters for the rain, ice, snow, graupel, and hail, d parameters for the rain, ice, snow, graupel, and hail,
a parameters for the rain, ice, snow, graupel, and hail, b parameters for the rain, ice, snow, graupel, and hail,
triple moment flag, and hail processes flag and initializes the Milbrandt scheme with the default parameters.*/
MilbrandtScheme::MilbrandtScheme(
    double qc0, double c_auto,
    double alpha_r, double alpha_i, double alpha_s, double alpha_g, double alpha_h,
    double c_r, double d_r, double c_i, double d_i, double c_s, double d_s,
    double c_g, double d_g, double c_h, double d_h,
    double a_r, double b_r, double a_i, double b_i, double a_s, double b_s,
    double a_g, double b_g, double a_h, double b_h,
    bool triple_moment, bool hail_processes
) : qc0_(qc0), c_auto_(c_auto),
    alpha_r_(alpha_r), alpha_i_(alpha_i), alpha_s_(alpha_s), alpha_g_(alpha_g), alpha_h_(alpha_h),
    c_r_(c_r), d_r_(d_r), c_i_(c_i), d_i_(d_i), c_s_(c_s), d_s_(d_s),
    c_g_(c_g), d_g_(d_g), c_h_(c_h), d_h_(d_h),
    a_r_(a_r), b_r_(b_r), a_i_(a_i), b_i_(b_i), a_s_(a_s), b_s_(b_s),
    a_g_(a_g), b_g_(b_g), a_h_(a_h), b_h_(b_h),
    triple_moment_(triple_moment), hail_processes_(hail_processes) {}

/*This function computes the tendencies for the Milbrandt scheme.
Takes in the pressure, potential temperature, vapor mixing ratio, cloud water mixing ratio, 
rainwater mixing ratio, ice mixing ratio, snow mixing ratio, graupel mixing ratio, hail mixing ratio,
and the time step and computes the tendencies for the Milbrandt scheme.*/
void MilbrandtScheme::compute_tendencies(
    const std::vector<std::vector<std::vector<float>>>& p,
    const std::vector<std::vector<std::vector<float>>>& theta,
    const std::vector<std::vector<std::vector<float>>>& qv,
    const std::vector<std::vector<std::vector<float>>>& qc,
    const std::vector<std::vector<std::vector<float>>>& qr,
    const std::vector<std::vector<std::vector<float>>>& qi,
    const std::vector<std::vector<std::vector<float>>>& qs,
    const std::vector<std::vector<std::vector<float>>>& qg,
    const std::vector<std::vector<std::vector<float>>>& qh,
    double dt,
    std::vector<std::vector<std::vector<float>>>& dtheta_dt,
    std::vector<std::vector<std::vector<float>>>& dqv_dt,
    std::vector<std::vector<std::vector<float>>>& dqc_dt,
    std::vector<std::vector<std::vector<float>>>& dqr_dt,
    std::vector<std::vector<std::vector<float>>>& dqi_dt,
    std::vector<std::vector<std::vector<float>>>& dqs_dt,
    std::vector<std::vector<std::vector<float>>>& dqg_dt,
    std::vector<std::vector<std::vector<float>>>& dqh_dt
) 

{
    // Get the number of rows, columns, and levels.
    int NR = p.size();
    if (NR == 0) return;
    int NTH = p[0].size();
    int NZ = p[0][0].size();

    // Initialize prognostic fields if not already done
    if (Nr_.empty()) 
    {
        initialize_prognostic_fields(NR, NTH, NZ);
    }

    // Convert theta to temperature
    std::vector<std::vector<std::vector<float>>> temperature;
    thermodynamics::convert_theta_to_temperature_field(theta, p, temperature);

    // Calculate air density
    std::vector<std::vector<std::vector<float>>> rho(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));

    // Iterate over all grid points to calculate the air density.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal points
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical points
            for (int k = 0; k < NZ; ++k) 
            {
                // Calculate the air density.
                rho[i][j][k] = p[i][j][k] / (287.0 * temperature[i][j][k] * (1.0 + 0.608 * qv[i][j][k]));
            }
        }
    }

    // Initialize tendency arrays
    dtheta_dt.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dqv_dt.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dqc_dt.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dqr_dt.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dqi_dt.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dqs_dt.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dqg_dt.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dqh_dt.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));

    dNr_dt_.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dNi_dt_.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dNs_dt_.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dNg_dt_.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dNh_dt_.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));

    // Compute microphysical processes
    compute_warm_rain_processes(temperature, p, rho, qv, qc, qr,
                               dqc_dt, dqr_dt, dqv_dt, dtheta_dt, dNr_dt_);

    compute_ice_processes(temperature, p, rho, qv, qc, qi, qs, qg, qh,
                         dqc_dt, dqi_dt, dqs_dt, dqg_dt, dqh_dt, dqv_dt, dtheta_dt,
                         dNi_dt_, dNs_dt_, dNg_dt_, dNh_dt_);

    compute_melting_processes(temperature, qs, qg, qh,
                            dqs_dt, dqg_dt, dqh_dt, dqr_dt, dtheta_dt, dNr_dt_,
                            dNs_dt_, dNg_dt_, dNh_dt_);

    compute_sedimentation(qr, qs, qg, qh, Nr_, Ns_, Ng_, Nh_,
                         dqr_dt, dqs_dt, dqg_dt, dqh_dt,
                         dNr_dt_, dNs_dt_, dNg_dt_, dNh_dt_);

    // Iterate over all grid points and update prognostic number concentrations
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal points
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical points
            for (int k = 0; k < NZ; ++k) 
            {
                // Update prognostic number concentrations
                Nr_[i][j][k] += dNr_dt_[i][j][k] * dt;
                Ni_[i][j][k] += dNi_dt_[i][j][k] * dt;
                Ns_[i][j][k] += dNs_dt_[i][j][k] * dt;
                Ng_[i][j][k] += dNg_dt_[i][j][k] * dt;
                Nh_[i][j][k] += dNh_dt_[i][j][k] * dt;

                // Ensure positive values
                Nr_[i][j][k] = std::max(Nr_[i][j][k], 1.0f);
                Ni_[i][j][k] = std::max(Ni_[i][j][k], 1.0f);
                Ns_[i][j][k] = std::max(Ns_[i][j][k], 1.0f);
                Ng_[i][j][k] = std::max(Ng_[i][j][k], 1.0f);
                Nh_[i][j][k] = std::max(Nh_[i][j][k], 1.0f);
            }
        }
    }

    // Iterate over all grid points and convert temperature tendencies to potential temperature tendencies
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal points 
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical points
            for (int k = 0; k < NZ; ++k) 
            {
                // Convert temperature tendencies to potential temperature tendencies
                dtheta_dt[i][j][k] = static_cast<float>(
                    thermodynamics::temperature_tendency_to_theta(
                        dtheta_dt[i][j][k], theta[i][j][k], p[i][j][k]
                    )
                );
            }
        }
    }
}

/*This function initializes the prognostic number concentrations for the Milbrandt scheme.
Takes in the number of rows, columns, and levels and initializes the prognostic number concentrations for the Milbrandt scheme.*/
void MilbrandtScheme::initialize_prognostic_fields(int NR, int NTH, int NZ) 
{
    // Initialize number concentrations with reasonable defaults
    Nr_.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 1e6)));   // rain
    Ni_.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 1e5)));   // ice
    Ns_.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 1e5)));   // snow
    Ng_.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 1e4)));   // graupel
    Nh_.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 1e3)));   // hail

    // If triple moment is enabled, initialize the triple moment fields
    if (triple_moment_) 
    {
        Zr_.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 1.0)));
        Zi_.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 1.0)));
        Zs_.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 1.0)));
        Zg_.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 1.0)));
        Zh_.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 1.0)));
    }
}

/*This function calculates the lambda parameter for the Milbrandt scheme.
Takes in the hydrometeor mixing ratio, number concentration, alpha parameter, c parameter, and d parameter
and calculates the lambda parameter for the Milbrandt scheme.*/
double MilbrandtScheme::calculate_lambda(double q, double Nt, double alpha, double c, double d) 
{
    if (q <= 0.0 || Nt <= 0.0) return 1.0;

    double rho = 1.0; // air density, simplified
    double lambda = pow(gamma_function(1.0 + d + alpha) /
                       (gamma_function(1.0 + alpha) * c * Nt / (rho * q)), 1.0 / d);
    return std::max(lambda, 1.0);
}

/* This is the calculate N0 function for the Milbrandt scheme
// It calculates the N0 parameter for the Milbrandt scheme
// by iterating over all grid points
// and calculating the N0 parameter using the calculate_N0 function*/
double MilbrandtScheme::calculate_N0(double Nt, double alpha, double lambda) 
{
    return Nt * pow(lambda, 1.0 + alpha) / gamma_function(1.0 + alpha);
}

/* This is the calculate moment function for the Milbrandt scheme
// It calculates the moment parameter for the Milbrandt scheme
// by iterating over all grid points
// and calculating the moment parameter using the calculate_moment function*/
double MilbrandtScheme::calculate_moment(double q, double Nt, double alpha, double c, double d, int moment_order) 
{
    double lambda = calculate_lambda(q, Nt, alpha, c, d);
    return Nt * pow(lambda, -moment_order) * gamma_function(1.0 + alpha + moment_order);
}

/* This is the calculate reflectivity function for the Milbrandt scheme
// It calculates the reflectivity parameter for the Milbrandt scheme
// by iterating over all grid points
// and calculating the reflectivity parameter using the calculate_reflectivity function*/
double MilbrandtScheme::calculate_reflectivity(double q, double Nt, double alpha, double c, double d) 
{
    // Z = integral N(D) * D^6 dD
    double moment6 = calculate_moment(q, Nt, alpha, c, d, 6);
    return moment6 * 1e-18; // Convert to mm^6/m^3
}


/* This is the compute warm rain processes function for the Milbrandt scheme
// It computes the warm rain processes by iterating over all grid points
// and computing the warm rain processes using the compute_warm_rain_processes function*/
void MilbrandtScheme::compute_warm_rain_processes(
    const std::vector<std::vector<std::vector<float>>>& temperature,
    const std::vector<std::vector<std::vector<float>>>& p,
    const std::vector<std::vector<std::vector<float>>>& rho,
    const std::vector<std::vector<std::vector<float>>>& qv,
    const std::vector<std::vector<std::vector<float>>>& qc,
    const std::vector<std::vector<std::vector<float>>>& qr,
    std::vector<std::vector<std::vector<float>>>& dqc_dt,
    std::vector<std::vector<std::vector<float>>>& dqr_dt,
    std::vector<std::vector<std::vector<float>>>& dqv_dt,
    std::vector<std::vector<std::vector<float>>>& dtheta_dt,
    std::vector<std::vector<std::vector<float>>>& dNr_dt
) 
{
    int NR = temperature.size();
    int NTH = temperature[0].size();
    int NZ = temperature[0][0].size();

    // Iterate over all grid points and compute warm rain processes
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal points
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical points
            for (int k = 0; k < NZ; ++k) 
            {
                // Get hydrometeor values
                float qc_val = qc[i][j][k];
                float qr_val = qr[i][j][k];
                float qv_val = qv[i][j][k];
                float T = temperature[i][j][k];
                float P = p[i][j][k];
                float rho_val = rho[i][j][k];
                float Nr_val = Nr_[i][j][k];

                // If the cloud water mixing ratio is greater than the autoconversion threshold, compute the autoconversion rate.
                if (qc_val > qc0_) 
                {
                    // Compute autoconversion rate
                    double auto_rate = autoconversion_rate(qc_val, rho_val);
                    dqc_dt[i][j][k] -= auto_rate;
                    dqr_dt[i][j][k] += auto_rate;
                    dtheta_dt[i][j][k] += microphysics_constants::L_v / microphysics_constants::cp * auto_rate / T;

                    // Number concentration from autoconversion
                    dNr_dt[i][j][k] += auto_rate / (c_r_ * pow(1e-3, d_r_)); // Rough estimate
                }

                // If the cloud water mixing ratio is greater than 0 and the rainwater mixing ratio is greater than 0, compute the accretion rate.
                if (qc_val > 0.0f && qr_val > 0.0f) 
                {
                    // Compute accretion rate
                    double accr_rate = accretion_rate(qc_val, qr_val, 1e8, Nr_val,
                                                    0.0, alpha_r_, c_r_, d_r_, c_r_, d_r_, rho_val);
                    dqc_dt[i][j][k] -= accr_rate;
                    dqr_dt[i][j][k] += accr_rate;
                    dtheta_dt[i][j][k] += microphysics_constants::L_v / microphysics_constants::cp * accr_rate / T;
                }

                // If the rainwater number concentration is greater than 0, compute the rain self-collection rate.
                if (Nr_val > 0.0f) 
                {
                    double self_coll = rain_selfcollection_rate(Nr_val, rho_val);
                    dNr_dt[i][j][k] -= self_coll;
                }

                // If the rainwater mixing ratio is greater than 0, compute the evaporation rate.
                if (qr_val > 0.0f) 
                {
                    // Compute saturation mixing ratio
                    float qvsat = thermodynamics::saturation_mixing_ratio_water(T, P);
                    float RH = (qv_val > 0.0f) ? qv_val / qvsat : 0.0f;

                    // If the relative humidity is less than 1, compute the evaporation rate.
                    if (RH < 1.0f) 
                    {
                        double evap_rate = 1.0e-3 * (1.0f - RH) * qr_val;
                        dqv_dt[i][j][k] += evap_rate;
                        dqr_dt[i][j][k] -= evap_rate;
                        dtheta_dt[i][j][k] -= microphysics_constants::L_v / microphysics_constants::cp * evap_rate / T;

                        // Number concentration loss from evaporation
                        dNr_dt[i][j][k] -= evap_rate / (c_r_ * pow(1e-3, d_r_));
                    }
                }
            }
        }
    }
}

/*This function computes the ice processes for the Milbrandt scheme.
Takes in the temperature, pressure, vapor mixing ratio, cloud water mixing ratio, 
ice mixing ratio, snow mixing ratio, graupel mixing ratio, hail mixing ratio,
and the tendencies for the cloud water, ice, snow, graupel, hail, and vapor mixing ratios
and computes the ice processes for the Milbrandt scheme.*/
void MilbrandtScheme::compute_ice_processes(
    const std::vector<std::vector<std::vector<float>>>& temperature,
    const std::vector<std::vector<std::vector<float>>>& p,
    const std::vector<std::vector<std::vector<float>>>& rho,
    const std::vector<std::vector<std::vector<float>>>& qv,
    const std::vector<std::vector<std::vector<float>>>& qc,
    const std::vector<std::vector<std::vector<float>>>& qi,
    const std::vector<std::vector<std::vector<float>>>& qs,
    const std::vector<std::vector<std::vector<float>>>& qg,
    const std::vector<std::vector<std::vector<float>>>& qh,
    std::vector<std::vector<std::vector<float>>>& dqc_dt,
    std::vector<std::vector<std::vector<float>>>& dqi_dt,
    std::vector<std::vector<std::vector<float>>>& dqs_dt,
    std::vector<std::vector<std::vector<float>>>& dqg_dt,
    std::vector<std::vector<std::vector<float>>>& dqh_dt,
    std::vector<std::vector<std::vector<float>>>& dqv_dt,
    std::vector<std::vector<std::vector<float>>>& dtheta_dt,
    std::vector<std::vector<std::vector<float>>>& dNi_dt,
    std::vector<std::vector<std::vector<float>>>& dNs_dt,
    std::vector<std::vector<std::vector<float>>>& dNg_dt,
    std::vector<std::vector<std::vector<float>>>& dNh_dt
) 
{
    int NR = temperature.size();
    int NTH = temperature[0].size();
    int NZ = temperature[0][0].size();

    // Iterate over all grid points and compute ice processes
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal points
        for (int j = 0; j < NTH; ++j) 
        {
            for (int k = 0; k < NZ; ++k) 
            {
                // Get hydrometeor values
                float qc_val = qc[i][j][k];
                float qi_val = qi[i][j][k];
                float qs_val = qs[i][j][k];
                float qg_val = qg[i][j][k];
                float qh_val = qh[i][j][k];
                float qv_val = qv[i][j][k];
                float T = temperature[i][j][k];
                float P = p[i][j][k];
                float rho_val = rho[i][j][k];

                float Ni_val = Ni_[i][j][k];
                float Ns_val = Ns_[i][j][k];
                float Ng_val = Ng_[i][j][k];
                float Nh_val = Nh_[i][j][k];

                // If the cloud water mixing ratio is greater than 0 and the temperature is less than the freezing temperature 
                // and the ice number concentration is less than 1e6, compute the ice nucleation rate.(Simplified)
                
                if (qc_val > 0.0f && T < microphysics_constants::T0 && Ni_val < 1e6) 
                {
                    double nuc_rate = 1e-3 * qc_val;
                    dqc_dt[i][j][k] -= nuc_rate;
                    dqi_dt[i][j][k] += nuc_rate;
                    dNi_dt[i][j][k] += nuc_rate / (c_i_ * pow(1e-4, d_i_)); // Rough estimate
                }

                // If the ice mixing ratio is greater than 0, compute the deposition/sublimation rate.
                if (qi_val > 0.0f) 
                {
                    // Compute saturation mixing ratio
                    float qvsat_ice = thermodynamics::saturation_mixing_ratio_ice(T, P);

                    // If the vapor mixing ratio is greater than the saturation mixing ratio and the 
                    // temperature is less than the freezing temperature, compute the deposition rate.
                    if (qv_val > qvsat_ice && T < microphysics_constants::T0) 
                    {
                        // Compute deposition rate
                        // Deposition
                        double dep_rate = 1e-3 * (qv_val - qvsat_ice);
                        dqv_dt[i][j][k] -= dep_rate;
                        dqi_dt[i][j][k] += dep_rate;
                        dtheta_dt[i][j][k] += microphysics_constants::L_s / microphysics_constants::cp * dep_rate / T;
                    } 

                    // If the vapor mixing ratio is less than the saturation mixing ratio, compute the sublimation rate.
                    else if (qv_val < qvsat_ice) 
                    {
                        // Compute sublimation rate
                        // Sublimation
                        double subl_rate = 1e-3 * (qvsat_ice - qv_val);
                        dqv_dt[i][j][k] += subl_rate;
                        dqi_dt[i][j][k] -= subl_rate;
                        dtheta_dt[i][j][k] -= microphysics_constants::L_s / microphysics_constants::cp * subl_rate / T;
                    }
                }

                // If the cloud water mixing ratio is greater than 0 and the ice mixing ratio is greater than 0, compute the riming rate.
                if (qc_val > 0.0f && qi_val > 0.0f) 
                {
                    double rime_rate = accretion_rate(qc_val, qi_val, 1e8, Ni_val,
                                                    0.0, alpha_i_, c_r_, d_r_, c_i_, d_i_, rho_val);
                    dqc_dt[i][j][k] -= rime_rate;
                    dqg_dt[i][j][k] += rime_rate;  // Assume riming produces graupel
                    dtheta_dt[i][j][k] += microphysics_constants::L_f / microphysics_constants::cp * rime_rate / T;
                }

                // If the ice mixing ratio is greater than 1e-6 and the ice number concentration is greater than 1e3, 
                // compute the aggregation rate.(Simplified)
                if (qi_val > 1e-6 && Ni_val > 1e3) 
                {
                    double agg_rate = 1e-3 * qi_val;
                    dqi_dt[i][j][k] -= agg_rate;
                    dqs_dt[i][j][k] += agg_rate;
                    dNi_dt[i][j][k] -= agg_rate / (c_i_ * pow(1e-4, d_i_));
                    dNs_dt[i][j][k] += agg_rate / (c_s_ * pow(1e-3, d_s_));
                }
            }
        }
    }
}

/* This is the compute melting processes function for the Milbrandt scheme
// It computes the melting processes by iterating over all grid points
// and computing the melting processes using the compute_melting_processes function*/
void MilbrandtScheme::compute_melting_processes(
    const std::vector<std::vector<std::vector<float>>>& temperature,
    const std::vector<std::vector<std::vector<float>>>& qs,
    const std::vector<std::vector<std::vector<float>>>& qg,
    const std::vector<std::vector<std::vector<float>>>& qh,
    std::vector<std::vector<std::vector<float>>>& dqs_dt,
    std::vector<std::vector<std::vector<float>>>& dqg_dt,
    std::vector<std::vector<std::vector<float>>>& dqh_dt,
    std::vector<std::vector<std::vector<float>>>& dqr_dt,
    std::vector<std::vector<std::vector<float>>>& dtheta_dt,
    std::vector<std::vector<std::vector<float>>>& dNr_dt,
    std::vector<std::vector<std::vector<float>>>& dNs_dt,
    std::vector<std::vector<std::vector<float>>>& dNg_dt,
    std::vector<std::vector<std::vector<float>>>& dNh_dt
) 
{
    int NR = temperature.size();
    int NTH = temperature[0].size();
    int NZ = temperature[0][0].size();

    // Iterate over all grid points and compute melting processes
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal points
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical points
            for (int k = 0; k < NZ; ++k)
            {
                // Get hydrometeor values
                float qs_val = qs[i][j][k];
                float qg_val = qg[i][j][k];
                float qh_val = qh[i][j][k];
                float T = temperature[i][j][k];
                float Ns_val = Ns_[i][j][k];
                float Ng_val = Ng_[i][j][k];
                float Nh_val = Nh_[i][j][k];

                // If temperature is above freezing, compute melting processes
                if (T > microphysics_constants::T0) 
                {
                    // If the snow mixing ratio is greater than 0, compute the melting rate.
                    // Melting: snow/graupel/hail → rain
                    if (qs_val > 0.0f) 
                    {
                        // Compute melting rate
                        double melt_rate = 1e-3 * qs_val;
                        dqs_dt[i][j][k] -= melt_rate;
                        dqr_dt[i][j][k] += melt_rate;
                        dtheta_dt[i][j][k] -= microphysics_constants::L_f / microphysics_constants::cp * melt_rate / T;
                        dNr_dt[i][j][k] += melt_rate / (c_r_ * pow(1e-3, d_r_));
                        dNs_dt[i][j][k] -= melt_rate / (c_s_ * pow(1e-3, d_s_));
                    }

                    // If the graupel mixing ratio is greater than 0, compute the melting rate.
                    if (qg_val > 0.0f) 
                    {
                        // Compute melting rate
                        double melt_rate = 1e-3 * qg_val;
                        dqg_dt[i][j][k] -= melt_rate;
                        dqr_dt[i][j][k] += melt_rate;
                        dtheta_dt[i][j][k] -= microphysics_constants::L_f / microphysics_constants::cp * melt_rate / T;
                        dNr_dt[i][j][k] += melt_rate / (c_r_ * pow(1e-3, d_r_));
                        dNg_dt[i][j][k] -= melt_rate / (c_g_ * pow(1e-3, d_g_));
                    }

                    // If the hail mixing ratio is greater than 0 and the hail processes flag is true, compute the melting rate.
                    // Melting: hail → rain
                    if (qh_val > 0.0f && hail_processes_) 
                    {
                        // Compute melting rate
                        double melt_rate = 1e-3 * qh_val;
                        dqh_dt[i][j][k] -= melt_rate;
                        dqr_dt[i][j][k] += melt_rate;
                        dtheta_dt[i][j][k] -= microphysics_constants::L_f / microphysics_constants::cp * melt_rate / T;
                        dNr_dt[i][j][k] += melt_rate / (c_r_ * pow(1e-3, d_r_));
                        dNh_dt[i][j][k] -= melt_rate / (c_h_ * pow(1e-3, d_h_));
                    }
                }
            }
        }
    }
}

/*This function computes the sedimentation for the Milbrandt scheme.
Takes in the rainwater mixing ratio, snow mixing ratio, graupel mixing ratio, hail mixing ratio, 
and the tendencies for the rainwater, snow, graupel, and hail mixing ratios
and computes the sedimentation for the Milbrandt scheme.*/
void MilbrandtScheme::compute_sedimentation(
    const std::vector<std::vector<std::vector<float>>>& qr,
    const std::vector<std::vector<std::vector<float>>>& qs,
    const std::vector<std::vector<std::vector<float>>>& qg,
    const std::vector<std::vector<std::vector<float>>>& qh,
    const std::vector<std::vector<std::vector<float>>>& Nr,
    const std::vector<std::vector<std::vector<float>>>& Ns,
    const std::vector<std::vector<std::vector<float>>>& Ng,
    const std::vector<std::vector<std::vector<float>>>& Nh,
    std::vector<std::vector<std::vector<float>>>& dqr_dt,
    std::vector<std::vector<std::vector<float>>>& dqs_dt,
    std::vector<std::vector<std::vector<float>>>& dqg_dt,
    std::vector<std::vector<std::vector<float>>>& dqh_dt,
    std::vector<std::vector<std::vector<float>>>& dNr_dt,
    std::vector<std::vector<std::vector<float>>>& dNs_dt,
    std::vector<std::vector<std::vector<float>>>& dNg_dt,
    std::vector<std::vector<std::vector<float>>>& dNh_dt
) 
{
    // Get the number of rows, columns, and levels.
    int NR = qr.size();
    int NTH = qr[0].size();
    int NZ = qr[0][0].size();
    const double dz = 100.0;

    // Iterate over all grid points and compute sedimentation
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal points
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical points
            for (int k = 0; k < NZ; ++k) 
            {
                // If the rainwater mixing ratio is greater than 0 and the vertical level is greater than 0, compute the sedimentation rate.
                if (qr[i][j][k] > 0.0f && k > 0) 
                {
                    // Compute terminal velocity
                    double Vt = sedimentation_rate(qr[i][j][k], Nr[i][j][k], alpha_r_, c_r_, d_r_, dz);
                    dqr_dt[i][j][k] -= Vt;

                    // If the vertical level is less than the number of levels minus 1, compute the sedimentation rate for the next level.
                    if (k < NZ - 1) dqr_dt[i][j][k + 1] += Vt;
                    dNr_dt[i][j][k] -= Vt / (c_r_ * pow(1e-3, d_r_));

                    // If the vertical level is less than the number of levels minus 1, compute the number concentration for the next level.
                    if (k < NZ - 1) dNr_dt[i][j][k + 1] += Vt / (c_r_ * pow(1e-3, d_r_));
                }

                // If the snow mixing ratio is greater than 0 and the vertical level is greater than 0, compute the sedimentation rate.
                if (qs[i][j][k] > 0.0f && k > 0) 
                {
                    // Compute terminal velocity
                    double Vt = sedimentation_rate(qs[i][j][k], Ns[i][j][k], alpha_s_, c_s_, d_s_, dz);
                    dqs_dt[i][j][k] -= Vt;

                    // If the vertical level is less than the number of levels minus 1, compute the sedimentation rate for the next level.
                    if (k < NZ - 1) dqs_dt[i][j][k + 1] += Vt;
                    dNs_dt[i][j][k] -= Vt / (c_s_ * pow(1e-3, d_s_));

                    // If the vertical level is less than the number of levels minus 1, compute the number concentration for the next level.
                    if (k < NZ - 1) dNs_dt[i][j][k + 1] += Vt / (c_s_ * pow(1e-3, d_s_));
                }

                // If the graupel mixing ratio is greater than 0 and the vertical level is greater than 0, compute the sedimentation rate.
                if (qg[i][j][k] > 0.0f && k > 0) 
                {
                    // Compute terminal velocity
                    double Vt = sedimentation_rate(qg[i][j][k], Ng[i][j][k], alpha_g_, c_g_, d_g_, dz);
                    dqg_dt[i][j][k] -= Vt;
                    if (k < NZ - 1) dqg_dt[i][j][k + 1] += Vt;
                    dNg_dt[i][j][k] -= Vt / (c_g_ * pow(1e-3, d_g_));
                    if (k < NZ - 1) dNg_dt[i][j][k + 1] += Vt / (c_g_ * pow(1e-3, d_g_));
                }

                // If the hail mixing ratio is greater than 0 and the vertical level is greater than 0 and the hail processes flag is true, 
                // compute the sedimentation rate.
                if (qh[i][j][k] > 0.0f && k > 0 && hail_processes_) 
                {
                    double Vt = sedimentation_rate(qh[i][j][k], Nh[i][j][k], alpha_h_, c_h_, d_h_, dz);
                    dqh_dt[i][j][k] -= Vt;

                    // If the vertical level is less than the number of levels minus 1, compute the sedimentation rate for the next level.
                    if (k < NZ - 1) dqh_dt[i][j][k + 1] += Vt;
                    dNh_dt[i][j][k] -= Vt / (c_h_ * pow(1e-3, d_h_));

                    // If the vertical level is less than the number of levels minus 1, compute the number concentration for the next level.
                    if (k < NZ - 1) dNh_dt[i][j][k + 1] += Vt / (c_h_ * pow(1e-3, d_h_));
                }
            }
        }
    }
}

/*This function computes the radar reflectivity for the Milbrandt scheme.
Takes in the cloud water mixing ratio, rainwater mixing ratio, ice mixing ratio, snow mixing ratio, 
graupel mixing ratio, hail mixing ratio, and the radar reflectivity field
and computes the radar reflectivity for the Milbrandt scheme.*/
void MilbrandtScheme::compute_radar_reflectivity(
    const std::vector<std::vector<std::vector<float>>>& qc,
    const std::vector<std::vector<std::vector<float>>>& qr,
    const std::vector<std::vector<std::vector<float>>>& qi,
    const std::vector<std::vector<std::vector<float>>>& qs,
    const std::vector<std::vector<std::vector<float>>>& qg,
    const std::vector<std::vector<std::vector<float>>>& qh,
    std::vector<std::vector<std::vector<float>>>& reflectivity_dbz
) 
{
    // Get the number of rows, columns, and levels.
    int NR = qc.size();
    if (NR == 0) return;
    int NTH = qc[0].size();
    int NZ = qc[0][0].size();

    reflectivity_dbz.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));

    // Iterate over all grid points and compute radar reflectivity and use prognostic 
    // reflectivity if triple-moment, otherwise calculate
    for (int i = 0; i < NR; ++i) 
    {

        // Iterate over all azimuthal points
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical points
            for (int k = 0; k < NZ; ++k) 
            {
                // Initialize total reflectivity
                float Z_total = 0.0f;

                // If the rainwater mixing ratio is greater than 0 and the rainwater number concentration is greater than 0, 
                // compute the reflectivity for the rain.
                if (qr[i][j][k] > 0.0f && Nr_[i][j][k] > 0.0f) 
                {
                    Z_total += calculate_reflectivity(qr[i][j][k], Nr_[i][j][k], alpha_r_, c_r_, d_r_);
                }

                // If the ice mixing ratio is greater than 0 and the ice number concentration is greater than 0, 
                // compute the reflectivity for the ice.
                if (qi[i][j][k] > 0.0f && Ni_[i][j][k] > 0.0f) 
                {
                    Z_total += calculate_reflectivity(qi[i][j][k], Ni_[i][j][k], alpha_i_, c_i_, d_i_);
                }

                // If the snow mixing ratio is greater than 0 and the snow number concentration is greater than 0, 
                // compute the reflectivity for the snow.
                if (qs[i][j][k] > 0.0f && Ns_[i][j][k] > 0.0f) 
                {
                    Z_total += calculate_reflectivity(qs[i][j][k], Ns_[i][j][k], alpha_s_, c_s_, d_s_);
                }

                // If the graupel mixing ratio is greater than 0 and the graupel number concentration is greater than 0, 
                // compute the reflectivity for the graupel.
                if (qg[i][j][k] > 0.0f && Ng_[i][j][k] > 0.0f) 
                {
                    Z_total += calculate_reflectivity(qg[i][j][k], Ng_[i][j][k], alpha_g_, c_g_, d_g_);
                }

                // If the hail mixing ratio is greater than 0 and the hail number concentration is greater than 0 and 
                // the hail processes flag is true, compute the reflectivity for the hail.
                if (qh[i][j][k] > 0.0f && Nh_[i][j][k] > 0.0f && hail_processes_) 
                {
                    Z_total += calculate_reflectivity(qh[i][j][k], Nh_[i][j][k], alpha_h_, c_h_, d_h_);
                }

                // If the total reflectivity is greater than 1e-10, compute the reflectivity in dBZ.
                if (Z_total > 1e-10) 
                {
                    reflectivity_dbz[i][j][k] = 10.0f * std::log10(Z_total);
                } 
                else 
                {
                    reflectivity_dbz[i][j][k] = -30.0f; // Below threshold
                }
            }
        }
    }
}

// Utility functions

/*This function computes the autoconversion rate for the Milbrandt scheme.
Takes in the cloud water mixing ratio and air density and computes the autoconversion rate for the Milbrandt scheme.*/
double MilbrandtScheme::autoconversion_rate(double qc, double rho) 
{
    // Simplified autoconversion following Kessler-type
    return c_auto_ * std::max(qc - qc0_, 0.0);
}

/*This function computes the accretion rate for the Milbrandt scheme.
Takes in the first hydrometeor mixing ratio, second hydrometeor mixing ratio, first hydrometeor number concentration, 
second hydrometeor number concentration, first hydrometeor alpha parameter, second hydrometeor alpha parameter, 
first hydrometeor c parameter, first hydrometeor d parameter, second hydrometeor c parameter, second hydrometeor d parameter,
and air density and computes the accretion rate for the Milbrandt scheme.*/
double MilbrandtScheme::accretion_rate(double q1, double q2, double Nt1, double Nt2,
                                     double alpha1, double alpha2, double c1, double d1,
                                     double c2, double d2, double rho) 
{
    // Simplified accretion calculation
    return 2.2 * q1 * q2;  // Basic form
}

/*This function computes the rain self-collection rate for the Milbrandt scheme.
Takes in the rainwater number concentration and air density and computes the rain self-collection 
rate for the Milbrandt scheme.*/
double MilbrandtScheme::rain_selfcollection_rate(double Nr, double rho) 
{
    // Raindrop self-collection
    return 5.78e-4 * Nr * Nr / rho;
}

/*This function computes the sedimentation rate for the Milbrandt scheme.
Takes in the hydrometeor mixing ratio, hydrometeor number concentration, hydrometeor alpha parameter, 
hydrometeor c parameter, hydrometeor d parameter, and the vertical grid spacing and computes the sedimentation 
rate for the Milbrandt scheme.*/
double MilbrandtScheme::sedimentation_rate(double q, double Nt, double alpha, double c, double d, int dz) 
{
    // Calculate mass-weighted terminal velocity
    double lambda = calculate_lambda(q, Nt, alpha, c, d);
    double Vt_mean = a_r_ * pow(lambda, -b_r_);  // Simplified - should be more sophisticated
    return Vt_mean * q / dz;  // Flux per unit height
}

/*This function computes the gamma function for the Milbrandt scheme.
Takes in the x value and computes the gamma function for the Milbrandt scheme.*/
double MilbrandtScheme::gamma_function(double x) 
{
    // Simplified gamma function approximation for common values
    if (std::abs(x - 1.0) < 1e-6) return 1.0;
    if (std::abs(x - 2.0) < 1e-6) return 1.0;
    if (std::abs(x - 3.0) < 1e-6) return 2.0;
    if (std::abs(x - 4.0) < 1e-6) return 6.0;
    if (std::abs(x - 5.0) < 1e-6) return 24.0;
    if (std::abs(x - 6.0) < 1e-6) return 120.0;
    // For other values, use approximation
    return std::tgamma(x);
}

/*This function computes the incomplete gamma function for the Milbrandt scheme.
Takes in the a value and x value and computes the incomplete gamma function for 
the Milbrandt scheme.*/
double MilbrandtScheme::incomplete_gamma(double a, double x) 
{
    // Placeholder - would need proper implementation for triple-moment
    return gamma_function(a) * std::exp(-x);  // Approximation
}
