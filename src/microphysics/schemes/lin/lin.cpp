#include "lin.hpp"
#include <algorithm>
#include <cmath>


/* This is the constructor for the Lin scheme. This is simplified for the purpose of this project
// In a full implementation, this would be a more complex function
// For example, the parameters would be configurable
// and the scheme would be more complex
// For example, the accretion rate would be a function of the cloud water and rain water mixing ratios
// and the evaporation rate would be a function of the rain water mixing ratio and the saturation mixing ratio
// and the homogeneous freezing rate would be a function of the cloud water mixing ratio
// and the riming rate would be a function of the cloud water and ice mixing ratios*/
LinScheme::LinScheme(
    double qc0, double c_auto, double c_accr, double c_evap,
    double c_ihom, double c_rime, double c_agg, double c_melt, double c_subl,
    double N0r, double N0s, double N0g, double N0h,
    double a_r, double b_r, double a_s, double b_s, double a_g, double b_g, double a_h, double b_h
) : qc0_(qc0), c_auto_(c_auto), c_accr_(c_accr), c_evap_(c_evap),
    c_ihom_(c_ihom), c_rime_(c_rime), c_agg_(c_agg), c_melt_(c_melt), c_subl_(c_subl),
    N0r_(N0r), N0s_(N0s), N0g_(N0g), N0h_(N0h),
    a_r_(a_r), b_r_(b_r), a_s_(a_s), b_s_(b_s),
    a_g_(a_g), b_g_(b_g), a_h_(a_h), b_h_(b_h) {}

/*This function computes the tendencies for the Lin scheme.
Takes in the pressure, potential temperature, vapor mixing ratio, 
cloud water mixing ratio, rainwater mixing ratio, ice mixing ratio, 
snow mixing ratio, graupel mixing ratio, hail mixing ratio,
and the time step and computes the tendencies for the Lin scheme.*/
void LinScheme::compute_tendencies(
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
    // Get grid size
    int NR = p.size();
    if (NR == 0) return;
    int NTH = p[0].size();
    int NZ = p[0][0].size();

    // Convert theta to temperature for microphysics calculations by iterating over all grid points
    std::vector<std::vector<std::vector<float>>> temperature;
    thermodynamics::convert_theta_to_temperature_field(theta, p, temperature);

    // Initialize tendency arrays to zero
    dtheta_dt.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dqv_dt.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dqc_dt.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dqr_dt.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dqi_dt.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dqs_dt.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dqg_dt.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));
    dqh_dt.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));

    // Saturation adjustment (cloud condensation)
    auto qv_temp = qv;
    auto qc_temp = qc;
    saturation_adjustment(temperature, p, qv_temp, qc_temp);

    // Compute microphysical processes
    compute_warm_rain_processes(temperature, p, qv_temp, qc_temp, qr, dqc_dt, dqr_dt, dqv_dt, dtheta_dt);
    compute_ice_processes(temperature, p, qv_temp, qc_temp, qi, qs, qg, qh,
                         dqc_dt, dqi_dt, dqs_dt, dqg_dt, dqh_dt, dqv_dt, dtheta_dt);
    compute_melting_processes(temperature, qs, qg, qh, dqs_dt, dqg_dt, dqh_dt, dqr_dt, dtheta_dt);
    compute_sedimentation(qr, qs, qg, qh, dqr_dt, dqs_dt, dqg_dt, dqh_dt);

    // Iterate over all grid points to convert the temperature tendencies to the potential temperature tendencies.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal points
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical points
            for (int k = 0; k < NZ; ++k) 
            {
                // Convert temperature tendency to potential temperature tendency
                dtheta_dt[i][j][k] = static_cast<float>(
                    thermodynamics::temperature_tendency_to_theta(
                        dtheta_dt[i][j][k], theta[i][j][k], p[i][j][k]
                    )
                );
            }
        }
    }
}


/*This function computes the saturation adjustment for the Lin scheme.
Takes in the temperature, pressure, vapor mixing ratio, and cloud water 
mixing ratio and computes the saturation adjustment for the Lin scheme.*/
void LinScheme::saturation_adjustment(
    const std::vector<std::vector<std::vector<float>>>& temperature,
    const std::vector<std::vector<std::vector<float>>>& p,
    std::vector<std::vector<std::vector<float>>>& qv,
    std::vector<std::vector<std::vector<float>>>& qc
) 
{
    // Get the number of rows, columns, and levels.
    int NR = temperature.size();
    int NTH = temperature[0].size();
    int NZ = temperature[0][0].size();

    // Iterate over all grid points to compute the saturation adjustment.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal points
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical points
            for (int k = 0; k < NZ; ++k) 
            {
                // Get the temperature, pressure, and vapor mixing ratio
                float T = temperature[i][j][k];
                float P = p[i][j][k];
                float qvsat = thermodynamics::saturation_mixing_ratio_water(T, P);

                // If the vapor mixing ratio is greater than the saturation vapor mixing ratio, compute the saturation adjustment.
                if (qv[i][j][k] > qvsat) 
                {
                    float excess = qv[i][j][k] - qvsat;
                    qc[i][j][k] += excess;
                    qv[i][j][k] = qvsat;
                }
            }
        }
    }
}


/*This function computes the warm rain processes for the Lin scheme.
Takes in the temperature, pressure, vapor mixing ratio, cloud water mixing ratio, 
rainwater mixing ratio, and the tendencies for the cloud water, rainwater, 
and vapor mixing ratios and computes the warm rain processes for the Lin scheme.*/
void LinScheme::compute_warm_rain_processes(
    const std::vector<std::vector<std::vector<float>>>& temperature,
    const std::vector<std::vector<std::vector<float>>>& p,
    const std::vector<std::vector<std::vector<float>>>& qv,
    const std::vector<std::vector<std::vector<float>>>& qc,
    const std::vector<std::vector<std::vector<float>>>& qr,
    std::vector<std::vector<std::vector<float>>>& dqc_dt,
    std::vector<std::vector<std::vector<float>>>& dqr_dt,
    std::vector<std::vector<std::vector<float>>>& dqv_dt,
    std::vector<std::vector<std::vector<float>>>& dtheta_dt
) 
{
    // Get grid size
    int NR = temperature.size();
    int NTH = temperature[0].size();
    int NZ = temperature[0][0].size();

    // Iterate over all grid points to compute the warm rain processes.
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

                // If the cloud water mixing ratio is greater than the critical cloud water mixing ratio, compute the autoconversion rate.
                // Autoconversion: qc → qr
                if (qc_val > qc0_) 
                {
                    // Compute autoconversion rate
                    float auto_rate = c_auto_ * (qc_val - qc0_);
                    dqc_dt[i][j][k] -= auto_rate;
                    dqr_dt[i][j][k] += auto_rate;
                    dtheta_dt[i][j][k] += microphysics_constants::L_v / microphysics_constants::cp * auto_rate / T;
                }

                // If the cloud water mixing ratio is greater than 0 and the rainwater mixing ratio is greater than 0, compute the accretion rate.
                // Accretion: qr collects qc
                if (qc_val > 0.0f && qr_val > 0.0f) 
                {
                    float accr_rate = c_accr_ * qc_val * qr_val;
                    dqc_dt[i][j][k] -= accr_rate;
                    dqr_dt[i][j][k] += accr_rate;
                    dtheta_dt[i][j][k] += microphysics_constants::L_v / microphysics_constants::cp * accr_rate / T;
                }

                // If the rainwater mixing ratio is greater than 0, compute the evaporation rate.
                // Evaporation: qr → qv
                if (qr_val > 0.0f) 
                {
                    // Compute saturation mixing ratio
                    float qvsat = thermodynamics::saturation_mixing_ratio_water(T, P);
                    float RH = (qv_val > 0.0f) ? qv_val / qvsat : 0.0f;

                    // If the relative humidity is less than 1, compute the evaporation rate.
                    if (RH < 1.0f) 
                    {
                        // Compute evaporation rate
                        float evap_rate = c_evap_ * (1.0f - RH) * qr_val;
                        dqv_dt[i][j][k] += evap_rate;
                        dqr_dt[i][j][k] -= evap_rate;
                        dtheta_dt[i][j][k] -= microphysics_constants::L_v / microphysics_constants::cp * evap_rate / T;
                    }
                }
            }
        }
    }
}


/*This function computes the ice processes for the Lin scheme.
Takes in the temperature, pressure, vapor mixing ratio, cloud water mixing ratio, 
rainwater mixing ratio, ice mixing ratio, snow mixing ratio, graupel mixing ratio, hail mixing ratio,
and the tendencies for the cloud water, graupel, hail, and vapor mixing ratios
and computes the ice processes for the Lin scheme.*/
void LinScheme::compute_ice_processes(
    const std::vector<std::vector<std::vector<float>>>& temperature,
    const std::vector<std::vector<std::vector<float>>>& p,
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
    std::vector<std::vector<std::vector<float>>>& dtheta_dt
) 
{
    // Get grid size
    int NR = temperature.size();
    int NTH = temperature[0].size();
    int NZ = temperature[0][0].size();

    // For each row, azimuthal point, and vertical point, compute the ice processes.
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
                float qi_val = qi[i][j][k];
                float qs_val = qs[i][j][k];
                float qg_val = qg[i][j][k];
                float qh_val = qh[i][j][k];
                float qv_val = qv[i][j][k];
                float T = temperature[i][j][k];
                float P = p[i][j][k];
                float T_celsius = T - microphysics_constants::T0;

                // If the temperature is less than the freezing temperature and the cloud water mixing ratio is greater than 0, compute the ice nucleation rate.
                if (T < microphysics_constants::T0 && qc_val > 0.0f) 
                {
                    // Compute ice nucleation rate
                    float nuc_rate = c_ihom_ * qc_val;
                    dqc_dt[i][j][k] -= nuc_rate;
                    dqi_dt[i][j][k] += nuc_rate;
                }

                // If the ice mixing ratio is greater than 0 and the temperature is less than the freezing temperature, compute the ice growth/deposition rate.
                // Ice growth/deposition
                if (qi_val > 0.0f && T < microphysics_constants::T0) 
                {
                    // Compute saturation mixing ratio
                    float qvsat_ice = thermodynamics::saturation_mixing_ratio_ice(T, P);

                    // If the vapor mixing ratio is greater than the saturation vapor mixing ratio, compute the deposition rate.
                    if (qv_val > qvsat_ice) 
                    {
                        // Compute deposition rate
                        float dep_rate = c_subl_ * (qv_val - qvsat_ice);
                        dqv_dt[i][j][k] -= dep_rate;
                        dqi_dt[i][j][k] += dep_rate;
                        dtheta_dt[i][j][k] += microphysics_constants::L_s / microphysics_constants::cp * dep_rate / T;
                    }
                }

                // Accretion of cloud water by ice (riming)
                double Ni = ice_number_concentration(qi_val, T_celsius);
                double E_ac = collection_efficiency(T_celsius);

                // If the cloud water mixing ratio is greater than 0 and the ice mixing ratio is greater than 0, compute the riming rate.
                if (qc_val > 0.0f && qi_val > 0.0f) 
                {
                    // Compute riming rate
                    // Riming to graupel (simplified)
                    double rime_rate = M_PI/4.0 * 0.1 * E_ac * Ni * qc_val;  // Simplified
                    dqc_dt[i][j][k] -= rime_rate;
                    dqg_dt[i][j][k] += rime_rate;
                    dtheta_dt[i][j][k] += microphysics_constants::L_f / microphysics_constants::cp * rime_rate / T;
                }

                // If the ice mixing ratio is greater than 0, compute the aggregation rate.
                // Aggregation: ice → snow
                if (qi_val > 0.0f) 
                {
                    // Compute ice threshold temperature
                    double q_I0 = ice_threshold_temperature(qi_val);

                    // If the ice mixing ratio is greater than the ice threshold temperature, compute the aggregation rate.
                    if (qi_val > q_I0) 
                    {
                        double agg_rate = c_agg_ * (qi_val - q_I0) * exp(0.025 * T_celsius);
                        dqi_dt[i][j][k] -= agg_rate;
                        dqs_dt[i][j][k] += agg_rate;
                    }
                }

                // If the ice mixing ratio is greater than 0, compute the autoconversion rate.
                if (qi_val > 0.0f) 
                {
                    double auto_rate = c_auto_ * qi_val * qi_val;  // Simplified
                    dqi_dt[i][j][k] -= auto_rate;
                    dqs_dt[i][j][k] += auto_rate;
                }

                // If the temperature is less than the freezing temperature, compute the sublimation rate.
                if (T < microphysics_constants::T0) 
                {
                    float qvsat_ice = thermodynamics::saturation_mixing_ratio_ice(T, P);
                    float RH_ice = (qv_val > 0.0f) ? qv_val / qvsat_ice : 0.0f;

                    // If the relative humidity is less than 1, compute the sublimation rate.
                    if (RH_ice < 1.0f) 
                    {
                        // If the ice mixing ratio is greater than 0, compute the sublimation rate.
                        if (qi_val > 0.0f) 
                        {
                            // Compute sublimation rate
                            float subl_rate = c_subl_ * (1.0f - RH_ice) * qi_val;
                            dqv_dt[i][j][k] += subl_rate;
                            dqi_dt[i][j][k] -= subl_rate;
                            dtheta_dt[i][j][k] -= microphysics_constants::L_s / microphysics_constants::cp * subl_rate / T;
                        }

                        // If the snow mixing ratio is greater than 0, compute the sublimation rate.
                        if (qs_val > 0.0f) 
                        {
                            // Compute sublimation rate
                            float subl_rate = c_subl_ * (1.0f - RH_ice) * qs_val;
                            dqv_dt[i][j][k] += subl_rate;
                            dqs_dt[i][j][k] -= subl_rate;
                            dtheta_dt[i][j][k] -= microphysics_constants::L_s / microphysics_constants::cp * subl_rate / T;
                        }
                    }
                }
            }
        }
    }
}


/*This function computes the melting processes for the Lin scheme.
Takes in the temperature, snow mixing ratio, graupel mixing ratio, hail mixing ratio, 
and the tendencies for the snow, graupel, hail, and rainwater mixing ratios
and computes the melting processes for the Lin scheme.*/
void LinScheme::compute_melting_processes(
    const std::vector<std::vector<std::vector<float>>>& temperature,
    const std::vector<std::vector<std::vector<float>>>& qs,
    const std::vector<std::vector<std::vector<float>>>& qg,
    const std::vector<std::vector<std::vector<float>>>& qh,
    std::vector<std::vector<std::vector<float>>>& dqs_dt,
    std::vector<std::vector<std::vector<float>>>& dqg_dt,
    std::vector<std::vector<std::vector<float>>>& dqh_dt,
    std::vector<std::vector<std::vector<float>>>& dqr_dt,
    std::vector<std::vector<std::vector<float>>>& dtheta_dt
) 
{
    int NR = temperature.size();
    int NTH = temperature[0].size();
    int NZ = temperature[0][0].size();

    // Compute melting processes by iterating over all grid points
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

                // If the temperature is greater than the freezing temperature, compute the melting rate.
                if (T > microphysics_constants::T0) 
                {
                    // If the snow mixing ratio is greater than 0, compute the melting rate.
                    if (qs_val > 0.0f) 
                    {
                        // Compute melting rate
                        float melt_rate = c_melt_ * qs_val;
                        dqs_dt[i][j][k] -= melt_rate;
                        dqr_dt[i][j][k] += melt_rate;
                        dtheta_dt[i][j][k] -= microphysics_constants::L_f / microphysics_constants::cp * melt_rate / T;
                    }

                    // If the graupel mixing ratio is greater than 0, compute the melting rate.
                    if (qg_val > 0.0f) 
                    {
                        // Compute melting rate
                        float melt_rate = c_melt_ * qg_val;
                        dqg_dt[i][j][k] -= melt_rate;
                        dqr_dt[i][j][k] += melt_rate;
                        dtheta_dt[i][j][k] -= microphysics_constants::L_f / microphysics_constants::cp * melt_rate / T;
                    }

                    // If the hail mixing ratio is greater than 0, compute the melting rate.
                    if (qh_val > 0.0f) 
                    {
                        // Compute melting rate
                        float melt_rate = c_melt_ * qh_val;
                        dqh_dt[i][j][k] -= melt_rate;
                        dqr_dt[i][j][k] += melt_rate;
                        dtheta_dt[i][j][k] -= microphysics_constants::L_f / microphysics_constants::cp * melt_rate / T;
                    }
                }
            }
        }
    }
}


/*This function computes the sedimentation for the Lin scheme.
Takes in the rainwater mixing ratio, snow mixing ratio, graupel mixing ratio, 
hail mixing ratio, and the tendencies for the rainwater, snow, graupel, 
and hail mixing ratios and computes the sedimentation for the Lin scheme.*/
void LinScheme::compute_sedimentation(
    const std::vector<std::vector<std::vector<float>>>& qr,
    const std::vector<std::vector<std::vector<float>>>& qs,
    const std::vector<std::vector<std::vector<float>>>& qg,
    const std::vector<std::vector<std::vector<float>>>& qh,
    std::vector<std::vector<std::vector<float>>>& dqr_dt,
    std::vector<std::vector<std::vector<float>>>& dqs_dt,
    std::vector<std::vector<std::vector<float>>>& dqg_dt,
    std::vector<std::vector<std::vector<float>>>& dqh_dt
) 
{
    // Get the number of rows, columns, and levels.
    int NR = qr.size();
    int NTH = qr[0].size();
    int NZ = qr[0][0].size();
    const double dz = 100.0;

    // Iterate over all grid points to compute the sedimentation.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal points
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical points
            for (int k = 0; k < NZ; ++k) 
            {
                // Rain sedimentation
                if (qr[i][j][k] > 0.0f && k > 0) 
                {
                    // Compute terminal velocity
                    float Vt = a_r_ * std::pow(std::max(qr[i][j][k], 1e-6f), b_r_);
                    float sed_flux = Vt * qr[i][j][k];
                    dqr_dt[i][j][k] -= sed_flux / dz;
                    if (k < NZ - 1) dqr_dt[i][j][k + 1] += sed_flux / dz;
                }

                // If the snow mixing ratio is greater than 0 and the vertical level is greater than 0, compute the sedimentation rate.
                if (qs[i][j][k] > 0.0f && k > 0) 
                {
                    // Compute terminal velocity
                    float Vt = a_s_ * std::pow(std::max(qs[i][j][k], 1e-6f), b_s_);
                    float sed_flux = Vt * qs[i][j][k];
                    dqs_dt[i][j][k] -= sed_flux / dz;
                    if (k < NZ - 1) dqs_dt[i][j][k + 1] += sed_flux / dz;
                }

                // If the graupel mixing ratio is greater than 0 and the vertical level is greater than 0, compute the sedimentation rate.
                if (qg[i][j][k] > 0.0f && k > 0) 
                {
                    // Compute terminal velocity
                    float Vt = a_g_ * std::pow(std::max(qg[i][j][k], 1e-6f), b_g_);
                    float sed_flux = Vt * qg[i][j][k];
                    dqg_dt[i][j][k] -= sed_flux / dz;
                    if (k < NZ - 1) dqg_dt[i][j][k + 1] += sed_flux / dz;
                }

                // If the hail mixing ratio is greater than 0 and the vertical level is greater than 0, compute the sedimentation rate.
                if (qh[i][j][k] > 0.0f && k > 0) 
                {
                    // Compute terminal velocity
                    float Vt = a_h_ * std::pow(std::max(qh[i][j][k], 1e-6f), b_h_);
                    float sed_flux = Vt * qh[i][j][k];
                    dqh_dt[i][j][k] -= sed_flux / dz;
                    if (k < NZ - 1) dqh_dt[i][j][k + 1] += sed_flux / dz;
                }
            }
        }
    }
}


/*This function computes the radar reflectivity for the Lin scheme.
Takes in the cloud water mixing ratio, rainwater mixing ratio, ice mixing ratio, snow mixing ratio, 
graupel mixing ratio, hail mixing ratio, and the radar reflectivity field
and computes the radar reflectivity for the Lin scheme.*/
void LinScheme::compute_radar_reflectivity(
    const std::vector<std::vector<std::vector<float>>>& qc,
    const std::vector<std::vector<std::vector<float>>>& qr,
    const std::vector<std::vector<std::vector<float>>>& qi,
    const std::vector<std::vector<std::vector<float>>>& qs,
    const std::vector<std::vector<std::vector<float>>>& qg,
    const std::vector<std::vector<std::vector<float>>>& qh,
    std::vector<std::vector<std::vector<float>>>& reflectivity_dbz
) 
{
    // Get grid size
    int NR = qc.size();
    if (NR == 0) return;
    int NTH = qc[0].size();
    int NZ = qc[0][0].size();

    reflectivity_dbz.assign(NR, std::vector<std::vector<float>>(NTH, std::vector<float>(NZ, 0.0f)));

    // Lin scheme reflectivity constants
    const float K_qc = 0.1e-3f;
    const float K_qr = 4.0e-3f;
    const float K_qi = 0.1e-3f;
    const float K_qs = 1.0e-3f;
    const float K_qg = 1.0e-2f;
    const float K_qh = 2.0e-2f;
    const float alpha = 1.5f;
    const float Z_min = 1e-10f;

    // Iterate over all grid points to compute the radar reflectivity.
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
                float qi_val = qi[i][j][k];
                float qs_val = qs[i][j][k];
                float qg_val = qg[i][j][k];
                float qh_val = qh[i][j][k];

                float Z_linear = K_qc * std::pow(std::max(qc_val, 0.0f), alpha) +
                                K_qr * std::pow(std::max(qr_val, 0.0f), alpha) +
                                K_qi * std::pow(std::max(qi_val, 0.0f), alpha) +
                                K_qs * std::pow(std::max(qs_val, 0.0f), alpha) +
                                K_qg * std::pow(std::max(qg_val, 0.0f), alpha) +
                                K_qh * std::pow(std::max(qh_val, 0.0f), alpha);

                Z_linear = std::max(Z_linear, Z_min);
                float Z_dBZ = 10.0f * std::log10(Z_linear);
                reflectivity_dbz[i][j][k] = Z_dBZ;
            }
        }
    }
}


/*This function computes the ice number concentration for the Lin scheme.
Takes in the ice mixing ratio and temperature and computes the ice number concentration for the Lin scheme.*/
double LinScheme::ice_number_concentration(double qi, double T_celsius) 
{
    // If the ice mixing ratio is less than 1e-6, compute the ice number concentration.
    if (qi < 1e-6) 
    {
        return 1e3;
    } 

    // If the ice mixing ratio is less than 1e-3, compute the ice number concentration.
    else if (qi < 1e-3) 
    {
        // Compute ice number concentration
        return 1e3 * std::pow(qi / 1e-6, 0.15);
    } 
    else 
    {
        // Compute ice number concentration
        return 1e3 * std::pow(10.0, 0.45);  // ~2.8e4
    }
}


/*This function computes the collection efficiency for the Lin scheme.
Takes in the temperature and computes the collection efficiency for the Lin scheme.*/
double LinScheme::collection_efficiency(double T_celsius) 
{
    // Temperature-dependent collection efficiency
    if (T_celsius < -5.0) 
    {
        // Collection efficiency
        return 1.0;
    } 
    else if (T_celsius <= 0.0) 
    {
        // Collection efficiency
        return 1.0 - (-5.0 - T_celsius) / 5.0;
    } 
    else 
    {
        // Collection efficiency
        return 0.0;
    }
}

// This is the ice threshold temperature function for the Lin scheme
double LinScheme::ice_threshold_temperature(double qi) 
{
    return 1e-3 * exp(0.1 * 0.0);  // Simplified - should depend on temperature
}
