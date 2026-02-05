#include "kessler.hpp"
#include <algorithm>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif

/* This is the constructor for the KesslerScheme class
// It initializes the parameters for the Kessler scheme
// This is simplified for the purpose of this project
// In a full implementation, these parameters would be configurable
// and the scheme would be more complex
// For example, the accretion rate would be a function of the cloud water and rain water mixing ratios
// and the evaporation rate would be a function of the rain water mixing ratio and the saturation mixing ratio
// and the homogeneous freezing rate would be a function of the cloud water mixing ratio
// and the riming rate would be a function of the cloud water and ice mixing ratios
// and the melting rate would be a function of the ice and graupel mixing ratios
// and the sublimation rate would be a function of the ice and hail mixing ratios*/


/*This is the constructor for the KesslerScheme class.
Takes in the parameters for the Kessler scheme and initializes the parameters.*/
KesslerScheme::KesslerScheme(
    double qc0, double c_auto, double c_accr, double c_evap,
    double c_freeze, double c_rime, double c_melt, double c_subl,
    double a_term, double b_term, double Vt_max,
    double a_hail, double b_hail, double a_grau, double b_grau,
    double Vt_max_hail, double Vt_max_grau
) : qc0_(qc0), c_auto_(c_auto), c_accr_(c_accr), c_evap_(c_evap),
    c_freeze_(c_freeze), c_rime_(c_rime), c_melt_(c_melt), c_subl_(c_subl),
    a_term_(a_term), b_term_(b_term), Vt_max_(Vt_max),
    a_hail_(a_hail), b_hail_(b_hail), a_grau_(a_grau), b_grau_(b_grau),
    Vt_max_hail_(Vt_max_hail), Vt_max_grau_(Vt_max_grau) {}


/*This function computes the tendencies for the Kessler scheme.
Takes in the pressure, potential temperature, vapor mixing ratio, cloud water mixing ratio, 
rainwater mixing ratio, ice mixing ratio, snow mixing ratio, graupel mixing ratio, and hail mixing ratio
and computes the tendencies for the Kessler scheme.*/
void KesslerScheme::compute_tendencies(
    const Field3D& p,
    const Field3D& theta,
    const Field3D& qv,
    const Field3D& qc,
    const Field3D& qr,
    const Field3D& qi,
    const Field3D& qs,
    const Field3D& qg,
    const Field3D& qh,
    double dt,
    Field3D& dtheta_dt,
    Field3D& dqv_dt,
    Field3D& dqc_dt,
    Field3D& dqr_dt,
    Field3D& dqi_dt,
    Field3D& dqs_dt,
    Field3D& dqg_dt,
    Field3D& dqh_dt
) 
{
    // Get grid size
    int NR = p.size_r();
    if (NR == 0) return;
    int NTH = p.size_th();
    if (NTH == 0) return;
    int NZ = p.size_z();

    // Convert theta to temperature for microphysics calculations
    Field3D temperature(NR, NTH, NZ);
    thermodynamics::convert_theta_to_temperature_field(theta, p, temperature);

    // Initialize tendency arrays to zero
    dtheta_dt.resize(NR, NTH, NZ, 0.0f);
    dqv_dt.resize(NR, NTH, NZ, 0.0f);
    dqc_dt.resize(NR, NTH, NZ, 0.0f);
    dqr_dt.resize(NR, NTH, NZ, 0.0f);
    dqi_dt.resize(NR, NTH, NZ, 0.0f);
    dqs_dt.resize(NR, NTH, NZ, 0.0f);
    dqg_dt.resize(NR, NTH, NZ, 0.0f);
    dqh_dt.resize(NR, NTH, NZ, 0.0f);

    // Compute microphysical processes
    compute_warm_rain_processes(temperature, qv, qc, qr, dqc_dt, dqr_dt, dqv_dt, dtheta_dt);
    compute_ice_processes(temperature, qv, qc, qr, qg, qh, dqc_dt, dqg_dt, dqh_dt, dqv_dt, dtheta_dt);
    compute_melting_processes(temperature, qg, qh, dqg_dt, dqh_dt, dqr_dt, dtheta_dt);
    compute_sedimentation(qr, qg, qh, dqr_dt, dqg_dt, dqh_dt);

    // Iterate over all rows.
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal points.
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical points.
            for (int k = 0; k < NZ; ++k) 
            {
                // Convert temperature tendency to potential temperature tendency
                dtheta_dt[i][j][k] = static_cast<float>(
                    thermodynamics::temperature_tendency_to_theta(
                        static_cast<float>(dtheta_dt[i][j][k]), static_cast<float>(theta[i][j][k]), static_cast<float>(p[i][j][k])
                    )
                );
            }
        }
    }
}

/*This function computes the warm rain processes for the Kessler scheme.
Takes in the temperature, vapor mixing ratio, cloud water mixing ratio, 
rainwater mixing ratio, and the tendencies for the cloud water, rainwater, and vapor mixing ratios
and computes the warm rain processes for the Kessler scheme.*/
void KesslerScheme::compute_warm_rain_processes(
    const Field3D& temperature,
    const Field3D& qv,
    const Field3D& qc,
    const Field3D& qr,
    Field3D& dqc_dt,
    Field3D& dqr_dt,
    Field3D& dqv_dt,
    Field3D& dtheta_dt
) 
{
    // Get grid size
    int NR = temperature.size_r();
    int NTH = temperature.size_th();
    int NZ = temperature.size_z();

    // Compute warm rain processes by iterating over all grid points
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal points
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical points
            for (int k = 0; k < NZ; ++k) 
            {
                // Get hydrometeor values
                float qc_val = static_cast<float>(qc[i][j][k]);
                float qr_val = static_cast<float>(qr[i][j][k]);
                float qv_val = static_cast<float>(qv[i][j][k]);
                float T = static_cast<float>(temperature[i][j][k]);

                // If the cloud water mixing ratio is greater than the critical cloud water mixing ratio, compute the autoconversion rate.
                // Autoconversion: qc → qr when qc > qc0
                if (qc_val > qc0_) 
                {
                    // Compute autoconversion rate
                    float auto_rate = c_auto_ * (qc_val - qc0_);
                    dqc_dt[i][j][k] -= auto_rate;
                    dqr_dt[i][j][k] += auto_rate;
                    dtheta_dt[i][j][k] += microphysics_constants::L_v / microphysics_constants::cp * auto_rate / T;
                }

                // if the cloud water mixing ratio is greater than 0 and the rainwater mixing ratio is greater than 0, compute the accretion rate.
                if (qc_val > 0.0f && qr_val > 0.0f) 
                {
                    // Compute accretion rate
                    float accr_rate = c_accr_ * qc_val * qr_val;
                    dqc_dt[i][j][k] -= accr_rate;
                    dqr_dt[i][j][k] += accr_rate;
                    dtheta_dt[i][j][k] += microphysics_constants::L_v / microphysics_constants::cp * accr_rate / T;
                }

                // Evaporation: qr → qv when subsaturated
                if (qr_val > 0.0f) 
                {
                    float qvsat = 0.001f;  // Simplified - should use proper saturation
                    float RH = (qv_val > 0.0f) ? qv_val / qvsat : 0.0f;

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

/*This function computes the ice processes for the Kessler scheme.
Takes in the temperature, vapor mixing ratio, cloud water mixing ratio, 
rainwater mixing ratio, ice mixing ratio, snow mixing ratio, graupel mixing ratio, and hail mixing ratio
and computes the ice processes for the Kessler scheme.*/
void KesslerScheme::compute_ice_processes(
    const Field3D& temperature,
    const Field3D& qv,
    const Field3D& qc,
    const Field3D& qr,
    const Field3D& qg,
    const Field3D& qh,
    Field3D& dqc_dt,
    Field3D& dqg_dt,
    Field3D& dqh_dt,
    Field3D& dqv_dt,
    Field3D& dtheta_dt
) 
{
    // Get grid size
    int NR = temperature.size_r();
    int NTH = temperature.size_th();
    int NZ = temperature.size_z();

    // Iterate over all grid points.
    #pragma omp parallel for collapse(2)
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
                float qg_val = qg[i][j][k];
                float qh_val = qh[i][j][k];
                float qv_val = qv[i][j][k];
                float T = temperature[i][j][k];

                bool is_cold = (T < microphysics_constants::T0);

                // If the temperature is less than the freezing temperature, compute the ice processes.
                if (is_cold) 
                {
                    // If the cloud water mixing ratio is greater than 0, compute the homogeneous freezing rate.
                    // Homogeneous freezing: qc → qg
                    if (qc_val > 0.0f) 
                    {
                        // Compute homogeneous freezing rate
                        float freeze_rate = c_freeze_ * qc_val;
                        dqc_dt[i][j][k] -= freeze_rate;
                        dqg_dt[i][j][k] += freeze_rate;
                        dtheta_dt[i][j][k] += microphysics_constants::L_f / microphysics_constants::cp * freeze_rate / T;
                    }

                    // If the cloud water mixing ratio is greater than 0 and the graupel mixing ratio is greater than 0, compute the riming rate.
                    // Riming: ice collects qc
                    if (qc_val > 0.0f && qg_val > 0.0f) 
                    {
                        float rime_rate = c_rime_ * qc_val * qg_val;
                        dqc_dt[i][j][k] -= rime_rate;
                        dqg_dt[i][j][k] += rime_rate;
                        dtheta_dt[i][j][k] += microphysics_constants::L_f / microphysics_constants::cp * rime_rate / T;
                    }

                    // If the cloud water mixing ratio is greater than 0 and the hail mixing ratio is greater than 0, compute the riming rate.
                    // Riming: ice collects qh
                    if (qc_val > 0.0f && qh_val > 0.0f) 
                    {
                        // Compute riming rate
                        float rime_rate = c_rime_ * qc_val * qh_val;
                        dqc_dt[i][j][k] -= rime_rate;
                        dqh_dt[i][j][k] += rime_rate;
                        dtheta_dt[i][j][k] += microphysics_constants::L_f / microphysics_constants::cp * rime_rate / T;
                    }

                    // If the graupel mixing ratio is greater than the critical graupel mixing ratio, compute the autoconversion rate.
                    // Conversion: graupel → hail
                    if (qg_val > qc0_) 
                    {
                        float conv_rate = c_auto_ * (qg_val - qc0_);
                        dqg_dt[i][j][k] -= conv_rate;
                        dqh_dt[i][j][k] += conv_rate;
                    }

                    // Sublimation: ice → qv when subsaturated
                    float qvsat_ice = 0.001f;  // Simplified
                    float RH_ice = (qv_val > 0.0f) ? qv_val / qvsat_ice : 0.0f;

                    // If the relative humidity is less than 1, compute the sublimation rate.
                    if (RH_ice < 1.0f) 
                    {
                        // If the graupel mixing ratio is greater than 0, compute the sublimation rate.
                        // Sublimation: graupel → qv when subsaturated
                        if (qg_val > 0.0f) 
                        {
                            // Compute sublimation rate
                            float subl_rate = c_subl_ * (1.0f - RH_ice) * qg_val;
                            dqv_dt[i][j][k] += subl_rate;
                            dqg_dt[i][j][k] -= subl_rate;
                            dtheta_dt[i][j][k] -= microphysics_constants::L_s / microphysics_constants::cp * subl_rate / T;
                        }

                        // If the hail mixing ratio is greater than 0, compute the sublimation rate.
                        // Sublimation: hail → qv when subsaturated
                        if (qh_val > 0.0f)
                        {
                            // Compute sublimation rate
                            float subl_rate = c_subl_ * (1.0f - RH_ice) * qh_val;
                            dqv_dt[i][j][k] += subl_rate;
                            dqh_dt[i][j][k] -= subl_rate;
                            dtheta_dt[i][j][k] -= microphysics_constants::L_s / microphysics_constants::cp * subl_rate / T;
                        }
                    }
                }
            }
        }
    }
}

/*This function computes the melting processes for the Kessler scheme.
Takes in the temperature, graupel mixing ratio, hail mixing ratio, and 
the tendencies for the graupel, hail, and rainwater mixing ratios
and computes the melting processes for the Kessler scheme.*/
void KesslerScheme::compute_melting_processes(
    const Field3D& temperature,
    const Field3D& qg,
    const Field3D& qh,
    Field3D& dqg_dt,
    Field3D& dqh_dt,
    Field3D& dqr_dt,
    Field3D& dtheta_dt
) 
{
    // Get grid size
    int NR = temperature.size_r();
    int NTH = temperature.size_th();
    int NZ = temperature.size_z();

    // Iterate over all grid points.
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal points
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical points
            for (int k = 0; k < NZ; ++k) 
            {
                // Get hydrometeor values
                float qg_val = static_cast<float>(qg[i][j][k]);
                float qh_val = static_cast<float>(qh[i][j][k]);
                float T = static_cast<float>(temperature[i][j][k]);

                // If the temperature is greater than the freezing temperature, compute the melting processes.
                if (T > microphysics_constants::T0) 
                {
                    // If the hail mixing ratio is greater than 0, compute the melting rate.
                    // Melting: ice → water
                    if (qh_val > 0.0f) 
                    {
                        // Compute melting rate
                        float melt_rate = c_melt_ * qh_val;
                        dqh_dt[i][j][k] -= melt_rate;
                        dqr_dt[i][j][k] += melt_rate;
                        dtheta_dt[i][j][k] -= microphysics_constants::L_f / microphysics_constants::cp * melt_rate / T;
                    }

                    // If the graupel mixing ratio is greater than 0, compute the melting rate.
                    // Melting: graupel → water
                    if (qg_val > 0.0f) 
                    {
                        float melt_rate = c_melt_ * qg_val;
                        dqg_dt[i][j][k] -= melt_rate;
                        dqr_dt[i][j][k] += melt_rate;
                        dtheta_dt[i][j][k] -= microphysics_constants::L_f / microphysics_constants::cp * melt_rate / T;
                    }
                }
            }
        }
    }
}

/*This function computes the sedimentation for the Kessler scheme.
Takes in the rainwater mixing ratio, graupel mixing ratio, hail mixing ratio, 
and the tendencies for the rainwater, graupel, and hail mixing ratios
and computes the sedimentation for the Kessler scheme.*/
void KesslerScheme::compute_sedimentation(
    const Field3D& qr,
    const Field3D& qg,
    const Field3D& qh,
    Field3D& dqr_dt,
    Field3D& dqg_dt,
    Field3D& dqh_dt
) 
{
    

    int NR = qr.size_r();
    int NTH = qr.size_th();
    int NZ = qr.size_z();

    // Assume dz = 100.0 m (should be passed as parameter)
    const double dz = 100.0;

    // Iterate over all grid points.
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal points
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical points
            for (int k = 0; k < NZ; ++k) 
            {
                // Rain sedimentation
                float qr_val = static_cast<float>(qr[i][j][k]);
                if (qr_val > 0.0f && k > 0) 
                {
                    // Compute terminal velocity
                    float Vt_val = a_term_ * std::pow(std::max(qr_val, 1e-6f), b_term_);
                    float Vt = std::min(static_cast<float>(Vt_max_), Vt_val);
                    float sed_flux = Vt * qr_val;
                    dqr_dt[i][j][k] -= sed_flux / dz;
                    if (k < NZ - 1) dqr_dt[i][j][k + 1] += sed_flux / dz;
                }

                // If the graupel mixing ratio is greater than 0 and the vertical level is greater than 0, compute the sedimentation rate.
                float qg_val = static_cast<float>(qg[i][j][k]);
                if (qg_val > 0.0f && k > 0) 
                {
                    // Compute terminal velocity
                    float Vt_val = a_grau_ * std::pow(std::max(qg_val, 1e-6f), b_grau_);
                    float Vt = std::min(static_cast<float>(Vt_max_grau_), Vt_val);
                    float sed_flux = Vt * qg_val;
                    dqg_dt[i][j][k] -= sed_flux / dz;
                    if (k < NZ - 1) dqg_dt[i][j][k + 1] += sed_flux / dz;
                }

                // If the hail mixing ratio is greater than 0 and the vertical level is greater than 0, compute the sedimentation rate.
                // Hail sedimentation
                float qh_val = static_cast<float>(qh[i][j][k]);
                if (qh_val > 0.0f && k > 0) 
                {
                    // Compute terminal velocity
                    float Vt_val = a_hail_ * std::pow(std::max(qh_val, 1e-6f), b_hail_);
                    float Vt = std::min(static_cast<float>(Vt_max_hail_), Vt_val);
                    float sed_flux = Vt * qh_val;
                    dqh_dt[i][j][k] -= sed_flux / dz;

                    // If the vertical level is less than the number of vertical levels minus 1, compute the sedimentation rate.
                    if (k < NZ - 1) dqh_dt[i][j][k + 1] += sed_flux / dz;
                }
            }
        }
    }
}

/*This function computes the radar reflectivity for the Kessler scheme.
Takes in the cloud water mixing ratio, rainwater mixing ratio, ice mixing ratio, 
snow mixing ratio, graupel mixing ratio, hail mixing ratio, and the radar 
reflectivity field and computes the radar reflectivity for the Kessler scheme.*/
void KesslerScheme::compute_radar_reflectivity(
    const Field3D& qc,
    const Field3D& qr,
    const Field3D& qi,
    const Field3D& qs,
    const Field3D& qg,
    const Field3D& qh,
    Field3D& reflectivity_dbz
) 
{
    // Get grid size
    int NR = qc.size_r();
    if (NR == 0) return;
    int NTH = qc.size_th();
    if (NTH == 0) return;
    int NZ = qc.size_z();

    reflectivity_dbz.resize(NR, NTH, NZ, 0.0f);

    // Radar reflectivity constants (S-band)
    const float K_qc = 0.1e-3f;    // Cloud water
    const float K_qr = 4.0e-3f;    // Rain
    const float K_qi = 0.1e-3f;    // Cloud ice
    const float K_qs = 1.0e-3f;    // Snow
    const float K_qg = 1.0e-2f;    // Graupel
    const float K_qh = 2.0e-2f;    // Hail
    const float alpha = 1.5f;      // Exponent
    const float Z_min = 1e-10f;    // Minimum reflectivity

    //Iterate over all grid points to compute the radar reflectivity.
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal points.
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical points.
            for (int k = 0; k < NZ; ++k) 
            {
                // Get hydrometeor values
                float qc_val = static_cast<float>(qc[i][j][k]);
                float qr_val = static_cast<float>(qr[i][j][k]);
                float qi_val = static_cast<float>(qi[i][j][k]);
                float qs_val = static_cast<float>(qs[i][j][k]);
                float qg_val = static_cast<float>(qg[i][j][k]);
                float qh_val = static_cast<float>(qh[i][j][k]);

                // Compute linear radar reflectivity
                float Z_linear = K_qc * std::pow(std::max(qc_val, 0.0f), alpha) +
                                K_qr * std::pow(std::max(qr_val, 0.0f), alpha) +
                                K_qi * std::pow(std::max(qi_val, 0.0f), alpha) +
                                K_qs * std::pow(std::max(qs_val, 0.0f), alpha) +
                                K_qg * std::pow(std::max(qg_val, 0.0f), alpha) +
                                K_qh * std::pow(std::max(qh_val, 0.0f), alpha);

                // Ensure reflectivity is not negative
                Z_linear = std::max(Z_linear, Z_min);
                float Z_dBZ = 10.0f * std::log10(Z_linear);
                reflectivity_dbz[i][j][k] = Z_dBZ;
            }
        }
    }
}
