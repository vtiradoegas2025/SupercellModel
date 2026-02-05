#include "thompson.hpp"
#include <algorithm>
#include <cmath>
#ifdef _OPENMP
#include <omp.h>
#endif


/*This constructor initializes the Thompson scheme with default parameters.
Takes in the autoconversion threshold, autoconversion rate, cloud number concentration, 
ice number concentration, a parameters for the rain, b parameters for the rain, a parameters for the snow, 
b parameters for the snow, a parameters for the graupel, b parameters for the graupel, a parameters for the hail, 
and b parameters for the hail and initializes the Thompson scheme with the default parameters.*/
ThompsonScheme::ThompsonScheme(
    double qc0, double c_auto, double ccn_conc, double in_conc,
    double a_r, double b_r, double a_s, double b_s, double a_g, double b_g, double a_h, double b_h
) : qc0_(qc0), c_auto_(c_auto), ccn_conc_(ccn_conc), in_conc_(in_conc),
    a_r_(a_r), b_r_(b_r), a_s_(a_s), b_s_(b_s),
    a_g_(a_g), b_g_(b_g), a_h_(a_h), b_h_(b_h) {}

/*This function computes the tendencies for the Thompson scheme.
Takes in the pressure, potential temperature, vapor mixing ratio, cloud water mixing ratio, 
rainwater mixing ratio, ice mixing ratio, snow mixing ratio, graupel mixing ratio, hail mixing ratio,
and the time step and computes the tendencies for the Thompson scheme.*/
void ThompsonScheme::compute_tendencies(
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

    int NR = p.size_r();
    if (NR == 0) return;
    int NTH = p.size_th();
    int NZ = p.size_z();

    // If the ice number concentration is not already initialized, initialize it.
    if (Ni_.empty()) 
    {
        Ni_.resize(NR, NTH, NZ, in_conc_ * 1e-6f);  // Convert m⁻³ to cm⁻³
    }

    // Convert theta to temperature
    Field3D temperature(NR, NTH, NZ);
    thermodynamics::convert_theta_to_temperature_field(theta, p, temperature);

    // Initialize tendency arrays
    dtheta_dt.resize(NR, NTH, NZ, 0.0f);
    dqv_dt.resize(NR, NTH, NZ, 0.0f);
    dqc_dt.resize(NR, NTH, NZ, 0.0f);
    dqr_dt.resize(NR, NTH, NZ, 0.0f);
    dqi_dt.resize(NR, NTH, NZ, 0.0f);
    dqs_dt.resize(NR, NTH, NZ, 0.0f);
    dqg_dt.resize(NR, NTH, NZ, 0.0f);
    dqh_dt.resize(NR, NTH, NZ, 0.0f);
    dNi_dt_.resize(NR, NTH, NZ, 0.0f);

    // Saturation adjustment
    Field3D qv_temp = qv;
    Field3D qc_temp = qc;
    saturation_adjustment(temperature, p, qv_temp, qc_temp);

    // Compute microphysical processes
    compute_warm_rain_processes(temperature, p, qv_temp, qc_temp, qr, dqc_dt, dqr_dt, dqv_dt, dtheta_dt);
    compute_ice_processes(temperature, p, qv_temp, qc_temp, qi, qs, qg, qh,
                         dqc_dt, dqi_dt, dqs_dt, dqg_dt, dqh_dt, dqv_dt, dtheta_dt, dNi_dt_);
    compute_melting_processes(temperature, qs, qg, qh, dqs_dt, dqg_dt, dqh_dt, dqr_dt, dtheta_dt);
    compute_sedimentation(qr, qs, qg, qh, dqr_dt, dqs_dt, dqg_dt, dqh_dt);

    // Iterate over all grid points and update Ni
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal angles
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical levels
            for (int k = 0; k < NZ; ++k) 
            {
                // Update Ni
                Ni_[i][j][k] += dNi_dt_[i][j][k] * dt;
                Ni_[i][j][k] = std::max(static_cast<float>(Ni_[i][j][k]), 1.0f);  // Minimum concentration
            }
        }
    }

    // Iterate over all grid points and update temperature tendencies to potential temperature tendencies
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal angles
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical levels
            for (int k = 0; k < NZ; ++k) 
            {
                // Update temperature tendencies to potential temperature tendencies
                dtheta_dt[i][j][k] = static_cast<float>(
                    thermodynamics::temperature_tendency_to_theta(
                        static_cast<float>(dtheta_dt[i][j][k]), static_cast<float>(theta[i][j][k]), static_cast<float>(p[i][j][k])
                    )
                );
            }
        }
    }
}

/*This function performs the saturation adjustment for the Thompson scheme.
Takes in the temperature, pressure, vapor mixing ratio, and cloud water mixing ratio
and performs the saturation adjustment for the Thompson scheme.*/
void ThompsonScheme::saturation_adjustment(
    const Field3D& temperature,
    const Field3D& p,
    Field3D& qv,
    Field3D& qc
) 
{
    int NR = temperature.size_r();
    int NTH = temperature.size_th();
    int NZ = temperature.size_z();

    // Iterate over all grid points and perform saturation adjustment
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal angles
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical levels
            for (int k = 0; k < NZ; ++k) 
            {
                // Get current temperature and pressure
                float T = static_cast<float>(temperature[i][j][k]);
                float P = static_cast<float>(p[i][j][k]);
                float qvsat = thermodynamics::saturation_mixing_ratio_water(T, P);

                // If the vapor mixing ratio is greater than the saturation mixing ratio, compute the excess
                float qv_val = static_cast<float>(qv[i][j][k]);
                if (qv_val > qvsat) 
                {
                    float excess = qv_val - qvsat;
                    qc[i][j][k] += excess;
                    qv[i][j][k] = qvsat;
                }
            }
        }
    }
}

/*This function computes the warm rain processes for the Thompson scheme.
Takes in the temperature, pressure, vapor mixing ratio, cloud water mixing ratio, 
rainwater mixing ratio, and the tendencies for the cloud water, rainwater, and vapor mixing ratios
and computes the warm rain processes for the Thompson scheme.*/
void ThompsonScheme::compute_warm_rain_processes(
    const Field3D& temperature,
    const Field3D& p,
    const Field3D& qv,
    const Field3D& qc,
    const Field3D& qr,
    Field3D& dqc_dt,
    Field3D& dqr_dt,
    Field3D& dqv_dt,
    Field3D& dtheta_dt
) 
{
    // Iterate over all grid points and compute warm rain processes
    int NR = temperature.size_r();
    int NTH = temperature.size_th();
    int NZ = temperature.size_z();

    // Iterate over all grid points and compute warm rain processes
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal angles
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical levels
            for (int k = 0; k < NZ; ++k) 
            {
                float qc_val = static_cast<float>(qc[i][j][k]);
                float qr_val = static_cast<float>(qr[i][j][k]);
                float qv_val = static_cast<float>(qv[i][j][k]);
                float T = static_cast<float>(temperature[i][j][k]);
                float P = static_cast<float>(p[i][j][k]);

                // If the cloud water mixing ratio is greater than the autoconversion threshold, compute the autoconversion rate
                if (qc_val > qc0_) 
                {
                    // Compute the autoconversion rate
                    float auto_rate = c_auto_ * (qc_val - qc0_);

                    // Update the tendencies
                    dqc_dt[i][j][k] -= auto_rate;
                    dqr_dt[i][j][k] += auto_rate;
                    dtheta_dt[i][j][k] += microphysics_constants::L_v / microphysics_constants::cp * auto_rate / T;
                }

                // If the cloud water mixing ratio and rain water mixing ratio are greater than zero, compute the accretion rate
                // Accretion: qr collects qc
                if (qc_val > 0.0f && qr_val > 0.0f) 
                {
                    float accr_rate = 2.2 * qc_val * qr_val;  // Simplified
                    dqc_dt[i][j][k] -= accr_rate;
                    dqr_dt[i][j][k] += accr_rate;
                    dtheta_dt[i][j][k] += microphysics_constants::L_v / microphysics_constants::cp * accr_rate / T;
                }

                // If the rain water mixing ratio is greater than zero, compute the evaporation rate
                // Evaporation: qr → qv
                if (qr_val > 0.0f) 
                {
                    // Compute the saturation mixing ratio
                    float qvsat = thermodynamics::saturation_mixing_ratio_water(T, P);
                    float RH = (qv_val > 0.0f) ? qv_val / qvsat : 0.0f;

                    // If the relative humidity is less than 100%, compute the evaporation rate
                    if (RH < 1.0f) 
                    {
                        float evap_rate = 1.0e-3 * (1.0f - RH) * qr_val;  // Simplified
                        dqv_dt[i][j][k] += evap_rate;
                        dqr_dt[i][j][k] -= evap_rate;
                        dtheta_dt[i][j][k] -= microphysics_constants::L_v / microphysics_constants::cp * evap_rate / T;
                    }
                }
            }
        }
    }
}

/*This function computes the ice processes for the Thompson scheme.
Takes in the temperature, pressure, vapor mixing ratio, cloud water mixing ratio, 
ice mixing ratio, snow mixing ratio, graupel mixing ratio, hail mixing ratio,
and the tendencies for the cloud water, ice, snow, graupel, hail, and vapor mixing ratios
and computes the ice processes for the Thompson scheme.*/
void ThompsonScheme::compute_ice_processes(
    const Field3D& temperature,
    const Field3D& p,
    const Field3D& qv,
    const Field3D& qc,
    const Field3D& qi,
    const Field3D& qs,
    const Field3D& qg,
    const Field3D& qh,
    Field3D& dqc_dt,
    Field3D& dqi_dt,
    Field3D& dqs_dt,
    Field3D& dqg_dt,
    Field3D& dqh_dt,
    Field3D& dqv_dt,
    Field3D& dtheta_dt,
    Field3D& dNi_dt
) 
{
    // Iterate over all grid points and compute ice processes
    int NR = temperature.size_r();
    int NTH = temperature.size_th();
    int NZ = temperature.size_z();

    // Iterate over all grid points and compute ice processes
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal angles
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical levels
            for (int k = 0; k < NZ; ++k) 
            {
                float qc_val = static_cast<float>(qc[i][j][k]);
                float qi_val = static_cast<float>(qi[i][j][k]);
                float qs_val = static_cast<float>(qs[i][j][k]);
                float qg_val = static_cast<float>(qg[i][j][k]);
                float qv_val = static_cast<float>(qv[i][j][k]);
                float T = static_cast<float>(temperature[i][j][k]);
                float P = static_cast<float>(p[i][j][k]);
                float T_celsius = T - microphysics_constants::T0;
                float Ni_val = Ni_[i][j][k];

                // Homogeneous freezing of cloud droplets (Thompson freezing probability)
                if (qc_val > 0.0f && T < microphysics_constants::T0) 
                {
                    // Simplified droplet volume (assuming 10 µm radius)
                    double qc_vol = qc_val / (4.0/3.0 * M_PI * 1e-15);  // m⁻³
                    double P_freeze = freezing_probability(T, qc_vol);

                    // If the freezing probability is greater than zero, compute the freezing rate
                    if (P_freeze > 0.0) 
                    {
                        double freeze_rate = P_freeze * qc_val;
                        dqc_dt[i][j][k] -= freeze_rate;
                        dqi_dt[i][j][k] += freeze_rate;
                        dNi_dt[i][j][k] += P_freeze * ccn_conc_ * 1e-6;  // Add to ice number
                        dtheta_dt[i][j][k] += microphysics_constants::L_f / microphysics_constants::cp * freeze_rate / T;
                    }
                }

                // If the ice mixing ratio is greater than zero and the temperature is less than the freezing point, compute the ice growth/deposition rate
                // Ice growth/deposition
                if (qi_val > 0.0f && T < microphysics_constants::T0) 
                {
                    float qvsat_ice = thermodynamics::saturation_mixing_ratio_ice(T, P);

                    // If the vapor mixing ratio is greater than the saturation mixing ratio, compute the deposition rate
                    if (qv_val > qvsat_ice) 
                    {
                        float dep_rate = 1.0e-3 * (qv_val - qvsat_ice);  // Simplified
                        dqv_dt[i][j][k] -= dep_rate;
                        dqi_dt[i][j][k] += dep_rate;
                        dtheta_dt[i][j][k] += microphysics_constants::L_s / microphysics_constants::cp * dep_rate / T;
                    }
                }

                // If the ice mixing ratio is greater than zero and the ice number concentration is greater than zero, compute the large ice conversion rate
                // 200 µm rule: convert large ice to snow
                if (qi_val > 0.0f && Ni_val > 0.0f) 
                {
                    // Compute the large ice conversion rate
                    double D_large = ice_size_threshold(Ni_val);

                    // If the ice size is greater than 200 µm, compute the large ice conversion rate
                    if (D_large > 200e-6) 
                    {  
                        double f_large = 0.1;  // Simplified fraction
                        double conv_rate = f_large * qi_val;
                        dqi_dt[i][j][k] -= conv_rate;
                        dqs_dt[i][j][k] += conv_rate;
                    }
                }

                // If the cloud water mixing ratio and ice mixing ratio are greater than zero, compute the riming rate
                // Riming: graupel formation
                if (qc_val > 0.0f && qi_val > 0.0f) 
                {
                    double E_rime = collection_efficiency(T_celsius, qi_val, qc_val);
                    double rime_rate = E_rime * qc_val * qi_val;  // Simplified
                    dqc_dt[i][j][k] -= rime_rate;
                    dqg_dt[i][j][k] += rime_rate;
                    dtheta_dt[i][j][k] += microphysics_constants::L_f / microphysics_constants::cp * rime_rate / T;
                }

                // If the temperature is less than the freezing point, compute the sublimation rate
                // Sublimation of ice/snow
                if (T < microphysics_constants::T0) 
                {
                    float qvsat_ice = thermodynamics::saturation_mixing_ratio_ice(T, P);
                    float RH_ice = (qv_val > 0.0f) ? qv_val / qvsat_ice : 0.0f;

                    // If the relative humidity is less than 100%, compute the sublimation rate
                    if (RH_ice < 1.0f) 
                    {
                        // If the ice mixing ratio is greater than zero, compute the sublimation rate
                        if (qi_val > 0.0f) 
                        {
                            float subl_rate = 1.0e-3 * (1.0f - RH_ice) * qi_val;
                            dqv_dt[i][j][k] += subl_rate;
                            dqi_dt[i][j][k] -= subl_rate;
                            dtheta_dt[i][j][k] -= microphysics_constants::L_s / microphysics_constants::cp * subl_rate / T;
                        }

                        // If the snow mixing ratio is greater than zero, compute the sublimation rate
                        if (qs_val > 0.0f) 
                        {
                            float subl_rate = 1.0e-3 * (1.0f - RH_ice) * qs_val;
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

/*This function computes the melting processes for the Thompson scheme.
Takes in the temperature, snow mixing ratio, graupel mixing ratio, hail mixing ratio, 
and the tendencies for the snow, graupel, hail, and rainwater mixing ratios
and computes the melting processes for the Thompson scheme.*/
void ThompsonScheme::compute_melting_processes(
    const Field3D& temperature,
    const Field3D& qs,
    const Field3D& qg,
    const Field3D& qh,
    Field3D& dqs_dt,
    Field3D& dqg_dt,
    Field3D& dqh_dt,
    Field3D& dqr_dt,
    Field3D& dtheta_dt
) 
{
    int NR = temperature.size_r();
    int NTH = temperature.size_th();
    int NZ = temperature.size_z();

    // Iterate over all grid points and compute melting processes
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal angles
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical levels and compute melting processes
            for (int k = 0; k < NZ; ++k) 
            {
                float qs_val = static_cast<float>(qs[i][j][k]);
                float qg_val = static_cast<float>(qg[i][j][k]);
                float qh_val = static_cast<float>(qh[i][j][k]);
                float T = static_cast<float>(temperature[i][j][k]);

                // If the temperature is greater than the freezing point, compute the melting rate
                if (T > microphysics_constants::T0) 
                {
                    // If the snow mixing ratio is greater than zero, compute the melting rate
                    if (qs_val > 0.0f) 
                    {
                        float melt_rate = 1.0e-3 * qs_val;
                        dqs_dt[i][j][k] -= melt_rate;
                        dqr_dt[i][j][k] += melt_rate;
                        dtheta_dt[i][j][k] -= microphysics_constants::L_f / microphysics_constants::cp * melt_rate / T;
                    }

                    // If the graupel mixing ratio is greater than zero, compute the melting rate
                    if (qg_val > 0.0f) 
                    {
                        float melt_rate = 1.0e-3 * qg_val;
                        dqg_dt[i][j][k] -= melt_rate;
                        dqr_dt[i][j][k] += melt_rate;
                        dtheta_dt[i][j][k] -= microphysics_constants::L_f / microphysics_constants::cp * melt_rate / T;
                    }

                    // If the hail mixing ratio is greater than zero, compute the melting rate
                    if (qh_val > 0.0f) 
                    {
                        float melt_rate = 1.0e-3 * qh_val;
                        dqh_dt[i][j][k] -= melt_rate;
                        dqr_dt[i][j][k] += melt_rate;
                        dtheta_dt[i][j][k] -= microphysics_constants::L_f / microphysics_constants::cp * melt_rate / T;
                    }
                }
            }
        }
    }
}

/*This function computes the sedimentation for the Thompson scheme.
Takes in the rainwater mixing ratio, snow mixing ratio, graupel mixing ratio, hail mixing ratio, 
and the tendencies for the rainwater, snow, graupel, and hail mixing ratios
and computes the sedimentation for the Thompson scheme.*/
void ThompsonScheme::compute_sedimentation(
    const Field3D& qr,
    const Field3D& qs,
    const Field3D& qg,
    const Field3D& qh,
    Field3D& dqr_dt,
    Field3D& dqs_dt,
    Field3D& dqg_dt,
    Field3D& dqh_dt
) 
{
    // Iterate over all grid points and compute sedimentation
    int NR = qr.size_r();
    int NTH = qr.size_th();
    int NZ = qr.size_z();
    const double dz = 100.0;

    // Iterate over all grid points and compute sedimentation
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal angles
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical levels
            for (int k = 0; k < NZ; ++k) 
            {
                // Rain sedimentation
                float qr_val = static_cast<float>(qr[i][j][k]);
                if (qr_val > 0.0f && k > 0)
                {
                    float Vt = a_r_ * std::pow(std::max(qr_val, 1e-6f), b_r_);
                    float sed_flux = Vt * qr_val;
                    dqr_dt[i][j][k] -= sed_flux / dz;
                    if (k < NZ - 1) dqr_dt[i][j][k + 1] += sed_flux / dz;
                }

                // If the snow mixing ratio is greater than zero and the vertical level is greater than zero, compute the sedimentation rate
                // Snow sedimentation
                float qs_val = static_cast<float>(qs[i][j][k]);
                if (qs_val > 0.0f && k > 0) 
                {
                    float Vt = a_s_ * std::pow(std::max(qs_val, 1e-6f), b_s_);
                    float sed_flux = Vt * qs_val;
                    dqs_dt[i][j][k] -= sed_flux / dz;
                    if (k < NZ - 1) dqs_dt[i][j][k + 1] += sed_flux / dz;
                }

                // If the graupel mixing ratio is greater than zero and the vertical level is greater than zero, compute the sedimentation rate
                // Graupel sedimentation
                float qg_val = static_cast<float>(qg[i][j][k]);
                if (qg_val > 0.0f && k > 0) 
                {
                    float Vt = a_g_ * std::pow(std::max(qg_val, 1e-6f), b_g_);
                    float sed_flux = Vt * qg_val;
                    dqg_dt[i][j][k] -= sed_flux / dz;
                    if (k < NZ - 1) dqg_dt[i][j][k + 1] += sed_flux / dz;
                }

                // If the hail mixing ratio is greater than zero and the vertical level is greater than zero, compute the sedimentation rate
                // Hail sedimentation
                float qh_val = static_cast<float>(qh[i][j][k]);
                if (qh_val > 0.0f && k > 0) 
                {
                    float Vt = a_h_ * std::pow(std::max(qh_val, 1e-6f), b_h_);
                    float sed_flux = Vt * qh_val;
                    dqh_dt[i][j][k] -= sed_flux / dz;
                    if (k < NZ - 1) dqh_dt[i][j][k + 1] += sed_flux / dz;
                }
            }
        }
    }
}

/*This function computes the radar reflectivity for the Thompson scheme.
Takes in the cloud water mixing ratio, rainwater mixing ratio, ice mixing ratio, snow mixing ratio, 
graupel mixing ratio, hail mixing ratio, and the radar reflectivity field
and computes the radar reflectivity for the Thompson scheme.*/
void ThompsonScheme::compute_radar_reflectivity(
    const Field3D& qc,
    const Field3D& qr,
    const Field3D& qi,
    const Field3D& qs,
    const Field3D& qg,
    const Field3D& qh,
    Field3D& reflectivity_dbz
) 
{
    // Iterate over all grid points and compute radar reflectivity
    int NR = qc.size_r();
    if (NR == 0) return;
    int NTH = qc.size_th();
    int NZ = qc.size_z();

    reflectivity_dbz.resize(NR, NTH, NZ, 0.0f);

    // Thompson reflectivity constants
    const float K_qc = 0.1e-3f;
    const float K_qr = 4.0e-3f;
    const float K_qi = 0.1e-3f;
    const float K_qs = 1.0e-3f;
    const float K_qg = 1.0e-2f;
    const float K_qh = 2.0e-2f;
    const float alpha = 1.5f;
    const float Z_min = 1e-10f;

    // Iterate over all grid points and compute radar reflectivity
    for (int i = 0; i < NR; ++i) 
    {
        // Iterate over all azimuthal angles
        for (int j = 0; j < NTH; ++j) 
        {
            // Iterate over all vertical levels
            for (int k = 0; k < NZ; ++k) 
            {
                // Get current values
                float qc_val = static_cast<float>(qc[i][j][k]);
                float qr_val = static_cast<float>(qr[i][j][k]);
                float qi_val = static_cast<float>(qi[i][j][k]);
                float qs_val = static_cast<float>(qs[i][j][k]);
                float qg_val = static_cast<float>(qg[i][j][k]);
                float qh_val = static_cast<float>(qh[i][j][k]);

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

//Helper functions


/*This function computes the freezing probability for the Thompson scheme.
Takes in the temperature and cloud water volume and computes the freezing probability 
for the Thompson scheme.*/
double ThompsonScheme::freezing_probability(double T, double qc_vol) 
{
    // Thompson freezing probability (simplified)
    double T_celsius = T - microphysics_constants::T0;
    if (T_celsius >= 0.0) return 0.0;

    double exp_arg = -120.0 * qc_vol * 5.2e-4 * (std::exp(0.2 * T_celsius) - 1.0);
    return 1.0 - std::exp(exp_arg);
}


/*This function computes the ice size threshold for the Thompson scheme.
Takes in the ice number concentration and computes the ice size threshold 
for the Thompson scheme.*/
double ThompsonScheme::ice_size_threshold(double Ni) 
{
    // Simplified 200 µm rule - should be based on PSD
    return 200e-6;  // Fixed threshold for now
}

/*This function computes the collection efficiency for the Thompson scheme.
Takes in the temperature, ice mixing ratio, and rainwater mixing ratio
and computes the collection efficiency for the Thompson scheme.*/
double ThompsonScheme::collection_efficiency(double T_celsius, double qi, double qr) 
{
    // Simplified collection efficiency
    if (T_celsius < -5.0) return 0.8;
    else if (T_celsius < 0.0) return 0.6;
    else return 0.0;
}
