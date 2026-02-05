#pragma once
#include "../../base/thermodynamics.hpp"
#include "microphysics_base.hpp"
#include "field3d.hpp"

/*This is the Thompson scheme class
It is a derived class from the MicrophysicsScheme base class
and implements the Thompson scheme for the microphysics*/
class ThompsonScheme : public MicrophysicsScheme 
{
private:
    // Thompson scheme parameters (Thompson et al. 2008)
    double qc0_;        // autoconversion threshold (kg/kg)
    double c_auto_;     // autoconversion rate (s⁻¹)

    // Ice nucleation parameters
    double ccn_conc_;   // CCN concentration (m⁻³)
    double in_conc_;    // ice nuclei concentration (m⁻³)

    // Terminal velocity parameters
    double a_r_, b_r_;   // rain: Vt = a_r * (ρq_r)^b_r
    double a_s_, b_s_;   // snow: Vt = a_s * (ρq_s)^b_s
    double a_g_, b_g_;   // graupel: Vt = a_g * (ρq_g)^b_g
    double a_h_, b_h_;   // hail: Vt = a_h * (ρq_h)^b_h

    // Prognostic ice number concentration
    Field3D Ni_;  // ice number concentration (m⁻³)

    // Internal arrays for Thompson processes
    Field3D dNi_dt_;  // ice number tendency

public:
    /*This constructor initializes the Thompson scheme with default parameters.
    Takes in the autoconversion threshold, autoconversion rate, cloud number concentration, 
    ice number concentration, a parameters for the rain, b parameters for the rain, a parameters for the snow, 
    b parameters for the snow, a parameters for the graupel, b parameters for the graupel, a parameters for the hail, 
    and b parameters for the hail and initializes the Thompson scheme with the default parameters.*/
    ThompsonScheme(
        double qc0 = 1.0e-3,
        double c_auto = 1.0e-3,
        double ccn_conc = 100e6,    // m⁻³
        double in_conc = 100e3,     // m⁻³
        double a_r = 65.0, double b_r = 0.125,
        double a_s = 40.0, double b_s = 0.37,
        double a_g = 114.0, double b_g = 0.5,
        double a_h = 114.0, double b_h = 0.5
    );

    /*This function computes the tendencies for the Thompson scheme.
    Takes in the pressure, potential temperature, vapor mixing ratio, cloud water mixing ratio, 
    rainwater mixing ratio, ice mixing ratio, snow mixing ratio, graupel mixing ratio, hail mixing ratio,
    and the time step and computes the tendencies for the Thompson scheme.*/
    void compute_tendencies(
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
    ) override;

    /*This function computes the radar reflectivity for the Thompson scheme.
    Takes in the cloud water mixing ratio, rainwater mixing ratio, ice mixing ratio, snow mixing ratio, 
    graupel mixing ratio, hail mixing ratio, and the radar reflectivity field
    and computes the radar reflectivity for the Thompson scheme.*/
    void compute_radar_reflectivity(
        const Field3D& qc,
        const Field3D& qr,
        const Field3D& qi,
        const Field3D& qs,
        const Field3D& qg,
        const Field3D& qh,
        Field3D& reflectivity_dbz
    ) override;

    std::string get_scheme_name() const override { return "thompson"; }
    int get_num_prognostic_vars() const override { return 8; }  // qv, qc, qr, qi, qs, qg, qh, Ni

private:
    // Helper functions for Thompson scheme processes

   /*This function performs the saturation adjustment for the Thompson scheme.
   Takes in the temperature, pressure, vapor mixing ratio, and cloud water mixing ratio
   and performs the saturation adjustment for the Thompson scheme.*/
    void saturation_adjustment(
        const Field3D& temperature,
        const Field3D& p,
        Field3D& qv,
        Field3D& qc
    );

    /*This function computes the warm rain processes for the Thompson scheme.
    Takes in the temperature, pressure, vapor mixing ratio, cloud water mixing ratio, 
    rainwater mixing ratio, and the tendencies for the cloud water, rainwater, and vapor mixing ratios
    and computes the warm rain processes for the Thompson scheme.*/
    void compute_warm_rain_processes(
        const Field3D& temperature,
        const Field3D& p,
        const Field3D& qv,
        const Field3D& qc,
        const Field3D& qr,
        Field3D& dqc_dt,
        Field3D& dqr_dt,
        Field3D& dqv_dt,
        Field3D& dtheta_dt
    );

    /*This function computes the ice processes for the Thompson scheme.
    Takes in the temperature, pressure, vapor mixing ratio, cloud water mixing ratio, 
    ice mixing ratio, snow mixing ratio, graupel mixing ratio, hail mixing ratio,
    and the tendencies for the cloud water, ice, snow, graupel, hail, and vapor mixing ratios
    and computes the ice processes for the Thompson scheme.*/
    void compute_ice_processes(
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
    );

    /*This function computes the melting processes for the Thompson scheme.
    Takes in the temperature, snow mixing ratio, graupel mixing ratio, hail mixing ratio, 
    and the tendencies for the snow, graupel, hail, and rainwater mixing ratios
    and computes the melting processes for the Thompson scheme.*/
    void compute_melting_processes(
        const Field3D& temperature,
        const Field3D& qs,
        const Field3D& qg,
        const Field3D& qh,
        Field3D& dqs_dt,
        Field3D& dqg_dt,
        Field3D& dqh_dt,
        Field3D& dqr_dt,
        Field3D& dtheta_dt
    );

    /*This function computes the sedimentation for the Thompson scheme.
    Takes in the rainwater mixing ratio, snow mixing ratio, graupel mixing ratio, hail mixing ratio, 
    and the tendencies for the rainwater, snow, graupel, and hail mixing ratios
    and computes the sedimentation for the Thompson scheme.*/
    void compute_sedimentation(
        const Field3D& qr,
        const Field3D& qs,
        const Field3D& qg,
        const Field3D& qh,
        Field3D& dqr_dt,
        Field3D& dqs_dt,
        Field3D& dqg_dt,
        Field3D& dqh_dt
    );

    // Thompson-specific utility functions
    double freezing_probability(double T, double qc_vol);  // P freezing
    double ice_size_threshold(double Ni);  // 200 µm rule for ice->snow conversion
    double collection_efficiency(double T_celsius, double qi, double qr);  // E_xy
};
