#pragma once
#include "../../base/thermodynamics.hpp"
#include "microphysics_base.hpp"

/*This class implements the Kessler scheme for microphysics.
This class implements the Kessler scheme for microphysics.*/
class KesslerScheme : public MicrophysicsScheme 
{
private:
    // Kessler microphysics parameters (Kessler 1969)
    double qc0_;        // autoconversion threshold (kg/kg)
    double c_auto_;     // autoconversion rate (s⁻¹)
    double c_accr_;     // accretion coefficient (s⁻¹)
    double c_evap_;     // evaporation rate (s⁻¹)

    // Ice microphysics parameters (extended Kessler)
    double c_freeze_;   // homogeneous freezing rate (s⁻¹)
    double c_rime_;     // riming efficiency
    double c_melt_;     // melting rate (s⁻¹)
    double c_subl_;     // sublimation rate (s⁻¹)

    // Terminal velocity parameters
    double a_term_;     // rain terminal velocity coefficient
    double b_term_;     // rain terminal velocity exponent
    double Vt_max_;     // maximum terminal velocity (m/s)
    double a_hail_;     // hail terminal velocity coefficient
    double b_hail_;     // hail terminal velocity exponent
    double a_grau_;     // graupel terminal velocity coefficient
    double b_grau_;     // graupel terminal velocity exponent
    double Vt_max_hail_; // max hail terminal velocity (m/s)
    double Vt_max_grau_; // max graupel terminal velocity (m/s)

public:
    // Constructor with default Kessler parameters
    KesslerScheme(
        double qc0 = 1.0e-3,
        double c_auto = 1.0e-3,
        double c_accr = 2.2,
        double c_evap = 3.0e-3,
        double c_freeze = 1.0e-3,
        double c_rime = 1.0,
        double c_melt = 1.0e-3,
        double c_subl = 1.0e-3,
        double a_term = 65.0,
        double b_term = 0.125,
        double Vt_max = 20.0,
        double a_hail = 114.0,
        double b_hail = 0.5,
        double a_grau = 40.0,
        double b_grau = 0.37,
        double Vt_max_hail = 40.0,
        double Vt_max_grau = 15.0
    );

    /*This function computes the tendencies for the Kessler scheme.
    Takes in the pressure, potential temperature, vapor mixing ratio, cloud water mixing ratio, 
    rainwater mixing ratio, ice mixing ratio, snow mixing ratio, graupel mixing ratio, hail mixing ratio,
    and the time step and computes the tendencies for the Kessler scheme.*/
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

    /*This function computes the radar reflectivity for the Kessler scheme.
    Takes in the cloud water mixing ratio, rainwater mixing ratio, ice mixing ratio, snow mixing ratio, 
    graupel mixing ratio, hail mixing ratio, and the radar reflectivity field
    and computes the radar reflectivity for the Kessler scheme.*/
    void compute_radar_reflectivity(
        const Field3D& qc,
        const Field3D& qr,
        const Field3D& qi,
        const Field3D& qs,
        const Field3D& qg,
        const Field3D& qh,
        Field3D& reflectivity_dbz
    ) override;

    std::string get_scheme_name() const override { return "kessler"; }
    int get_num_prognostic_vars() const override { return 5; }  // qv, qc, qr, qg, qh

private:
    
    /*This function computes the warm rain processes for the Kessler scheme.
    Takes in the temperature, vapor mixing ratio, cloud water mixing ratio, 
    rainwater mixing ratio, and the tendencies for the cloud water, rainwater, and vapor mixing ratios
    and computes the warm rain processes for the Kessler scheme.*/
    void compute_warm_rain_processes(
        const Field3D& temperature,
        const Field3D& qv,
        const Field3D& qc,
        const Field3D& qr,
        Field3D& dqc_dt,
        Field3D& dqr_dt,
        Field3D& dqv_dt,
        Field3D& dtheta_dt
    );

    /*This function computes the ice processes for the Kessler scheme.
    Takes in the temperature, vapor mixing ratio, cloud water mixing ratio, 
    rainwater mixing ratio, ice mixing ratio, snow mixing ratio, graupel mixing ratio, hail mixing ratio,
    and the tendencies for the cloud water, graupel, hail, and vapor mixing ratios
    and computes the ice processes for the Kessler scheme.*/
    void compute_ice_processes(
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
    );

    /*This function computes the melting processes for the Kessler scheme.
    Takes in the temperature, graupel mixing ratio, hail mixing ratio, and 
    the tendencies for the graupel, hail, and rainwater mixing ratios
    and computes the melting processes for the Kessler scheme.*/
    void compute_melting_processes(
        const Field3D& temperature,
        const Field3D& qg,
        const Field3D& qh,
        Field3D& dqg_dt,
        Field3D& dqh_dt,
        Field3D& dqr_dt,
        Field3D& dtheta_dt
    );

    /*This function computes the sedimentation for the Kessler scheme.
    Takes in the rainwater mixing ratio, graupel mixing ratio, hail mixing ratio, 
    and the tendencies for the rainwater, graupel, and hail mixing ratios
    and computes the sedimentation for the Kessler scheme.*/    
    void compute_sedimentation(
        const Field3D& qr,
        const Field3D& qg,
        const Field3D& qh,
        Field3D& dqr_dt,
        Field3D& dqg_dt,
        Field3D& dqh_dt
    );
};
