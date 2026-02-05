#pragma once
#include "../../base/thermodynamics.hpp"
#include "microphysics_base.hpp"
#include "field3d.hpp"


/* This is the Lin scheme class
// It is a derived class from the MicrophysicsScheme base class
// and implements the Lin scheme for the microphysics*/

class LinScheme : public MicrophysicsScheme 
{
private:
    // Lin scheme parameters (Lin et al. 1983)
    double qc0_;        // autoconversion threshold (kg/kg)
    double c_auto_;     // autoconversion rate (s⁻¹)
    double c_accr_;     // accretion coefficient
    double c_evap_;     // evaporation rate (s⁻¹)

    // Ice process parameters
    double c_ihom_;     // homogeneous ice nucleation rate (s⁻¹)
    double c_rime_;     // riming efficiency
    double c_agg_;      // aggregation rate
    double c_melt_;     // melting rate (s⁻¹)
    double c_subl_;     // sublimation rate (s⁻¹)

    // Size distribution parameters (fixed intercept)
    double N0r_;        // rain intercept parameter (m⁻⁴)
    double N0s_;        // snow intercept parameter (m⁻⁴)
    double N0g_;        // graupel intercept parameter (m⁻⁴)
    double N0h_;        // hail intercept parameter (m⁻⁴)

    // Terminal velocity parameters
    double a_r_, b_r_;   // rain: Vt = a_r * (ρq_r)^b_r
    double a_s_, b_s_;   // snow: Vt = a_s * (ρq_s)^b_s
    double a_g_, b_g_;   // graupel: Vt = a_g * (ρq_g)^b_g
    double a_h_, b_h_;   // hail: Vt = a_h * (ρq_h)^b_h

public:

    /*This constructor initializes the Lin scheme with default parameters.
    Takes in the autoconversion threshold, autoconversion rate, accretion coefficient, 
    evaporation rate, homogeneous ice nucleation rate, riming efficiency, aggregation rate, 
    melting rate, sublimation rate, rain intercept parameter, snow intercept parameter, 
    graupel intercept parameter, hail intercept parameter, rain terminal velocity coefficient, 
    rain terminal velocity exponent, snow terminal velocity coefficient, snow terminal velocity exponent, 
    graupel terminal velocity coefficient, graupel terminal velocity exponent, hail terminal velocity coefficient, 
    and hail terminal velocity exponent and initializes the Lin scheme with the default parameters.*/
    LinScheme(
        double qc0 = 1.0e-3,
        double c_auto = 1.0e-3,
        double c_accr = 2.2,
        double c_evap = 1.0e-3,
        double c_ihom = 1.0e-3,
        double c_rime = 1.0,
        double c_agg = 1.0e-3,
        double c_melt = 1.0e-3,
        double c_subl = 1.0e-3,
        double N0r = 8.0e6,
        double N0s = 3.0e6,
        double N0g = 4.0e5,
        double N0h = 1.0e4,
        double a_r = 65.0, double b_r = 0.125,
        double a_s = 40.0, double b_s = 0.37,
        double a_g = 114.0, double b_g = 0.5,
        double a_h = 114.0, double b_h = 0.5
    );

    /*This function computes the tendencies for the Lin scheme.
    Takes in the pressure, potential temperature, vapor mixing ratio, cloud water mixing ratio, 
    rainwater mixing ratio, ice mixing ratio, snow mixing ratio, graupel mixing ratio, hail mixing ratio,
    and the time step and computes the tendencies for the Lin scheme.*/
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

    /*This function computes the radar reflectivity for the Lin scheme.
    Takes in the cloud water mixing ratio, rainwater mixing ratio, ice mixing ratio, snow mixing ratio, 
    graupel mixing ratio, hail mixing ratio, and the radar reflectivity field
    and computes the radar reflectivity for the Lin scheme.*/
    void compute_radar_reflectivity(
        const Field3D& qc,
        const Field3D& qr,
        const Field3D& qi,
        const Field3D& qs,
        const Field3D& qg,
        const Field3D& qh,
        Field3D& reflectivity_dbz
    ) override;

    std::string get_scheme_name() const override { return "lin"; }
    int get_num_prognostic_vars() const override { return 7; }  // qv, qc, qr, qi, qs, qg, qh

private:
    // Helper functions for Lin scheme processes

    /*This function computes the saturation adjustment for the Lin scheme.
    Takes in the temperature, pressure, vapor mixing ratio, and cloud water mixing ratio
    and computes the saturation adjustment for the Lin scheme. Simplified for now.*/
    void saturation_adjustment(
        const Field3D& temperature,
        const Field3D& p,
        Field3D& qv,
        Field3D& qc
    );

    /*This function computes the warm rain processes for the Lin scheme.
    Takes in the temperature, pressure, vapor mixing ratio, cloud water mixing ratio, 
    rainwater mixing ratio, and the tendencies for the cloud water, rainwater, and vapor mixing ratios
    and computes the warm rain processes for the Lin scheme.*/
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

    /*This function computes the ice processes for the Lin scheme.
    Takes in the temperature, pressure, vapor mixing ratio, cloud water mixing ratio, 
    ice mixing ratio, snow mixing ratio, graupel mixing ratio, hail mixing ratio,
    and the tendencies for the cloud water, ice, snow, graupel, hail, and vapor mixing ratios
    and computes the ice processes for the Lin scheme.*/
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
        Field3D& dtheta_dt
    );

    /*This function computes the melting processes for the Lin scheme.
    Takes in the temperature, snow mixing ratio, graupel mixing ratio, hail mixing ratio, 
    and the tendencies for the snow, graupel, hail, and rainwater mixing ratios
    and computes the melting processes for the Lin scheme.*/
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

    /*This function computes the sedimentation for the Lin scheme.
    Takes in the rainwater mixing ratio, snow mixing ratio, graupel mixing ratio, hail mixing ratio, 
    and the tendencies for the rainwater, snow, graupel, and hail mixing ratios
    and computes the sedimentation for the Lin scheme.*/
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

    // Utility functions
    double ice_number_concentration(double qi, double T_celsius);
    double collection_efficiency(double T_celsius);
    double ice_threshold_temperature(double qi);
};
