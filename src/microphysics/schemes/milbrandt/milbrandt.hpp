#pragma once
#include "../../base/thermodynamics.hpp"
#include "microphysics_base.hpp"
#include "field3d.hpp"

/*This header file contains the declaration of the MilbrandtScheme class.
This class implements the Milbrandt scheme.
It is a subclass of the MicrophysicsScheme class.
It implements the compute_tendencies method.
It implements the compute_radar_reflectivity method.
*/

class MilbrandtScheme : public MicrophysicsScheme 
{
private:
    // Milbrandt-Yau scheme parameters (Milbrandt & Yau 2005)
    double qc0_;        // autoconversion threshold (kg/kg)
    double c_auto_;     // autoconversion rate (s⁻¹)

    // PSD parameters (fixed shape parameters)
    double alpha_r_;    // rain shape parameter
    double alpha_i_;    // cloud ice shape parameter
    double alpha_s_;    // snow shape parameter
    double alpha_g_;    // graupel shape parameter
    double alpha_h_;    // hail shape parameter

    // Mass-diameter relationships: m(D) = c_x * D^{d_x}
    double c_r_, d_r_;  // rain
    double c_i_, d_i_;  // cloud ice
    double c_s_, d_s_;  // snow
    double c_g_, d_g_;  // graupel
    double c_h_, d_h_;  // hail

    // Terminal velocity parameters: Vt(D) = a_x * D^{b_x}
    double a_r_, b_r_;  // rain
    double a_i_, b_i_;  // cloud ice
    double a_s_, b_s_;  // snow
    double a_g_, b_g_;  // graupel
    double a_h_, b_h_;  // hail

    // Triple-moment option (predict Z for graupel/hail)
    bool triple_moment_;
    bool hail_processes_;  // include hail category

    // Prognostic arrays for number concentrations and reflectivity
    Field3D Nr_;  // rain number (m⁻³)
    Field3D Ni_;  // ice number (m⁻³)
    Field3D Ns_;  // snow number (m⁻³)
    Field3D Ng_;  // graupel number (m⁻³)
    Field3D Nh_;  // hail number (m⁻³)

    Field3D Zr_;  // rain reflectivity (mm⁶/m³)
    Field3D Zi_;  // ice reflectivity (mm⁶/m³)
    Field3D Zs_;  // snow reflectivity (mm⁶/m³)
    Field3D Zg_;  // graupel reflectivity (mm⁶/m³)
    Field3D Zh_;  // hail reflectivity (mm⁶/m³)

    // Tendency arrays
    Field3D dNr_dt_, dNi_dt_, dNs_dt_, dNg_dt_, dNh_dt_;
    Field3D dZr_dt_, dZi_dt_, dZs_dt_, dZg_dt_, dZh_dt_;

public:
    /*This constructor initializes the Milbrandt scheme with default parameters.
    Takes in the autoconversion threshold, autoconversion rate, alpha parameters for the rain, ice, snow, graupel, and hail,
    c parameters for the rain, ice, snow, graupel, and hail, d parameters for the rain, ice, snow, graupel, and hail,
    a parameters for the rain, ice, snow, graupel, and hail, b parameters for the rain, ice, snow, graupel, and hail,
    triple moment flag, and hail processes flag and initializes the Milbrandt scheme with the default parameters.*/
    MilbrandtScheme(
        double qc0 = 1.0e-3,
        double c_auto = 1.0e-3,
        double alpha_r = 2.5, double alpha_i = 0.0, double alpha_s = 0.0,
        double alpha_g = 1.0, double alpha_h = 1.0,
        double c_r = 524.0, double d_r = 3.0,
        double c_i = 110.8, double d_i = 2.93,
        double c_s = 2.77e-2, double d_s = 2.0,
        double c_g = 24.0, double d_g = 2.8,
        double c_h = 114.0, double d_h = 2.8,
        double a_r = 65.0, double b_r = 0.125,
        double a_i = 40.0, double b_i = 0.37,
        double a_s = 2.77, double b_s = 0.1,
        double a_g = 114.0, double b_g = 0.5,
        double a_h = 114.0, double b_h = 0.5,
        bool triple_moment = false,
        bool hail_processes = true
    );

    /*This function computes the tendencies for the Milbrandt scheme.
    Takes in the pressure, potential temperature, vapor mixing ratio, cloud water mixing ratio, 
    rainwater mixing ratio, ice mixing ratio, snow mixing ratio, graupel mixing ratio, hail mixing ratio,
    and the time step and computes the tendencies for the Milbrandt scheme.*/
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

    /*This function computes the radar reflectivity for the Milbrandt scheme.
    Takes in the cloud water mixing ratio, rainwater mixing ratio, ice mixing ratio, 
    snow mixing ratio, graupel mixing ratio, hail mixing ratio, and the radar reflectivity 
    fieldand computes the radar reflectivity for the Milbrandt scheme.*/
    void compute_radar_reflectivity(
        const Field3D& qc,
        const Field3D& qr,
        const Field3D& qi,
        const Field3D& qs,
        const Field3D& qg,
        const Field3D& qh,
        Field3D& reflectivity_dbz
    ) override;

    /*This function returns the name of the Milbrandt scheme.
    Returns the name of the Milbrandt scheme.*/
    std::string get_scheme_name() const override { return "milbrandt"; }
    int get_num_prognostic_vars() const override 
    {
        int vars = 7; // qv, qc, qr, qi, qs, qg, qh
        vars += 5;    // Nr, Ni, Ns, Ng, Nh
        if (triple_moment_) vars += 5; // Zr, Zi, Zs, Zg, Zh
        return vars;
    }

    /*This function returns the number of prognostic variables for the Milbrandt scheme.
    It is a subclass of the MicrophysicsScheme class.
    It implements the get_num_prognostic_vars method.
    */
private:
    // PSD and moment calculation functions
    double calculate_lambda(double q, double Nt, double alpha, double c, double d);
    double calculate_N0(double Nt, double alpha, double lambda);
    double calculate_moment(double q, double Nt, double alpha, double c, double d, int moment_order);
    double calculate_reflectivity(double q, double Nt, double alpha, double c, double d);

    // Process rate calculations
    double autoconversion_rate(double qc, double rho);
    double accretion_rate(double q1, double q2, double Nt1, double Nt2, double alpha1, double alpha2,
                         double c1, double d1, double c2, double d2, double rho);
    double sedimentation_rate(double q, double Nt, double alpha, double c, double d, int dz);

    // Collection efficiency calculations
    double rain_selfcollection_rate(double Nr, double rho);
    double ice_collection_efficiency(double T_celsius, double D1, double D2);

    // Helper functions for microphysical processes
    void initialize_prognostic_fields(int NR, int NTH, int NZ);

    /*This function computes the warm rain processes for the Milbrandt scheme.
    Takes in the temperature, pressure, air density, vapor mixing ratio, cloud water mixing ratio, 
    rainwater mixing ratio, and the tendencies for the cloud water, rainwater, and vapor mixing ratios
    and computes the warm rain processes for the Milbrandt scheme.*/
    void compute_warm_rain_processes(
        const Field3D& temperature,
        const Field3D& p,
        const Field3D& rho,
        const Field3D& qv,
        const Field3D& qc,
        const Field3D& qr,
        Field3D& dqc_dt,
        Field3D& dqr_dt,
        Field3D& dqv_dt,
        Field3D& dtheta_dt,
        Field3D& dNr_dt
    );

    /*This function computes the ice processes for the Milbrandt scheme.
    Takes in the temperature, pressure, air density, vapor mixing ratio, cloud water mixing ratio, 
    ice mixing ratio, snow mixing ratio, graupel mixing ratio, hail mixing ratio,
    and the tendencies for the cloud water, ice, snow, graupel, hail, and vapor mixing ratios
    and computes the ice processes for the Milbrandt scheme.*/
    void compute_ice_processes(
        const Field3D& temperature,
        const Field3D& p,
        const Field3D& rho,
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
        Field3D& dNi_dt,
        Field3D& dNs_dt,
        Field3D& dNg_dt,
        Field3D& dNh_dt
    );

    /*This function computes the melting processes for the Milbrandt scheme.
    Takes in the temperature, snow mixing ratio, graupel mixing ratio, hail mixing ratio, 
    and the tendencies for the snow, graupel, hail, and rainwater mixing ratios
    and computes the melting processes for the Milbrandt scheme.*/
    void compute_melting_processes(
        const Field3D& temperature,
        const Field3D& qs,
        const Field3D& qg,
        const Field3D& qh,
        Field3D& dqs_dt,
        Field3D& dqg_dt,
        Field3D& dqh_dt,
        Field3D& dqr_dt,
        Field3D& dtheta_dt,
        Field3D& dNr_dt,
        Field3D& dNs_dt,
        Field3D& dNg_dt,
        Field3D& dNh_dt
    );

    /*This function computes the sedimentation for the Milbrandt scheme.
    Takes in the rainwater mixing ratio, snow mixing ratio, graupel mixing ratio, 
    hail mixing ratio, and the tendencies for the rainwater, snow, graupel, 
    and hail mixing ratios and computes the sedimentation for the Milbrandt scheme.*/
    void compute_sedimentation(
        const Field3D& qr,
        const Field3D& qs,
        const Field3D& qg,
        const Field3D& qh,
        const Field3D& Nr,
        const Field3D& Ns,
        const Field3D& Ng,
        const Field3D& Nh,
        Field3D& dqr_dt,
        Field3D& dqs_dt,
        Field3D& dqg_dt,
        Field3D& dqh_dt,
        Field3D& dNr_dt,
        Field3D& dNs_dt,
        Field3D& dNg_dt,
        Field3D& dNh_dt
    );

    // Gamma function for moment calculations
    double gamma_function(double x);
    double incomplete_gamma(double a, double x);
};
