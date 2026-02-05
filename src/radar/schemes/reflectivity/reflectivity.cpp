#include "reflectivity.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

/*This function initializes the reflectivity scheme.
Takes in the configuration, the number of rows, the number of columns, 
and the number of vertical levels and initializes the reflectivity scheme.*/

const ReflectivityScheme::HydrometeorProps ReflectivityScheme::rain_props = {
    1000.0,  // density (kg/m³)
    842.0,   // alpha (intercept parameter, m⁻⁴)
    0.8,     // c (size distribution parameter)
    0.0      // d (not used for rain)
};

/*This function initializes the snow properties.
Takes in the snow properties.*/
const ReflectivityScheme::HydrometeorProps ReflectivityScheme::snow_props = {
    100.0,   // density (kg/m³)
    2.0e6,   // alpha (intercept parameter)
    0.25,    // c
    0.0      // d
};

/*This function initializes the graupel properties.
Takes in the graupel properties.*/
const ReflectivityScheme::HydrometeorProps ReflectivityScheme::graupel_props = {
    400.0,   // density (kg/m³)
    4.0e5,   // alpha
    0.5,     // c
    0.0      // d
};

/*This function initializes the hail properties.
Takes in the hail properties.*/
const ReflectivityScheme::HydrometeorProps ReflectivityScheme::hail_props = {
    900.0,   // density (kg/m³)
    2.0e4,   // alpha
    0.5,     // c
    0.0      // d
};

/*This function initializes the ice properties.
Takes in the ice properties.*/
const ReflectivityScheme::HydrometeorProps ReflectivityScheme::ice_props = {
    917.0,   // density (kg/m³)
    2.0e6,   // alpha
    0.25,    // c
    0.0      // d
};

/*This function initializes the reflectivity scheme.
Takes in the configuration, the number of rows, the number of columns, 
and the number of vertical levels and initializes the reflectivity scheme.*/
void ReflectivityScheme::initialize(const RadarConfig& config, int NR, int NTH, int NZ) 
{
    NR_ = NR;
    NTH_ = NTH;
    NZ_ = NZ;

    std::cout << "Initialized reflectivity radar scheme (" << NR << "x" << NTH << "x" << NZ << ")" << std::endl;
    std::cout << "Operator tier: " << config.operator_tier << std::endl;
}

/*This function computes the reflectivity.
Takes in the configuration, the state, and the output and computes the reflectivity.*/
void ReflectivityScheme::compute(const RadarConfig& config, const RadarStateView& state, RadarOut& out) 
{
    // If the operator tier is fast DA, compute the fast DA reflectivity.
    if (config.operator_tier == "fast_da") 
    {
        compute_fast_da(config, state, out);
    }
    // If the operator tier is PSD moment, compute the PSD moment reflectivity.
     else if (config.operator_tier == "psd_moment") {
        compute_psd_moment(config, state, out);
    } 
    else 
    {
        throw std::runtime_error("Unknown reflectivity operator tier: " + config.operator_tier);
    }
}

/*This function computes the fast DA reflectivity.
Takes in the configuration, the state, and the output and computes the fast DA reflectivity.*/
void ReflectivityScheme::compute_fast_da(const RadarConfig& config, const RadarStateView& state, RadarOut& out) {
    // Fast DA-style operator: assumes size distributions and computes Z_e directly
    // This is similar to what's currently in the microphysics schemes

    // If the rain is present, compute the rain reflectivity.
    if (state.qr && config.has_qr) 
    {
        const auto* Nr_ptr = config.has_Nr ? state.Nr : nullptr;
        compute_species_reflectivity(*state.qr, Nr_ptr, rain_props, out.Ze_rain);
    }

    // If the snow is present, compute the snow reflectivity.
    if (state.qs && config.has_qs) 
    {
        const auto* Ns_ptr = config.has_Ns ? state.Ns : nullptr;
        compute_species_reflectivity(*state.qs, Ns_ptr, snow_props, out.Ze_snow);
    }

    // If the graupel is present, compute the graupel reflectivity.
    if (state.qg && config.has_qg) 
    {
        const auto* Ng_ptr = config.has_Ng ? state.Ng : nullptr;
        compute_species_reflectivity(*state.qg, Ng_ptr, graupel_props, out.Ze_graupel);
    }

    // If the hail is present, compute the hail reflectivity.
    if (state.qh && config.has_qh) 
    {
        const auto* Nh_ptr = config.has_Nh ? state.Nh : nullptr;
        compute_species_reflectivity(*state.qh, Nh_ptr, hail_props, out.Ze_hail);
    }

    // If the ice is present, compute the ice reflectivity.
    if (state.qi && config.has_qi) 
    {
        const auto* Ni_ptr = config.has_Ni ? state.Ni : nullptr;
        compute_species_reflectivity(*state.qi, Ni_ptr, ice_props, out.Ze_ice);
    }

    // Iterate over the rows, columns, and vertical levels and compute the total linear reflectivity.
    for (int i = 0; i < NR_; ++i) 
    {
        // Iterate over the columns and vertical levels and compute the total linear reflectivity.
        for (int j = 0; j < NTH_; ++j) 
        {
            // Iterate over the vertical levels and compute the total linear reflectivity.
            for (int k = 0; k < NZ_; ++k) 
            {
                float Ze_total = out.Ze_rain[i][j][k] + out.Ze_snow[i][j][k] +
                                out.Ze_graupel[i][j][k] + out.Ze_hail[i][j][k] +
                                out.Ze_ice[i][j][k];

                out.Ze_linear[i][j][k] = Ze_total;
                out.Z_dBZ[i][j][k] = linear_to_dbz(Ze_total);

                // For now, assume Z_H ≈ Z_e (isotropic scatterers)
                out.ZH_dBZ[i][j][k] = out.Z_dBZ[i][j][k];
            }
        }
    }
}

/*This function computes the PSD moment reflectivity.
Takes in the configuration, the state, and the output and computes the PSD moment reflectivity.*/
void ReflectivityScheme::compute_psd_moment(const RadarConfig& config, const RadarStateView& state, RadarOut& out) 
{
    // PSD-moment operator: uses actual particle size distribution moments
    // This provides physically consistent reflectivity when number concentrations are available

    // If the rain is present, compute the rain reflectivity.
    if (state.qr && config.has_qr && config.has_Nr) 
    {
        const auto* Nr_ptr = state.Nr;
        compute_moment_reflectivity(*state.qr, Nr_ptr, rain_props, out.Ze_rain, "rain");
    }
 
    // If the snow is present, compute the snow reflectivity.
    if (state.qs && config.has_qs && config.has_Ns) 
    {
        const auto* Ns_ptr = state.Ns;
        compute_moment_reflectivity(*state.qs, Ns_ptr, snow_props, out.Ze_snow, "snow");
    }

    // If the graupel is present, compute the graupel reflectivity.
    if (state.qg && config.has_qg && config.has_Ng) 
    {
        const auto* Ng_ptr = state.Ng;
        compute_moment_reflectivity(*state.qg, Ng_ptr, graupel_props, out.Ze_graupel, "graupel");
    }

    // If the hail is present, compute the hail reflectivity.
    if (state.qh && config.has_qh && config.has_Nh) 
    {
        const auto* Nh_ptr = state.Nh;
        compute_moment_reflectivity(*state.qh, Nh_ptr, hail_props, out.Ze_hail, "hail");
    }

    // If the ice is present, compute the ice reflectivity.
    if (state.qi && config.has_qi && config.has_Ni) 
    {
        const auto* Ni_ptr = state.Ni;
        compute_moment_reflectivity(*state.qi, Ni_ptr, ice_props, out.Ze_ice, "ice");
    }

    // If the rain is present, but no number concentrations are available, compute the rain reflectivity.       
    if (state.qr && config.has_qr && !config.has_Nr) 
    {
        compute_species_reflectivity(*state.qr, nullptr, rain_props, out.Ze_rain);
    }

    // If the snow is present, but no number concentrations are available, compute the snow reflectivity.   
    if (state.qs && config.has_qs && !config.has_Ns) 
    {
        compute_species_reflectivity(*state.qs, nullptr, snow_props, out.Ze_snow);
    }

    // If the graupel is present, but no number concentrations are available, compute the graupel reflectivity.   
    if (state.qg && config.has_qg && !config.has_Ng) 
    {
        compute_species_reflectivity(*state.qg, nullptr, graupel_props, out.Ze_graupel);
    }

    // If the hail is present, but no number concentrations are available, compute the hail reflectivity.   
    if (state.qh && config.has_qh && !config.has_Nh) 
    {
        compute_species_reflectivity(*state.qh, nullptr, hail_props, out.Ze_hail);
    }

    // If the ice is present, but no number concentrations are available, compute the ice reflectivity.   
    if (state.qi && config.has_qi && !config.has_Ni) 
    {
        compute_species_reflectivity(*state.qi, nullptr, ice_props, out.Ze_ice);
    }

    // Iterate over the rows, columns, and vertical levels and compute the total linear reflectivity.
    for (int i = 0; i < NR_; ++i)
    {
        // Iterate over the columns and vertical levels and compute the total linear reflectivity.
        for (int j = 0; j < NTH_; ++j) 
        {
            // Iterate over the vertical levels and compute the total linear reflectivity.
            for (int k = 0; k < NZ_; ++k) 
            {
                float Ze_total = out.Ze_rain[i][j][k] + out.Ze_snow[i][j][k] +
                                out.Ze_graupel[i][j][k] + out.Ze_hail[i][j][k] +
                                out.Ze_ice[i][j][k];

                out.Ze_linear[i][j][k] = Ze_total;
                out.Z_dBZ[i][j][k] = linear_to_dbz(Ze_total);

                // For now, assume Z_H ≈ Z_e (isotropic scatterers)
                out.ZH_dBZ[i][j][k] = out.Z_dBZ[i][j][k];
            }
        }
    }
}

/*This function computes the species reflectivity.
Takes in the q, the number concentrations, the properties, 
and the output and computes the species reflectivity.*/
void ReflectivityScheme::compute_species_reflectivity(
    const Field3D& q,
    const Field3D* Nt,
    const HydrometeorProps& props,
    Field3D& Ze_out) {

    // If we have number concentrations, use moment-based calculation
    if (Nt) 
    {
        std::string species_name = "unknown";
        if (&props == &rain_props) species_name = "rain";
        else if (&props == &snow_props) species_name = "snow";
        else if (&props == &graupel_props) species_name = "graupel";
        else if (&props == &hail_props) species_name = "hail";
        else if (&props == &ice_props) species_name = "ice";

        compute_moment_reflectivity(q, Nt, props, Ze_out, species_name);
        return;
    }

    // Iterate over the rows, columns, and vertical levels and compute the species reflectivity.
    for (int i = 0; i < NR_; ++i) 
    {
        // Iterate over the columns and vertical levels and compute the species reflectivity.
        for (int j = 0; j < NTH_; ++j) 
        {
            // Iterate over the vertical levels and compute the species reflectivity.
            for (int k = 0; k < NZ_; ++k) 
            {
                float q_val = q[i][j][k];

                if (q_val <= 0.0f) {
                    Ze_out[i][j][k] = 0.0f;
                    continue;
                }

                // Fast DA-style reflectivity calculation
                // Z_e = (N_0 * |K|^2 / lambda^4) * (rho / rho_water)^2 * q^2
                // Simplified from AMS literature for DA applications

                double rho_ratio = props.density / 1000.0;  // normalized by water density
                double lambda_4 = std::pow(0.1, 4.0);       // lambda^4 for S-band

                // Simplified reflectivity calculation (from Kessler/Lin schemes)
                double Ze = 3.63e9 * rho_ratio * rho_ratio * q_val * q_val;

                // Apply minimum threshold
                Ze = std::max(Ze, 1e-10);

                Ze_out[i][j][k] = static_cast<float>(Ze);
            }
        }
    }
}

/*This function computes the moment reflectivity.
Takes in the q, the number concentrations, the properties, 
and the output and computes the moment reflectivity.*/
void ReflectivityScheme::compute_moment_reflectivity(
    const Field3D& q,
    const Field3D* Nt,
    const HydrometeorProps& props,
    Field3D& Ze_out,
    const std::string& species_name) {

    // If the number concentrations are not available, compute the species reflectivity.
    if (!Nt) 
    {
        std::cerr << "Warning: compute_moment_reflectivity called without number concentrations for " << species_name << std::endl;
        compute_species_reflectivity(q, nullptr, props, Ze_out);
        return;
    }

    // Iterate over the rows, columns, and vertical levels and compute the moment reflectivity.
    for (int i = 0; i < NR_; ++i) 
    {
        // Iterate over the columns and vertical levels and compute the moment reflectivity.
        for (int j = 0; j < NTH_; ++j) 
        {
            // Iterate over the vertical levels and compute the moment reflectivity.
            for (int k = 0; k < NZ_; ++k) 
            {
                float q_val = q[i][j][k];
                float Nt_val = (*Nt)[i][j][k];

                // If the quantity or number concentration is less than or equal to 0, set the reflectivity to 0.
                if (q_val <= 0.0f || Nt_val <= 0.0f) 
                {
                    Ze_out[i][j][k] = 0.0f;
                    continue;
                }

                // PSD-moment reflectivity calculation
                // Based on AMS literature for double-moment schemes
                // Z = (720 * lambda^4 / (pi^5 * |K|^2)) * (rho * q / rho_water)^2 / Nt
                // Simplified for computational efficiency

                // Convert mixing ratio to mass content per unit volume
                // Assume air density ~ 1.2 kg/m³ for simplicity (could be improved with actual rho)
                double rho_air = 1.2;  // kg/m³
                double q_mass = q_val * rho_air;  // kg/m³

                // Convert number concentration to m⁻³ (assuming input is in consistent units)
                double Nt_m3 = Nt_val;  // Assume already in m⁻³

                // Rayleigh reflectivity for exponential PSD: Z ∝ (q^2 / Nt) * f(ρ_hydrometeor)
                // From AMS literature (e.g., Milbrandt and Yau, 2005)
                double rho_ratio = props.density / 1000.0;  // hydrometeor density / water density

                // Moment-based reflectivity calculation
                // Z_e = C * (ρ_hyd * q)^2 / N_t * |K|^2 / λ^4
                // where C is a constant depending on PSD assumptions
                double C = 3.6e9;  // Calibration constant for S-band (can be tuned)

                double Ze = C * rho_ratio * rho_ratio * q_mass * q_mass / Nt_m3;

                // Apply minimum threshold
                Ze = std::max(Ze, 1e-10);

                Ze_out[i][j][k] = static_cast<float>(Ze);
            }
        }
    }
}
