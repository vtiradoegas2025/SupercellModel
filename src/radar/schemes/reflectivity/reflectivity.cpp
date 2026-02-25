/**
 * @file reflectivity.cpp
 * @brief Implementation for the radar module.
 *
 * Provides executable logic for the radar runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/radar subsystem.
 */

#include "reflectivity.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>

namespace
{
constexpr float kZeMinLinear = 1.0e-10f;
constexpr float kZeMaxLinear = 1.0e12f;

void ensure_field_shape(Field3D& field, int nr, int nth, int nz, float init_value = 0.0f)
{
    if (field.size_r() != nr || field.size_th() != nth || field.size_z() != nz)
    {
        field.resize(nr, nth, nz, init_value);
        return;
    }
    field.fill(init_value);
}

void prepare_reflectivity_output(RadarOut& out, int nr, int nth, int nz)
{
    out.NR = nr;
    out.NTH = nth;
    out.NZ = nz;

    ensure_field_shape(out.Ze_linear, nr, nth, nz, 0.0f);
    ensure_field_shape(out.Ze_rain, nr, nth, nz, 0.0f);
    ensure_field_shape(out.Ze_snow, nr, nth, nz, 0.0f);
    ensure_field_shape(out.Ze_graupel, nr, nth, nz, 0.0f);
    ensure_field_shape(out.Ze_hail, nr, nth, nz, 0.0f);
    ensure_field_shape(out.Ze_ice, nr, nth, nz, 0.0f);
    ensure_field_shape(out.Z_dBZ, nr, nth, nz, 0.0f);
    ensure_field_shape(out.ZH_dBZ, nr, nth, nz, 0.0f);
    ensure_field_shape(out.ZV_dBZ, nr, nth, nz, 0.0f);
}

bool has_matching_dimensions(const RadarStateView& state, int nr, int nth, int nz)
{
    return state.NR == nr && state.NTH == nth && state.NZ == nz;
}
}

/**
 * @brief Initializes the reflectivity scheme.
 */

const ReflectivityScheme::HydrometeorProps ReflectivityScheme::rain_props = {
    1000.0,
    842.0,
    0.8,
    0.0
};

/**
 * @brief Initializes the snow properties.
 */
const ReflectivityScheme::HydrometeorProps ReflectivityScheme::snow_props = {
    100.0,
    2.0e6,
    0.25,
    0.0
};

/**
 * @brief Initializes the graupel properties.
 */
const ReflectivityScheme::HydrometeorProps ReflectivityScheme::graupel_props = {
    400.0,
    4.0e5,
    0.5,
    0.0
};

/**
 * @brief Initializes the hail properties.
 */
const ReflectivityScheme::HydrometeorProps ReflectivityScheme::hail_props = {
    900.0,
    2.0e4,
    0.5,
    0.0
};

/**
 * @brief Initializes the ice properties.
 */
const ReflectivityScheme::HydrometeorProps ReflectivityScheme::ice_props = {
    917.0,
    2.0e6,
    0.25,
    0.0
};

/**
 * @brief Initializes the reflectivity scheme.
 */
void ReflectivityScheme::initialize(const RadarConfig& config, int NR, int NTH, int NZ) 
{
    NR_ = NR;
    NTH_ = NTH;
    NZ_ = NZ;

    std::cout << "Initialized reflectivity radar scheme (" << NR << "x" << NTH << "x" << NZ << ")" << std::endl;
    std::cout << "Operator tier: " << config.operator_tier << std::endl;
}

/**
 * @brief Computes the reflectivity.
 */
void ReflectivityScheme::compute(const RadarConfig& config, const RadarStateView& state, RadarOut& out) 
{
    prepare_reflectivity_output(out, NR_, NTH_, NZ_);

    if (!has_matching_dimensions(state, NR_, NTH_, NZ_))
    {
        static bool warned_dimension_mismatch = false;
        if (!warned_dimension_mismatch)
        {
            std::cerr << "Warning: Reflectivity scheme state/output dimensions mismatch. "
                      << "Returning zeroed reflectivity output." << std::endl;
            warned_dimension_mismatch = true;
        }
        return;
    }

    if (config.operator_tier == "fast_da") 
    {
        compute_fast_da(config, state, out);
    }
     else if (config.operator_tier == "psd_moment") {
        compute_psd_moment(config, state, out);
    } 
    else 
    {
        throw std::runtime_error("Unknown reflectivity operator tier: " + config.operator_tier);
    }
}

/**
 * @brief Computes the fast DA reflectivity.
 */
void ReflectivityScheme::compute_fast_da(const RadarConfig& config, const RadarStateView& state, RadarOut& out) {

    if (state.qr && config.has_qr) 
    {
        const auto* Nr_ptr = config.has_Nr ? state.Nr : nullptr;
        compute_species_reflectivity(*state.qr, Nr_ptr, rain_props, out.Ze_rain);
    }

    if (state.qs && config.has_qs) 
    {
        const auto* Ns_ptr = config.has_Ns ? state.Ns : nullptr;
        compute_species_reflectivity(*state.qs, Ns_ptr, snow_props, out.Ze_snow);
    }

    if (state.qg && config.has_qg) 
    {
        const auto* Ng_ptr = config.has_Ng ? state.Ng : nullptr;
        compute_species_reflectivity(*state.qg, Ng_ptr, graupel_props, out.Ze_graupel);
    }

    if (state.qh && config.has_qh) 
    {
        const auto* Nh_ptr = config.has_Nh ? state.Nh : nullptr;
        compute_species_reflectivity(*state.qh, Nh_ptr, hail_props, out.Ze_hail);
    }

    if (state.qi && config.has_qi) 
    {
        const auto* Ni_ptr = config.has_Ni ? state.Ni : nullptr;
        compute_species_reflectivity(*state.qi, Ni_ptr, ice_props, out.Ze_ice);
    }

    for (int i = 0; i < NR_; ++i) 
    {
        for (int j = 0; j < NTH_; ++j) 
        {
            for (int k = 0; k < NZ_; ++k) 
            {
                float Ze_total = out.Ze_rain[i][j][k] + out.Ze_snow[i][j][k] +
                                out.Ze_graupel[i][j][k] + out.Ze_hail[i][j][k] +
                                out.Ze_ice[i][j][k];
                if (!std::isfinite(static_cast<double>(Ze_total)))
                {
                    Ze_total = kZeMinLinear;
                }
                Ze_total = std::clamp(Ze_total, kZeMinLinear, kZeMaxLinear);

                out.Ze_linear[i][j][k] = Ze_total;
                out.Z_dBZ[i][j][k] = linear_to_dbz(Ze_total);

                out.ZH_dBZ[i][j][k] = out.Z_dBZ[i][j][k];
                out.ZV_dBZ[i][j][k] = out.Z_dBZ[i][j][k];
            }
        }
    }
}

/**
 * @brief Computes the PSD moment reflectivity.
 */
void ReflectivityScheme::compute_psd_moment(const RadarConfig& config, const RadarStateView& state, RadarOut& out) 
{

    if (state.qr && config.has_qr && config.has_Nr) 
    {
        const auto* Nr_ptr = state.Nr;
        compute_moment_reflectivity(*state.qr, Nr_ptr, rain_props, out.Ze_rain, "rain");
    }
 
    if (state.qs && config.has_qs && config.has_Ns) 
    {
        const auto* Ns_ptr = state.Ns;
        compute_moment_reflectivity(*state.qs, Ns_ptr, snow_props, out.Ze_snow, "snow");
    }

    if (state.qg && config.has_qg && config.has_Ng) 
    {
        const auto* Ng_ptr = state.Ng;
        compute_moment_reflectivity(*state.qg, Ng_ptr, graupel_props, out.Ze_graupel, "graupel");
    }

    if (state.qh && config.has_qh && config.has_Nh) 
    {
        const auto* Nh_ptr = state.Nh;
        compute_moment_reflectivity(*state.qh, Nh_ptr, hail_props, out.Ze_hail, "hail");
    }

    if (state.qi && config.has_qi && config.has_Ni) 
    {
        const auto* Ni_ptr = state.Ni;
        compute_moment_reflectivity(*state.qi, Ni_ptr, ice_props, out.Ze_ice, "ice");
    }

    if (state.qr && config.has_qr && !config.has_Nr) 
    {
        compute_species_reflectivity(*state.qr, nullptr, rain_props, out.Ze_rain);
    }

    if (state.qs && config.has_qs && !config.has_Ns) 
    {
        compute_species_reflectivity(*state.qs, nullptr, snow_props, out.Ze_snow);
    }

    if (state.qg && config.has_qg && !config.has_Ng) 
    {
        compute_species_reflectivity(*state.qg, nullptr, graupel_props, out.Ze_graupel);
    }

    if (state.qh && config.has_qh && !config.has_Nh) 
    {
        compute_species_reflectivity(*state.qh, nullptr, hail_props, out.Ze_hail);
    }

    if (state.qi && config.has_qi && !config.has_Ni) 
    {
        compute_species_reflectivity(*state.qi, nullptr, ice_props, out.Ze_ice);
    }

    for (int i = 0; i < NR_; ++i)
    {
        for (int j = 0; j < NTH_; ++j) 
        {
            for (int k = 0; k < NZ_; ++k) 
            {
                float Ze_total = out.Ze_rain[i][j][k] + out.Ze_snow[i][j][k] +
                                out.Ze_graupel[i][j][k] + out.Ze_hail[i][j][k] +
                                out.Ze_ice[i][j][k];
                if (!std::isfinite(static_cast<double>(Ze_total)))
                {
                    Ze_total = kZeMinLinear;
                }
                Ze_total = std::clamp(Ze_total, kZeMinLinear, kZeMaxLinear);

                out.Ze_linear[i][j][k] = Ze_total;
                out.Z_dBZ[i][j][k] = linear_to_dbz(Ze_total);

                out.ZH_dBZ[i][j][k] = out.Z_dBZ[i][j][k];
                out.ZV_dBZ[i][j][k] = out.Z_dBZ[i][j][k];
            }
        }
    }
}

/**
 * @brief Computes the species reflectivity.
 */
void ReflectivityScheme::compute_species_reflectivity(
    const Field3D& q,
    const Field3D* Nt,
    const HydrometeorProps& props,
    Field3D& Ze_out) {

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

    for (int i = 0; i < NR_; ++i) 
    {
        for (int j = 0; j < NTH_; ++j) 
        {
            for (int k = 0; k < NZ_; ++k) 
            {
                float q_val = q[i][j][k];
                if (!std::isfinite(static_cast<double>(q_val)))
                {
                    q_val = 0.0f;
                }

                if (q_val <= 0.0f) {
                    Ze_out[i][j][k] = 0.0f;
                    continue;
                }


                double rho_ratio = props.density / 1000.0;
                double lambda_4 = std::pow(0.1, 4.0);

                double Ze = 3.63e9 * rho_ratio * rho_ratio * q_val * q_val;

                if (!std::isfinite(Ze))
                {
                    Ze = static_cast<double>(kZeMinLinear);
                }
                Ze = std::clamp(Ze, static_cast<double>(kZeMinLinear), static_cast<double>(kZeMaxLinear));

                Ze_out[i][j][k] = static_cast<float>(Ze);
            }
        }
    }
}

/**
 * @brief Computes the moment reflectivity.
 */
void ReflectivityScheme::compute_moment_reflectivity(
    const Field3D& q,
    const Field3D* Nt,
    const HydrometeorProps& props,
    Field3D& Ze_out,
    const std::string& species_name) {

    if (!Nt) 
    {
        std::cerr << "Warning: compute_moment_reflectivity called without number concentrations for " << species_name << std::endl;
        compute_species_reflectivity(q, nullptr, props, Ze_out);
        return;
    }

    for (int i = 0; i < NR_; ++i) 
    {
        for (int j = 0; j < NTH_; ++j) 
        {
            for (int k = 0; k < NZ_; ++k) 
            {
                float q_val = q[i][j][k];
                float Nt_val = (*Nt)[i][j][k];
                if (!std::isfinite(static_cast<double>(q_val)))
                {
                    q_val = 0.0f;
                }
                if (!std::isfinite(static_cast<double>(Nt_val)))
                {
                    Nt_val = 0.0f;
                }

                if (q_val <= 0.0f || Nt_val <= 0.0f) 
                {
                    Ze_out[i][j][k] = 0.0f;
                    continue;
                }


                double rho_air = 1.2;
                double q_mass = q_val * rho_air;

                double Nt_m3 = Nt_val;

                double rho_ratio = props.density / 1000.0;

                double C = 3.6e9;

                double Ze = C * rho_ratio * rho_ratio * q_mass * q_mass / Nt_m3;

                if (!std::isfinite(Ze))
                {
                    Ze = static_cast<double>(kZeMinLinear);
                }
                Ze = std::clamp(Ze, static_cast<double>(kZeMinLinear), static_cast<double>(kZeMaxLinear));

                Ze_out[i][j][k] = static_cast<float>(Ze);
            }
        }
    }
}
