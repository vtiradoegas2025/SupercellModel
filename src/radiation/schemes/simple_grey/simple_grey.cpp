#include "simple_grey.hpp"
#include <iostream>
#include <algorithm>

/*This file contains the implementation of the simple grey radiation scheme.*/

/*This constructor initializes the simple grey radiation scheme.
Takes in the reference optical depths and the exponents and initializes the simple grey radiation scheme.*/
SimpleGreyScheme::SimpleGreyScheme()
    : tau_lw_ref_(6.0), tau_sw_ref_(0.22), n_lw_(4.0), n_sw_(2.0) 
{}

/*This function returns the required fields.
Takes in the required fields and returns the required fields.*/
int SimpleGreyScheme::required_fields() const 
{
    return static_cast<int>(RadiationRequirements::TEMPERATURE) |
           static_cast<int>(RadiationRequirements::SOLAR_GEOMETRY);
}

/*This function initializes the simple grey radiation scheme.
Takes in the configuration and initializes the simple grey radiation scheme.*/
void SimpleGreyScheme::initialize(const RadiationConfig& cfg)
{
    tau_lw_ref_ = cfg.tau_lw_ref;
    tau_sw_ref_ = cfg.tau_sw_ref;
    n_lw_ = cfg.n_lw;
    n_sw_ = cfg.n_sw;

    std::cout << "Initialized Simple Grey Radiation:" << std::endl;
    std::cout << "  LW optical depth: " << tau_lw_ref_ << " (exponent: " << n_lw_ << ")" << std::endl;
    std::cout << "  SW optical depth: " << tau_sw_ref_ << " (exponent: " << n_sw_ << ")" << std::endl;
}

/*This function computes the column of the simple grey radiation scheme.
Takes in the configuration, the column state, the tendencies, and the fluxes 
and computes the column of the simple grey radiation scheme.*/
void SimpleGreyScheme::compute_column(
    const RadiationConfig& cfg,
    const RadiationColumnStateView& col,
    RadiationColumnTendencies& tend,
    RadiationColumnFluxes* fluxes_opt
) 
{
    const size_t nz = col.rho->size();

    // Resize output arrays
    tend.dTdt_rad.resize(nz, 0.0);
    tend.dTdt_lw.resize(nz, 0.0);
    tend.dTdt_sw.resize(nz, 0.0);

    // If the fluxes are requested, resize the output arrays.
    if (fluxes_opt) 
    {
        fluxes_opt->Fup_lw.resize(nz + 1, 0.0);
        fluxes_opt->Fdn_lw.resize(nz + 1, 0.0);
        fluxes_opt->Fup_sw.resize(nz + 1, 0.0);
        fluxes_opt->Fdn_sw.resize(nz + 1, 0.0);
        fluxes_opt->Fnet_lw.resize(nz + 1, 0.0);
        fluxes_opt->Fnet_sw.resize(nz + 1, 0.0);
        fluxes_opt->Fnet_total.resize(nz + 1, 0.0);
    }

    // If the longwave is requested, compute the longwave heating.
    if (cfg.do_lw) 
    {
        compute_lw_heating(col, tend.dTdt_lw, fluxes_opt);
    }

    // If the shortwave is requested, compute the shortwave heating.
    if (cfg.do_sw) 
    {
        compute_sw_heating(col, tend.dTdt_sw, fluxes_opt);
    }

    //Iterate over the vertical levels and compute the total heating.
    for (size_t k = 0; k < nz; ++k) 
    {
        tend.dTdt_rad[k] = tend.dTdt_lw[k] + tend.dTdt_sw[k];
    }

    // If the fluxes are requested, compute the net fluxes.
    if (fluxes_opt) 
    {
        // Iterate over the vertical levels and compute the net fluxes.
        for (size_t k = 0; k <= nz; ++k) 
        {
            fluxes_opt->Fnet_lw[k] = fluxes_opt->Fup_lw[k] - fluxes_opt->Fdn_lw[k];
            fluxes_opt->Fnet_sw[k] = fluxes_opt->Fup_sw[k] - fluxes_opt->Fdn_sw[k];
            fluxes_opt->Fnet_total[k] = fluxes_opt->Fnet_lw[k] + fluxes_opt->Fnet_sw[k];
        }
    }
}

/*This function computes the longwave heating.
Takes in the column state, the tendencies, and the fluxes and computes the longwave heating.*/
void SimpleGreyScheme::compute_lw_heating(
    const RadiationColumnStateView& col,
    std::vector<double>& dTdt_lw,
    RadiationColumnFluxes* fluxes_opt
) 
{
    // Get the number of vertical levels.
    const size_t nz = col.rho->size();

    /*This function computes the optical depth profile.
    Takes in the pressure, the reference pressure, and the optical depth and computes the optical depth profile.*/
    tau_lw_profile_ = radiative_transfer::compute_optical_depth_profile(
        *col.p, tau_lw_ref_, n_lw_
    );

    /*This function computes the layer optical depths.
    Takes in the optical depth profile, and the vertical
    grid spacing and computes the layer optical depths.*/
    auto tau_layers = radiative_transfer::compute_layer_optical_depths(
        tau_lw_profile_, *col.dz
    );

    /*This function computes the grey longwave two stream solution.
    Takes in the temperature, the optical depth layers, the surface temperature, 
    the emissivity, the upward flux, and the downward flux and computes the
    grey longwave two stream solution.*/
    std::vector<double> Fup_lw, Fdn_lw;
    radiative_transfer::grey_lw_two_stream(
        *col.T, tau_layers, col.Tsfc, col.emissivity_lw,
        Fup_lw, Fdn_lw
    );

    // Compute net flux
    std::vector<double> Fnet_lw(nz + 1);

    // Iterate over the vertical levels and compute the net flux.
    for (size_t k = 0; k <= nz; ++k) 
    {
        Fnet_lw[k] = Fup_lw[k] - Fdn_lw[k];
    }

    
    /*This function computes the heating rate from the flux divergence.
    Takes in the density, the vertical grid spacing, the net flux, and the heating rate and computes the heating rate from the flux divergence.*/
    radiative_transfer::heating_rate_from_flux_divergence(
        *col.rho, *col.dz, Fnet_lw, dTdt_lw
    );

    //If the fluxes are requested, store the fluxes.
    if (fluxes_opt) 
    {
        fluxes_opt->Fup_lw = Fup_lw;
        fluxes_opt->Fdn_lw = Fdn_lw;
    }
}

/*This function computes the shortwave heating.
Takes in the column state, the tendencies, and 
the fluxes and computes the shortwave heating.*/
void SimpleGreyScheme::compute_sw_heating(
    const RadiationColumnStateView& col,
    std::vector<double>& dTdt_sw,
    RadiationColumnFluxes* fluxes_opt
) 
{
    const size_t nz = col.rho->size();

    // If the cosine of the zenith angle is less than or equal to 0, set the heating rate to 0(sun is below the horizon).
    if (col.mu0 <= 0.0) 
    {
        std::fill(dTdt_sw.begin(), dTdt_sw.end(), 0.0);
        return;
    }

    /*This function computes the optical depth profile.
    Takes in the pressure, the reference pressure, and the optical depth and computes the optical depth profile.*/
    tau_sw_profile_ = radiative_transfer::compute_optical_depth_profile(
        *col.p, tau_sw_ref_, n_sw_
    );

    /*This function computes the layer optical depths.
    Takes in the optical depth profile, and the vertical
    grid spacing and computes the layer optical depths.*/
    auto tau_layers = radiative_transfer::compute_layer_optical_depths(
        tau_sw_profile_, *col.dz
    );

    /*This function computes the grey shortwave beer-lambert solution.
    Takes in the optical depth layers, the cosine of the zenith angle, the solar constant, 
    the surface albedo, the upward flux, and the downward flux and computes the
    grey shortwave beer-lambert solution.*/
    std::vector<double> Fup_sw, Fdn_sw;
    radiative_transfer::grey_sw_beer_lambert(
        tau_layers, col.mu0, col.S0, col.albedo_sw,
        Fup_sw, Fdn_sw
    );

    /*This function computes the net flux.
    Takes in the optical depth layers, the cosine of the zenith angle, the solar constant, 
    the surface albedo, the upward flux, and the downward flux and computes the
    net flux.*/
    std::vector<double> Fnet_sw(nz + 1);

    // Iterate over the vertical levels and compute the net flux.
    for (size_t k = 0; k <= nz; ++k) 
    {
        Fnet_sw[k] = Fup_sw[k] - Fdn_sw[k];
    }

    /*This function computes the heating rate from the flux divergence.
    Takes in the density, the vertical grid spacing, the net flux, and the heating rate and computes the heating rate from the flux divergence.*/
    radiative_transfer::heating_rate_from_flux_divergence(
        *col.rho, *col.dz, Fnet_sw, dTdt_sw
    );

    // If the fluxes are requested, store the fluxes.
    if (fluxes_opt)
    {
        fluxes_opt->Fup_sw = Fup_sw;
        fluxes_opt->Fdn_sw = Fdn_sw;
    }
}
