#pragma once
#include "../../base/radiative_transfer.hpp"
#include "radiation_base.hpp"

/*This class implements the simple grey radiation scheme.*/ 
class SimpleGreyScheme : public RadiationSchemeBase 
{
private:
    // Configuration parameters
    double tau_lw_ref_;
    double tau_sw_ref_;
    double n_lw_;
    double n_sw_;

    // Cached optical depth profiles
    std::vector<double> tau_lw_profile_;
    std::vector<double> tau_sw_profile_;

public:
    SimpleGreyScheme();

    std::string name() const override { return "simple_grey"; }
    int required_fields() const override;

    /*This function initializes the simple grey radiation scheme.
    Takes in the configuration and initializes the simple grey radiation scheme.*/
    void initialize(const RadiationConfig& cfg) override;

    /*This function computes the column of the simple grey radiation scheme.
    Takes in the configuration, the column state, the tendencies, and the fluxes 
    and computes the column of the simple grey radiation scheme.*/
    void compute_column(
        const RadiationConfig& cfg,
        const RadiationColumnStateView& col,
        RadiationColumnTendencies& tend,
        RadiationColumnFluxes* fluxes_opt = nullptr) override;

private:

    /*This function computes the longwave heating.
    Takes in the column state, the tendencies, and the fluxes and computes the longwave heating.*/
    void compute_lw_heating(
        const RadiationColumnStateView& col,
        std::vector<double>& dTdt_lw,
        RadiationColumnFluxes* fluxes_opt);

    /*This function computes the shortwave heating.
    Takes in the column state, the tendencies, and the fluxes and computes the shortwave heating.*/
    void compute_sw_heating(
        const RadiationColumnStateView& col,
        std::vector<double>& dTdt_sw,
        RadiationColumnFluxes* fluxes_opt);
};
