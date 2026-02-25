/**
 * @file field_contract.cpp
 * @brief Implementation for the validation module.
 *
 * Provides executable logic for the validation runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/validation subsystem.
 */

#include "field_contract.hpp"
#include "string_utils.hpp"

#include <algorithm>

namespace tmv {
namespace {

/**
 * @brief Builds an inclusive min/max bounds descriptor.
 */
FieldBounds bounds(double min_value, double max_value) {
    FieldBounds out;
    out.has_min = true;
    out.has_max = true;
    out.min_value = min_value;
    out.max_value = max_value;
    return out;
}

const std::vector<FieldContract> kContracts = {
    {"u", "m/s", "Radial wind component", FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow,
     bounds(-180.0, 180.0), {}, {"ur", "u_r"}, FieldRequirementTier::RequiredNow},
    {"v", "m/s", "Azimuthal wind component", FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow,
     bounds(-180.0, 180.0), {}, {"v_theta", "uth"}, FieldRequirementTier::RequiredNow},
    {"w", "m/s", "Vertical wind component", FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow,
     bounds(-120.0, 120.0), {}, {"wz", "vertical_velocity"}, FieldRequirementTier::RequiredNow},
    {"rho", "kg/m^3", "Air density", FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow,
     bounds(0.05, 2.0), {}, {}, FieldRequirementTier::RequiredNow},
    {"p", "Pa", "Pressure", FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow,
     bounds(1000.0, 200000.0), {}, {}, FieldRequirementTier::RequiredNow},
    {"theta", "K", "Potential temperature", FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow,
     bounds(200.0, 500.0), {}, {"potential_temperature"}, FieldRequirementTier::RequiredNow},
    {"qv", "kg/kg", "Water vapor mixing ratio", FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow,
     bounds(0.0, 0.05), {}, {}, FieldRequirementTier::RequiredNow},
    {"qc", "kg/kg", "Cloud water mixing ratio", FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow,
     bounds(0.0, 0.05), {}, {}, FieldRequirementTier::RequiredNow},
    {"qr", "kg/kg", "Rain water mixing ratio", FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow,
     bounds(0.0, 0.05), {}, {}, FieldRequirementTier::RequiredNow},
    {"qh", "kg/kg", "Hail mixing ratio", FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow,
     bounds(0.0, 0.05), {}, {}, FieldRequirementTier::RequiredNow},
    {"qg", "kg/kg", "Graupel mixing ratio", FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow,
     bounds(0.0, 0.05), {}, {}, FieldRequirementTier::RequiredNow},
    {"radar", "mm^6/m^3", "Radar reflectivity (linear units)", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(0.0, 1.0e12), {}, {"radar_reflectivity", "reflectivity_linear"},
     FieldRequirementTier::RequiredNow},
    {"tracer", "1", "Passive tracer", FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow,
     bounds(0.0, 1.5), {}, {}, FieldRequirementTier::RequiredNow},

    {"qi", "kg/kg", "Cloud ice mixing ratio", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(0.0, 0.05), {}, {},
     FieldRequirementTier::RequiredNow},
    {"qs", "kg/kg", "Snow mixing ratio", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(0.0, 0.05), {}, {"snow_mixing_ratio"},
     FieldRequirementTier::RequiredNow},
    {"vorticity_r", "1/s", "Radial vorticity", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(-10.0, 10.0), {}, {}},
    {"vorticity_theta", "1/s", "Azimuthal vorticity", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(-10.0, 10.0), {}, {}},
    {"vorticity_z", "1/s", "Vertical vorticity", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(-10.0, 10.0), {}, {"zeta", "vertical_vorticity"},
     FieldRequirementTier::RequiredNow},
    {"stretching_term", "1/s^2", "Vorticity stretching term", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(-50.0, 50.0), {}, {}, FieldRequirementTier::RequiredNow},
    {"tilting_term", "1/s^2", "Vorticity tilting term", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(-50.0, 50.0), {}, {}, FieldRequirementTier::RequiredNow},
    {"baroclinic_term", "1/s^2", "Baroclinic vorticity generation", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(-50.0, 50.0), {}, {}, FieldRequirementTier::RequiredNow},
    {"p_prime", "Pa", "Perturbation pressure", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(-100000.0, 100000.0), {}, {"pressure_perturbation"},
     FieldRequirementTier::RequiredNow},
    {"dynamic_pressure", "Pa", "Dynamic pressure component", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(-100000.0, 100000.0), {}, {}},
    {"buoyancy_pressure", "Pa", "Buoyancy pressure component", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(-100000.0, 100000.0), {}, {}},
    {"angular_momentum", "m^2/s", "Angular momentum", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, {}, {}, {}},
    {"angular_momentum_tendency", "m^2/s^2", "Angular momentum tendency", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, {}, {}, {}},

    {"horizontal_vorticity_streamwise", "1/s", "Horizontal vorticity (streamwise)",
     FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow, bounds(-50.0, 50.0), {}, {}},
    {"horizontal_vorticity_crosswise", "1/s", "Horizontal vorticity (crosswise)",
     FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow, bounds(-50.0, 50.0), {}, {}},
    {"vorticity_magnitude", "1/s", "Vorticity magnitude", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(0.0, 50.0), {}, {}},
    {"divergence", "1/s", "Flow divergence", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(-5.0, 5.0), {}, {}},
    {"pressure_gradient_force_x", "m/s^2", "Pressure gradient force x-component",
     FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow, bounds(-200.0, 200.0), {}, {}},
    {"pressure_gradient_force_y", "m/s^2", "Pressure gradient force y-component",
     FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow, bounds(-200.0, 200.0), {}, {}},
    {"pressure_gradient_force_z", "m/s^2", "Pressure gradient force z-component",
     FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow, bounds(-200.0, 200.0), {}, {}},
    {"buoyancy", "m/s^2", "Buoyancy", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(-10.0, 10.0), {}, {}},
    {"theta_prime", "K", "Perturbation potential temperature", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(-50.0, 50.0), {}, {}},
    {"theta_v", "K", "Virtual potential temperature", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(180.0, 550.0), {}, {}},
    {"theta_e", "K", "Equivalent potential temperature", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(180.0, 800.0), {}, {}},
    {"theta_w", "K", "Wet-bulb potential temperature", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(180.0, 800.0), {}, {}},
    {"temperature", "K", "Temperature", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(180.0, 360.0), {}, {"t"}},
    {"dewpoint", "K", "Dewpoint temperature", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(150.0, 330.0), {}, {"td"}},
    {"relative_humidity", "%", "Relative humidity", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(0.0, 200.0), {}, {"rh"}},
    {"saturation_mixing_ratio", "kg/kg", "Saturation mixing ratio", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(0.0, 0.10), {}, {"qvs", "qsat"}},
    {"total_condensate", "kg/kg", "Total condensate", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(0.0, 0.10), {}, {"qt"}},
    {"reflectivity_dbz", "dBZ", "Radar reflectivity dBZ", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(-30.0, 120.0), {}, {"z"}},
    {"zdr", "dB", "Differential reflectivity", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(-5.0, 10.0), {}, {"zdr_db"}},
    {"kdp", "deg/km", "Specific differential phase", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(0.0, 60.0), {}, {}},
    {"rhohv", "1", "Correlation coefficient", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(0.0, 1.1), {}, {}},
    {"streamlines", "vector", "Streamline diagnostics", FieldDimensionality::Volume3D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"trajectory_paths", "path", "Trajectory path diagnostics", FieldDimensionality::Trajectory,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"storm_relative_winds", "m/s", "Storm-relative horizontal wind speed (magnitude proxy)",
     FieldDimensionality::Volume3D, FieldImplementationStatus::ExportedNow, bounds(0.0, 200.0), {}, {}},
    {"helicity_density", "m/s^2", "Helicity density", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(-200.0, 200.0), {}, {}},
    {"q_vectors", "vector", "Q-vector diagnostic", FieldDimensionality::Volume3D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"okubo_weiss", "1/s^2", "Okubo-Weiss parameter", FieldDimensionality::Volume3D,
     FieldImplementationStatus::ExportedNow, bounds(-500.0, 500.0), {}, {}},
    {"turbulent_diffusion_term", "varies", "Turbulent diffusion term", FieldDimensionality::Volume3D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},

    {"u10", "m/s", "10 m u wind", FieldDimensionality::Surface2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"v10", "m/s", "10 m v wind", FieldDimensionality::Surface2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"t2", "K", "2 m temperature", FieldDimensionality::Surface2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"td2", "K", "2 m dewpoint", FieldDimensionality::Surface2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"surface_pressure_perturbation", "Pa", "Surface pressure perturbation", FieldDimensionality::Surface2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"surface_sensible_heat_flux", "W/m^2", "Surface sensible heat flux", FieldDimensionality::Surface2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"surface_latent_heat_flux", "W/m^2", "Surface latent heat flux", FieldDimensionality::Surface2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"surface_moisture_flux", "kg/m^2/s", "Surface moisture flux", FieldDimensionality::Surface2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"skin_temperature", "K", "Surface skin temperature", FieldDimensionality::Surface2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"cold_pool_boundary", "K", "Cold pool boundary diagnostic", FieldDimensionality::Surface2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},

    {"composite_reflectivity", "dBZ", "Column max reflectivity", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"column_max_w", "m/s", "Column max vertical velocity", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"column_max_vorticity", "1/s", "Column max vorticity", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"accumulated_rainfall", "mm", "Accumulated rainfall", FieldDimensionality::Surface2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"precip_rate", "mm/h", "Instantaneous precipitation rate", FieldDimensionality::Surface2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"vil", "kg/m^2", "Vertically integrated liquid", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"cloud_top_height", "m", "Cloud top height", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"cloud_base_height", "m", "Cloud base height", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"lcl", "m", "Lifted condensation level", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"lfc", "m", "Level of free convection", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"el", "m", "Equilibrium level", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},

    {"cape", "J/kg", "Convective available potential energy", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"cin", "J/kg", "Convective inhibition", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"lifted_index", "K", "Lifted index", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"k_index", "K", "K-index", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"showalter_index", "K", "Showalter index", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"total_totals", "K", "Total totals index", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"srh_0_1km", "m^2/s^2", "Storm-relative helicity 0-1km", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"srh_0_3km", "m^2/s^2", "Storm-relative helicity 0-3km", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"ehi", "1", "Energy-helicity index", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"scp", "1", "Supercell composite parameter", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"stp", "1", "Significant tornado parameter", FieldDimensionality::Column2D,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},

    {"cross_section", "varies", "Arbitrary cross-section diagnostics", FieldDimensionality::CrossSection,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"rhi_slice", "varies", "RHI-style slice", FieldDimensionality::CrossSection,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"hodograph_aligned_cross_section", "varies", "Along-hodograph cross section", FieldDimensionality::CrossSection,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},

    {"ppi_sweep", "radar", "Simulated PPI sweep", FieldDimensionality::RadarSynthetic,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"rhi_sweep", "radar", "Simulated RHI sweep", FieldDimensionality::RadarSynthetic,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"bwer", "1", "Bounded weak echo region diagnostic", FieldDimensionality::RadarSynthetic,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"mesocyclone_diagnostic", "1", "Mesocyclone detection diagnostic", FieldDimensionality::RadarSynthetic,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"vrot", "m/s", "Rotational velocity", FieldDimensionality::RadarSynthetic,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},

    {"forward_trajectories", "path", "Forward trajectories", FieldDimensionality::Trajectory,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"backward_trajectories", "path", "Backward trajectories", FieldDimensionality::Trajectory,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"parcel_buoyancy_trajectory", "m/s^2", "Parcel buoyancy along trajectory", FieldDimensionality::Trajectory,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"vorticity_trajectory", "1/s", "Vorticity along trajectory", FieldDimensionality::Trajectory,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
    {"circulation_material_surface", "m^2/s", "Circulation on material surfaces", FieldDimensionality::Trajectory,
     FieldImplementationStatus::NotImplemented, {}, {}, {}},
};

}

/**
 * @brief Returns the canonical field contract inventory.
 */
const std::vector<FieldContract>& cm1_field_contracts() {
    return kContracts;
}

/**
 * @brief Finds a field contract by canonical id or known alias.
 */
const FieldContract* find_field_contract(std::string_view id_or_alias) {
    const std::string requested = strutil::lower_copy(id_or_alias);

    for (const auto& contract : kContracts) {
        if (strutil::lower_copy(contract.id) == requested) {
            return &contract;
        }

        for (const auto& alias : contract.aliases) {
            if (strutil::lower_copy(alias) == requested) {
                return &contract;
            }
        }
    }

    return nullptr;
}

/**
 * @brief Returns all contracts matching the provided implementation status.
 */
std::vector<const FieldContract*> contracts_with_status(FieldImplementationStatus status) {
    std::vector<const FieldContract*> out;
    out.reserve(kContracts.size());

    for (const auto& contract : kContracts) {
        if (contract.status == status) {
            out.push_back(&contract);
        }
    }

    return out;
}

/**
 * @brief Returns all contracts matching the provided requirement tier.
 */
std::vector<const FieldContract*> contracts_with_requirement(FieldRequirementTier requirement) {
    std::vector<const FieldContract*> out;
    out.reserve(kContracts.size());

    for (const auto& contract : kContracts) {
        if (contract.requirement == requirement) {
            out.push_back(&contract);
        }
    }

    return out;
}

/**
 * @brief Lists required-now fields that are still marked not implemented.
 */
std::vector<std::string> known_not_implemented_required_now_inventory() {
    std::vector<std::string> out;
    for (const auto& contract : kContracts) {
        if (contract.status == FieldImplementationStatus::NotImplemented &&
            contract.requirement == FieldRequirementTier::RequiredNow) {
            out.push_back(contract.id);
        }
    }
    return out;
}

/**
 * @brief Lists report-only backlog fields that are still not implemented.
 */
std::vector<std::string> known_not_implemented_backlog_inventory() {
    std::vector<std::string> out;
    for (const auto& contract : kContracts) {
        if (contract.status == FieldImplementationStatus::NotImplemented &&
            contract.requirement == FieldRequirementTier::ReportOnly) {
            out.push_back(contract.id);
        }
    }
    return out;
}

/**
 * @brief Converts field dimensionality enum to stable string id.
 */
const char* to_string(FieldDimensionality value) {
    switch (value) {
        case FieldDimensionality::Volume3D:
            return "3d";
        case FieldDimensionality::Surface2D:
            return "surface_2d";
        case FieldDimensionality::Column2D:
            return "column_2d";
        case FieldDimensionality::CrossSection:
            return "cross_section";
        case FieldDimensionality::Trajectory:
            return "trajectory";
        case FieldDimensionality::RadarSynthetic:
            return "radar_synthetic";
        default:
            return "unknown";
    }
}

/**
 * @brief Converts implementation status enum to stable string id.
 */
const char* to_string(FieldImplementationStatus value) {
    switch (value) {
        case FieldImplementationStatus::ExportedNow:
            return "exported_now";
        case FieldImplementationStatus::ComputedNotExported:
            return "computed_not_exported";
        case FieldImplementationStatus::NotImplemented:
            return "not_implemented";
        default:
            return "unknown";
    }
}

/**
 * @brief Converts requirement-tier enum to stable string id.
 */
const char* to_string(FieldRequirementTier value) {
    switch (value) {
        case FieldRequirementTier::RequiredNow:
            return "required_now";
        case FieldRequirementTier::ReportOnly:
            return "report_only";
        default:
            return "unknown";
    }
}

}
