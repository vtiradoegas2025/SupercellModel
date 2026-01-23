#include "radar.hpp"
#include "radar/factory.hpp"

/*This file contains the implementation of the radar system.
It manages the initialization of the radar system and the computation of the radar system.*/


/*This function computes all the radar observables.
Takes in the radar state and output and computes all the radar observables.*/
void RadarSystem::compute_all_observables(const RadarStateView& state, RadarOut& output,
                                         double radar_x, double radar_y, double radar_z) {
    // Create and configure all radar schemes
    auto reflectivity = RadarFactory::create(RadarSchemes::REFLECTIVITY);
    auto velocity = RadarFactory::create(RadarSchemes::VELOCITY);
    auto zdr = RadarFactory::create(RadarSchemes::ZDR);

    // Configure reflectivity
    RadarConfig refl_config;
    refl_config.scheme_id = RadarSchemes::REFLECTIVITY;
    refl_config.operator_tier = RadarTiers::Reflectivity::FAST_DA;
    refl_config.has_qr = (state.qr != nullptr);
    refl_config.has_qs = (state.qs != nullptr);
    refl_config.has_qg = (state.qg != nullptr);
    refl_config.has_qh = (state.qh != nullptr);
    refl_config.has_qi = (state.qi != nullptr);

    // Configure velocity
    RadarConfig vel_config;
    vel_config.scheme_id = RadarSchemes::VELOCITY;
    vel_config.radar_x = radar_x;
    vel_config.radar_y = radar_y;
    vel_config.radar_z = radar_z;

    // Configure ZDR
    RadarConfig zdr_config;
    zdr_config.scheme_id = RadarSchemes::ZDR;
    zdr_config.operator_tier = RadarTiers::ZDR::POLARIMETRIC_FO;
    zdr_config.has_qr = (state.qr != nullptr);
    zdr_config.has_qs = (state.qs != nullptr);
    zdr_config.has_qg = (state.qg != nullptr);
    zdr_config.has_qh = (state.qh != nullptr);
    zdr_config.has_qi = (state.qi != nullptr);

    // Initialize schemes
    reflectivity->initialize(refl_config, state.NR, state.NTH, state.NZ);
    velocity->initialize(vel_config, state.NR, state.NTH, state.NZ);
    zdr->initialize(zdr_config, state.NR, state.NTH, state.NZ);

    // Compute observables
    reflectivity->compute(refl_config, state, output);
    velocity->compute(vel_config, state, output);
    zdr->compute(zdr_config, state, output);
}

/*This function computes the radar reflectivity.
Takes in the radar state and output and computes the radar reflectivity.*/
void RadarSystem::compute_reflectivity(const RadarStateView& state, RadarOut& output) 
{
    auto reflectivity = RadarFactory::create(RadarSchemes::REFLECTIVITY);

    RadarConfig config;
    config.scheme_id = RadarSchemes::REFLECTIVITY;
    config.operator_tier = RadarTiers::Reflectivity::FAST_DA;
    config.has_qr = (state.qr != nullptr);
    config.has_qs = (state.qs != nullptr);
    config.has_qg = (state.qg != nullptr);
    config.has_qh = (state.qh != nullptr);
    config.has_qi = (state.qi != nullptr);

    reflectivity->initialize(config, state.NR, state.NTH, state.NZ);
    reflectivity->compute(config, state, output);
}

/*This function computes the radar velocity.
Takes in the radar state and output and computes the radar velocity.*/
void RadarSystem::compute_velocity(const RadarStateView& state, RadarOut& output,
                                  double radar_x, double radar_y, double radar_z) 
                                  {
    auto velocity = RadarFactory::create(RadarSchemes::VELOCITY);

    RadarConfig config;
    config.scheme_id = RadarSchemes::VELOCITY;
    config.radar_x = radar_x;
    config.radar_y = radar_y;
    config.radar_z = radar_z;

    velocity->initialize(config, state.NR, state.NTH, state.NZ);
    velocity->compute(config, state, output);
}

/*This function computes the radar ZDR.
Takes in the radar state and output and computes the radar ZDR.*/
void RadarSystem::compute_zdr(const RadarStateView& state, RadarOut& output) 
{
    auto zdr = RadarFactory::create(RadarSchemes::ZDR);

    RadarConfig config;
    config.scheme_id = RadarSchemes::ZDR;
    config.operator_tier = RadarTiers::ZDR::POLARIMETRIC_FO;
    config.has_qr = (state.qr != nullptr);
    config.has_qs = (state.qs != nullptr);
    config.has_qg = (state.qg != nullptr);
    config.has_qh = (state.qh != nullptr);
    config.has_qi = (state.qi != nullptr);

    zdr->initialize(config, state.NR, state.NTH, state.NZ);
    zdr->compute(config, state, output);
}
