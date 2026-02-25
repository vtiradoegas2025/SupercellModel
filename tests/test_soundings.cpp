#include "soundings/base/soundings_base.hpp"
#include "soundings/factory.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

namespace
{

SoundingData make_required_profile(std::size_t levels)
{
    SoundingData data;
    data.height_m.reserve(levels);
    data.pressure_hpa.reserve(levels);
    data.temperature_k.reserve(levels);

    for (std::size_t i = 0; i < levels; ++i)
    {
        data.height_m.push_back(static_cast<double>(i) * 250.0);
        data.pressure_hpa.push_back(1000.0 - static_cast<double>(i) * 35.0);
        data.temperature_k.push_back(300.0 - static_cast<double>(i) * 2.0);
    }
    return data;
}

void test_missing_dewpoint_drops_moisture_derived_fields()
{
    SoundingData data = make_required_profile(6);
    calculate_derived_quantities(data);

    assert(data.potential_temperature_k.size() == data.num_levels());
    assert(data.mixing_ratio_kgkg.empty());
    assert(data.equivalent_potential_temperature_k.empty());
}

void test_quality_control_discards_incomplete_optional_profiles()
{
    SoundingData data = make_required_profile(6);
    data.dewpoint_k = {294.0, 292.0, 290.0};
    data.wind_speed_ms = {5.0, 7.0, 9.0, 11.0, 13.0, 15.0};
    data.wind_direction_deg = {180.0, 190.0};

    SoundingConfig config;
    const bool ok = quality_control_sounding(data, config);
    assert(ok);
    assert(data.dewpoint_k.empty());
    assert(data.wind_speed_ms.empty());
    assert(data.wind_direction_deg.empty());
}

void test_interpolation_handles_missing_optional_profiles()
{
    const SoundingData source = make_required_profile(6);
    const std::vector<double> target_heights = {0.0, 125.0, 250.0, 375.0, 500.0};

    SoundingConfig config;
    SoundingData result = interpolate_sounding_linear(source, target_heights, config);

    assert(result.is_valid());
    assert(result.height_m.size() == target_heights.size());
    assert(result.pressure_hpa.size() == target_heights.size());
    assert(result.temperature_k.size() == target_heights.size());
    assert(result.potential_temperature_k.size() == target_heights.size());
    assert(result.dewpoint_k.empty());
    assert(result.wind_speed_ms.empty());
    assert(result.wind_direction_deg.empty());
    assert(result.mixing_ratio_kgkg.empty());
    assert(result.equivalent_potential_temperature_k.empty());
}

void test_interpolation_single_level_is_stable()
{
    SoundingData source;
    source.height_m = {100.0};
    source.pressure_hpa = {910.0};
    source.temperature_k = {294.0};

    SoundingConfig config;
    config.extrapolate_below_ground = true;
    config.extrapolate_above_top = true;

    const std::vector<double> target_heights = {-50.0, 100.0, 800.0};
    const SoundingData result = interpolate_sounding_linear(source, target_heights, config);
    assert(result.pressure_hpa.size() == target_heights.size());
    assert(result.temperature_k.size() == target_heights.size());
    for (std::size_t i = 0; i < target_heights.size(); ++i)
    {
        assert(std::abs(result.pressure_hpa[i] - 910.0) < 1e-9);
        assert(std::abs(result.temperature_k[i] - 294.0) < 1e-9);
    }
}

void test_factory_normalizes_scheme_id()
{
    {
        std::unique_ptr<SoundingScheme> sharpy = create_sounding_scheme("  ShArPy  ");
        assert(static_cast<bool>(sharpy));
    }
    {
        std::unique_ptr<SoundingScheme> none = create_sounding_scheme(" NONE ");
        assert(!none);
    }

    bool threw = false;
    try
    {
        (void)create_sounding_scheme("unknown-scheme");
    }
    catch (const std::runtime_error&)
    {
        threw = true;
    }
    assert(threw);
}

} // namespace

int main()
{
    test_missing_dewpoint_drops_moisture_derived_fields();
    test_quality_control_discards_incomplete_optional_profiles();
    test_interpolation_handles_missing_optional_profiles();
    test_interpolation_single_level_is_stable();
    test_factory_normalizes_scheme_id();

    std::cout << "Soundings regression tests passed." << std::endl;
    return 0;
}
