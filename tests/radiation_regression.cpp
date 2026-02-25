#include "radiation/base/radiative_transfer.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace
{

bool nearly_equal(double a, double b, double tol = 1.0e-12)
{
    return std::abs(a - b) <= tol;
}

int expect_true(bool cond, const std::string& message)
{
    if (!cond)
    {
        std::cerr << "[radiation-regression] FAIL: " << message << std::endl;
        return 1;
    }
    return 0;
}

int expect_close(double actual, double expected, const std::string& label, double tol = 1.0e-12)
{
    if (!nearly_equal(actual, expected, tol))
    {
        std::cerr << "[radiation-regression] FAIL: " << label
                  << " actual=" << actual
                  << " expected=" << expected
                  << " tol=" << tol << std::endl;
        return 1;
    }
    return 0;
}

int test_optical_depth_profile_and_layers()
{
    int failures = 0;

    const std::vector<double> p = {100000.0, 80000.0, 60000.0, 40000.0};
    const std::vector<double> dz = {1000.0, 1000.0, 1000.0, 1000.0};
    const double tau_ref = 6.0;
    const double n = 1.0;

    const std::vector<double> tau_profile =
        radiative_transfer::compute_optical_depth_profile(p, tau_ref, n);
    failures += expect_true(tau_profile.size() == p.size() + 1, "tau_profile must be nz+1");

    // Expected interfaces from implementation:
    // [p0, avg(p0,p1), avg(p1,p2), avg(p2,p3), 0]
    failures += expect_close(tau_profile[0], 6.0, "tau_profile[0]");
    failures += expect_close(tau_profile[1], 5.4, "tau_profile[1]");
    failures += expect_close(tau_profile[2], 4.2, "tau_profile[2]");
    failures += expect_close(tau_profile[3], 3.0, "tau_profile[3]");
    failures += expect_close(tau_profile[4], 0.0, "tau_profile[4]");

    for (size_t k = 1; k < tau_profile.size(); ++k)
    {
        failures += expect_true(
            tau_profile[k] <= tau_profile[k - 1] + 1.0e-12,
            "tau_profile must be monotonic non-increasing");
    }

    const std::vector<double> dtau =
        radiative_transfer::compute_layer_optical_depths(tau_profile, dz);
    failures += expect_true(dtau.size() == dz.size(), "dtau must be nz");
    failures += expect_close(dtau[0], 0.6, "dtau[0]");
    failures += expect_close(dtau[1], 1.2, "dtau[1]");
    failures += expect_close(dtau[2], 1.2, "dtau[2]");
    failures += expect_close(dtau[3], 3.0, "dtau[3]");

    double dtau_sum = 0.0;
    for (const double value : dtau)
    {
        failures += expect_true(value >= -1.0e-12, "dtau must be non-negative");
        dtau_sum += value;
    }
    failures += expect_close(
        dtau_sum,
        tau_profile.front() - tau_profile.back(),
        "sum(dtau) must equal tau(surface)-tau(top)");

    return failures;
}

int test_shortwave_layer_attenuation()
{
    int failures = 0;

    const std::vector<double> tau_layers = {0.2, 0.3};
    const double mu0 = 0.5;
    const double S0 = 1000.0;
    const double albedo = 0.2;
    std::vector<double> Fup;
    std::vector<double> Fdn;

    radiative_transfer::grey_sw_beer_lambert(tau_layers, mu0, S0, albedo, Fup, Fdn);

    failures += expect_true(Fup.size() == 3, "Fup size must be nz+1");
    failures += expect_true(Fdn.size() == 3, "Fdn size must be nz+1");

    const double expected_Fdn2 = mu0 * S0;
    const double expected_Fdn1 = expected_Fdn2 * std::exp(-tau_layers[1] / mu0);
    const double expected_Fdn0 = expected_Fdn1 * std::exp(-tau_layers[0] / mu0);
    failures += expect_close(Fdn[2], expected_Fdn2, "Fdn[2]");
    failures += expect_close(Fdn[1], expected_Fdn1, "Fdn[1]");
    failures += expect_close(Fdn[0], expected_Fdn0, "Fdn[0]");

    const double expected_Fup0 = albedo * expected_Fdn0;
    const double expected_Fup1 = expected_Fup0 * std::exp(-tau_layers[0] / mu0);
    const double expected_Fup2 = expected_Fup1 * std::exp(-tau_layers[1] / mu0);
    failures += expect_close(Fup[0], expected_Fup0, "Fup[0]");
    failures += expect_close(Fup[1], expected_Fup1, "Fup[1]");
    failures += expect_close(Fup[2], expected_Fup2, "Fup[2]");

    return failures;
}

int test_flux_divergence_shape_guard()
{
    int failures = 0;
    const std::vector<double> rho = {1.0, 0.9};
    const std::vector<double> dz = {100.0, 100.0};
    const std::vector<double> bad_Fnet = {100.0, 90.0};  // wrong size: should be nz+1
    std::vector<double> dTdt;
    radiative_transfer::heating_rate_from_flux_divergence(rho, dz, bad_Fnet, dTdt);
    failures += expect_true(dTdt.size() == rho.size(), "dTdt size must match rho size on guard path");
    failures += expect_close(dTdt[0], 0.0, "guard dTdt[0]");
    failures += expect_close(dTdt[1], 0.0, "guard dTdt[1]");
    return failures;
}

} // namespace

int main()
{
    int failures = 0;
    failures += test_optical_depth_profile_and_layers();
    failures += test_shortwave_layer_attenuation();
    failures += test_flux_divergence_shape_guard();

    if (failures > 0)
    {
        std::cerr << "[radiation-regression] FAILED with " << failures << " check(s)." << std::endl;
        return 1;
    }

    std::cout << "[radiation-regression] all checks passed" << std::endl;
    return 0;
}
