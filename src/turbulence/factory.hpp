#pragma once
#include <memory>
#include <string>
#include "../../include/turbulence_base.hpp"

// Forward declarations of scheme classes
class SmagorinskyScheme;
class TKEScheme;

// Factory function to create turbulence schemes
std::unique_ptr<TurbulenceSchemeBase> create_turbulence_scheme(const std::string& scheme_name);

// Registry of available schemes
std::vector<std::string> get_available_turbulence_schemes();
