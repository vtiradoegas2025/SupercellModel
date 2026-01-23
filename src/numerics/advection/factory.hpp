#pragma once
#include <memory>
#include <string>
#include "../../../include/advection_base.hpp"

// Forward declarations of scheme classes
class TVDScheme;
class WENO5Scheme;

// Factory function to create advection schemes
std::unique_ptr<AdvectionSchemeBase> create_advection_scheme(const std::string& scheme_name);

// Registry of available schemes
std::vector<std::string> get_available_advection_schemes();
