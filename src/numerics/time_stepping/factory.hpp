#pragma once
#include <memory>
#include <string>
#include "../../../include/time_stepping_base.hpp"

// Forward declarations of scheme classes
class RK3Scheme;
class RK4Scheme;

// Factory function to create time stepping schemes
std::unique_ptr<TimeSteppingSchemeBase> create_time_stepping_scheme(const std::string& scheme_name);

// Registry of available schemes
std::vector<std::string> get_available_time_stepping_schemes();
