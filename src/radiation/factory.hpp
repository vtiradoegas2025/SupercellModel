#pragma once
#include <memory>
#include <string>
#include "../../include/radiation_base.hpp"

// Forward declarations of scheme classes
class SimpleGreyScheme;
// class RRTMGScheme;  // Future implementation

// Factory function to create radiation schemes
std::unique_ptr<RadiationSchemeBase> create_radiation_scheme(const std::string& scheme_name);

// Registry of available schemes
std::vector<std::string> get_available_radiation_schemes();
