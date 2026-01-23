#pragma once
#include <memory>
#include <string>
#include "../../include/microphysics_base.hpp"

// Forward declarations of scheme classes
class KesslerScheme;
class LinScheme;
class ThompsonScheme;
class MilbrandtScheme;

// Factory function to create microphysics schemes
std::unique_ptr<MicrophysicsScheme> create_microphysics_scheme(const std::string& scheme_name);

// Registry of available schemes
std::vector<std::string> get_available_schemes();
