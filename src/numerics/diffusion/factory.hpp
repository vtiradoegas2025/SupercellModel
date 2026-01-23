#pragma once
#include <memory>
#include <string>
#include "../../../include/diffusion_base.hpp"

// Forward declarations of scheme classes
class ExplicitDiffusionScheme;
class ImplicitDiffusionScheme;

// Factory function to create diffusion schemes
std::unique_ptr<DiffusionSchemeBase> create_diffusion_scheme(const std::string& scheme_name);

// Registry of available schemes
std::vector<std::string> get_available_diffusion_schemes();
