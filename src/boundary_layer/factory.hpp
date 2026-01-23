#pragma once
#include <memory>
#include <string>
#include "../../include/boundary_layer_base.hpp"

// Forward declarations of scheme classes
class SlabScheme;
class YSUScheme;
class MYNNScheme;

// Factory function to create boundary layer schemes
std::unique_ptr<BoundaryLayerSchemeBase> create_boundary_layer_scheme(const std::string& scheme_name);

// Registry of available schemes
std::vector<std::string> get_available_boundary_layer_schemes();
