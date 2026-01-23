#pragma once
#include <memory>
#include <string>
#include "../../include/terrain_base.hpp"

// Forward declarations of scheme classes
class BellScheme;
class ScharScheme;
class NoneScheme;

// Factory function to create terrain schemes
std::unique_ptr<TerrainSchemeBase> create_terrain_scheme(const std::string& scheme_name);

// Registry of available schemes
std::vector<std::string> get_available_terrain_schemes();
