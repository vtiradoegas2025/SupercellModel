/**
 * @file factory.cpp
 * @brief Implementation for the boundary_layer module.
 *
 * Provides executable logic for the boundary_layer runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/boundary_layer subsystem.
 */

#include "factory.hpp"
#include "schemes/slab/slab.hpp"
#include "schemes/ysu/ysu.hpp"
#include "schemes/mynn/mynn.hpp"
#include "string_utils.hpp"

namespace
{
std::string normalize_boundary_layer_scheme_name(std::string scheme_name)
{
    scheme_name = tmv::strutil::lower_copy(scheme_name);

    if (scheme_name == "yonsei" ||
        scheme_name == "yonsei_university")
    {
        return "ysu";
    }
    if (scheme_name == "mixed_layer" ||
        scheme_name == "mixed-layer")
    {
        return "slab";
    }
    if (scheme_name == "mynn2" ||
        scheme_name == "mellor_yamada_nakanishi_niino")
    {
        return "mynn";
    }
    return scheme_name;
}
}

/**
 * @brief Creates a boundary layer scheme based on the scheme name.
 */
std::unique_ptr<BoundaryLayerSchemeBase> create_boundary_layer_scheme(const std::string& scheme_name) 
{
    const std::string normalized_name = normalize_boundary_layer_scheme_name(scheme_name);

    if (normalized_name == "slab") 
    {
        return std::make_unique<SlabScheme>();
    }

    else if (normalized_name == "ysu") 
    {
        return std::make_unique<YSUScheme>();
    }

    else if (normalized_name == "mynn") 
    {
        return std::make_unique<MYNNScheme>();
    }
    else 
    {
        throw std::runtime_error("Unknown boundary layer scheme: " + scheme_name +
                                 " (normalized: " + normalized_name + ")");
    }
}

/**
 * @brief Returns the available boundary layer schemes.
 */
std::vector<std::string> get_available_boundary_layer_schemes() 
{
    return {"slab", "ysu", "mynn"};
}
