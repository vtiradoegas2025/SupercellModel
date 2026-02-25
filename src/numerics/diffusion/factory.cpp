/**
 * @file factory.cpp
 * @brief Implementation for the numerics module.
 *
 * Provides executable logic for the numerics runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/numerics subsystem.
 */

#include "factory.hpp"
#include "schemes/explicit/explicit.hpp"
#include "schemes/implicit/implicit.hpp"
#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace
{
std::string canonicalize_scheme_name(std::string value)
{
    std::string canonical;
    canonical.reserve(value.size());
    for (unsigned char c : value)
    {
        if (c == '_' || c == '-' || std::isspace(c))
        {
            continue;
        }
        canonical.push_back(static_cast<char>(std::tolower(c)));
    }
    return canonical;
}
}

/**
 * @brief Creates the diffusion scheme.
 */
std::unique_ptr<DiffusionSchemeBase> create_diffusion_scheme(const std::string& scheme_name) {
    const std::string canonical = canonicalize_scheme_name(scheme_name);

    if (canonical == "explicit") 
    {
        return std::make_unique<ExplicitDiffusionScheme>();
    }
    else if (canonical == "implicit") 
    {
        return std::make_unique<ImplicitDiffusionScheme>();
    }
    else
    {
        throw std::runtime_error("Unknown diffusion scheme: " + scheme_name);
    }
}

/**
 * @brief Gets the available diffusion schemes.
 */
std::vector<std::string> get_available_diffusion_schemes() 
{
    return {"explicit", "implicit"};
}
