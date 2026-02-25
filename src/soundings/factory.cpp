/**
 * @file factory.cpp
 * @brief Implementation for the soundings module.
 *
 * Provides executable logic for the soundings runtime path,
 * including initialization, stepping, and diagnostics helpers.
 * This file is part of the src/soundings subsystem.
 */

#include "factory.hpp"
#include "schemes/sharpy/sharpy_sounding.hpp"
#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace
{

std::string normalize_scheme_id(std::string value)
{
    auto is_space = [](unsigned char ch) { return std::isspace(ch) != 0; };
    value.erase(value.begin(), std::find_if(value.begin(), value.end(), [&](char ch) {
        return !is_space(static_cast<unsigned char>(ch));
    }));
    value.erase(std::find_if(value.rbegin(), value.rend(), [&](char ch) {
        return !is_space(static_cast<unsigned char>(ch));
    }).base(), value.end());

    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

}

/**
 * @brief Creates the SHARPY sounding scheme.
 */
std::unique_ptr<SoundingScheme> create_sharpy_sounding_scheme() 
{
    return std::make_unique<SharpySoundingScheme>();
}

std::unique_ptr<SoundingScheme> create_sounding_scheme(const std::string& scheme_id) 
{
    const std::string normalized_id = normalize_scheme_id(scheme_id);

    if (normalized_id == "sharpy")
    {
        return create_sharpy_sounding_scheme();
    }

    else if (normalized_id == "none" || normalized_id.empty()) {
        return nullptr;
    } 
    else 
    {
        throw std::runtime_error("Unknown sounding scheme: " + scheme_id +
                               ". Available schemes: 'sharpy', 'none'");
    }
}
