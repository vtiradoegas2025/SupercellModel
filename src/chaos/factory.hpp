/**
 * @file factory.hpp
 * @brief Declarations for the chaos module.
 *
 * Defines interfaces, data structures, and contracts used by
 * the chaos runtime and scheme implementations.
 * This file is part of the src/chaos subsystem.
 */

#pragma once
#include "chaos_base.hpp"
#include <vector>
#include <string>
#include <memory>


#include "schemes/none/none.hpp"
#include "schemes/initial_conditions/initial_conditions.hpp"
#include "schemes/boundary_layer/boundary_layer.hpp"
#include "schemes/full_stochastic/full_stochastic.hpp"

/**
 * @brief Creates a chaos scheme by configured name.
 */
std::unique_ptr<chaos::ChaosScheme> create_chaos_scheme(const std::string& scheme_name);

/**
 * @brief Returns names of available chaos schemes.
 */
std::vector<std::string> get_available_chaos_schemes();
