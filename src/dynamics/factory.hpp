/**
 * @file factory.hpp
 * @brief Declarations for the dynamics module.
 *
 * Defines interfaces, data structures, and contracts used by
 * the dynamics runtime and scheme implementations.
 * This file is part of the src/dynamics subsystem.
 */

#pragma once
#include "dynamics_base.hpp"
#include <vector>
#include <string>


#include "schemes/supercell/supercell.hpp"
#include "schemes/tornado/tornado.hpp"

/**
 * @brief Creates a dynamics scheme by configured name.
 */
std::unique_ptr<DynamicsScheme> create_dynamics_scheme(const std::string& scheme_name);

/**
 * @brief Returns names of available dynamics schemes.
 */
std::vector<std::string> get_available_dynamics_schemes();
