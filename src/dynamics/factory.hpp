#pragma once
#include "../../include/dynamics_base.hpp"
#include <vector>
#include <string>

/*This header file contains the declaration of the create_dynamics_scheme function.
This function creates a dynamics scheme based on the name of the scheme.
It is used to create the dynamics scheme for the simulation.
*/

#include "schemes/supercell/supercell.hpp"
#include "schemes/tornado/tornado.hpp"

// Factory function declaration
std::unique_ptr<DynamicsScheme> create_dynamics_scheme(const std::string& scheme_name);

// Registry of available schemes
std::vector<std::string> get_available_dynamics_schemes();
