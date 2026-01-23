#pragma once
#include "../../include/chaos_base.hpp"
#include <vector>
#include <string>
#include <memory>

/*This header file contains the declaration of the create_chaos_scheme function.
This function creates a chaos perturbation scheme based on the name of the scheme.
It is used to create the chaos scheme for the simulation.
*/

#include "schemes/none/none.hpp"
#include "schemes/initial_conditions/initial_conditions.hpp"
#include "schemes/boundary_layer/boundary_layer.hpp"
#include "schemes/full_stochastic/full_stochastic.hpp"

// Factory function declaration
std::unique_ptr<chaos::ChaosScheme> create_chaos_scheme(const std::string& scheme_name);

// Registry of available schemes
std::vector<std::string> get_available_chaos_schemes();
