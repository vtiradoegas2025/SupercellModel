#include "factory.hpp"

/*This file contains the implementation of the dynamics scheme factory.
This file contains the implementation of the create_dynamics_scheme and get_available_dynamics_schemes functions.*/
std::unique_ptr<DynamicsScheme> create_dynamics_scheme(const std::string& scheme_name) 
{
    if (scheme_name == "supercell") 
    {
        return std::make_unique<SupercellScheme>();
    } 
    else if (scheme_name == "tornado") 
    {
        return std::make_unique<TornadoScheme>();
    } 
    else 
    {
        throw std::runtime_error("Unknown dynamics scheme: " + scheme_name);
    }
}

// Registry of available schemes
std::vector<std::string> get_available_dynamics_schemes() 
{
    return {"supercell", "tornado"};
}
