#include "factory.hpp"

/*This file contains the implementation of the factory function.
This file contains the implementation of the create_chaos_scheme function.*/
namespace chaos 
{

// Factory function to create chaos schemes
std::unique_ptr<ChaosScheme> create_chaos_scheme(const std::string& scheme_name) 
{
    // If the scheme name is none, return a none scheme.
    if (scheme_name == "none") 
    {
        return std::make_unique<NoneScheme>();
    } 

    // If the scheme name is initial_conditions, return an initial conditions scheme.
    else if (scheme_name == "initial_conditions") 
    {
        return std::make_unique<InitialConditionsScheme>();
    } 

    // If the scheme name is boundary_layer, return a boundary layer scheme.
    else if (scheme_name == "boundary_layer") 
    {
        return std::make_unique<BoundaryLayerScheme>();
    } 

    // If the scheme name is full_stochastic, return a full stochastic scheme.
    else if (scheme_name == "full_stochastic") 
    {
        return std::make_unique<FullStochasticScheme>();
    } 

    // If the scheme name is unknown, throw an error.
    else 
    {
        throw std::runtime_error("Unknown chaos scheme: " + scheme_name);
    }
}

} // namespace chaos

/*This function creates the chaos scheme.
Takes in the scheme name and creates the chaos scheme. Global factory function (for external use)*/
std::unique_ptr<chaos::ChaosScheme> create_chaos_scheme(const std::string& scheme_name) 
{
    return chaos::create_chaos_scheme(scheme_name);
}

/*This function gets the available chaos schemes.
Takes in the available chaos schemes and returns the available chaos schemes.*/
std::vector<std::string> get_available_chaos_schemes() 
{
    return {"none", "initial_conditions", "boundary_layer", "full_stochastic"};
}
