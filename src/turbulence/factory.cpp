#include "factory.hpp"
#include "schemes/smagorinsky/smagorinsky.hpp"
#include "schemes/tke/tke.hpp"


/*This file contains the implementation of the turbulence scheme factory.
It manages the creation of the turbulence schemes.*/
std::unique_ptr<TurbulenceSchemeBase> create_turbulence_scheme(const std::string& scheme_name) 
{
    if (scheme_name == "smagorinsky") 
    {
        return std::make_unique<SmagorinskyScheme>();
    }
    else if (scheme_name == "tke") 
    {
        return std::make_unique<TKEScheme>();
    }
    else 
    {
        throw std::runtime_error("Unknown turbulence scheme: " + scheme_name);
    }
}

/*This function gets the available turbulence schemes.
Takes in the available turbulence schemes and returns the available turbulence schemes.*/
std::vector<std::string> get_available_turbulence_schemes() 
{
    return {"smagorinsky", "tke"};
}
