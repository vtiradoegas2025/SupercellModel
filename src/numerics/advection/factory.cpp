#include "factory.hpp"
#include "schemes/tvd/tvd.hpp"
#include "schemes/weno5/weno5.hpp"

/*This function creates the advection scheme.
Takes in the scheme name and creates the advection scheme.*/
std::unique_ptr<AdvectionSchemeBase> create_advection_scheme(const std::string& scheme_name) 
{
    if (scheme_name == "tvd") 
    {
        return std::make_unique<TVDScheme>();
    }
    else if (scheme_name == "weno5") 
    {
        return std::make_unique<WENO5Scheme>();
    }
    else 
    {
        throw std::runtime_error("Unknown advection scheme: " + scheme_name);
    }
}

std::vector<std::string> get_available_advection_schemes() 
{
    return {"tvd", "weno5"};
}
