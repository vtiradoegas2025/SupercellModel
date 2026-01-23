#include "factory.hpp"
#include "schemes/simple_grey/simple_grey.hpp"
// #include "schemes/rrtmg/rrtmg.hpp"  // Future

/*This function creates the radiation scheme.
Takes in the scheme name and creates the radiation scheme. Comments out the RRTMG scheme for now.*/
std::unique_ptr<RadiationSchemeBase> create_radiation_scheme(const std::string& scheme_name) 
{
    if (scheme_name == "simple_grey") {
        return std::make_unique<SimpleGreyScheme>();
    }
    // else if (scheme_name == "rrtmg") {
    //     return std::make_unique<RRTMGScheme>();
    // }
    else 
    {
        throw std::runtime_error("Unknown radiation scheme: " + scheme_name);
    }
}

/*This function gets the available radiation schemes.
Takes in the available radiation schemes and returns the available radiation schemes.*/
std::vector<std::string> get_available_radiation_schemes() 
{
    return {"simple_grey"};  // , "rrtmg"
}
