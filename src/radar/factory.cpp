#include "factory.hpp"

// Include scheme implementations
#include "schemes/reflectivity/reflectivity.hpp"
#include "schemes/velocity/velocity.hpp"
#include "schemes/zdr/zdr.hpp"
/*This function creates the radar scheme.
Takes in the scheme id and creates the radar scheme.*/

std::unique_ptr<RadarSchemeBase> RadarFactory::create(const std::string& scheme_id) 
{
    if (scheme_id == "reflectivity") 
    {
        return std::make_unique<ReflectivityScheme>();
    } 
    else if (scheme_id == "velocity") 
    {
        return std::make_unique<VelocityScheme>();
    } 
    else if (scheme_id == "zdr") 
    {
        return std::make_unique<ZDRScheme>();
    } 
    else 
    {
        throw std::runtime_error("Unknown radar scheme: " + scheme_id);
    }
}

/*This function gets the available radar schemes.
Takes in the available radar schemes.*/
std::vector<std::string> RadarFactory::available_schemes() 
{
    return {"reflectivity", "velocity", "zdr"};
}
