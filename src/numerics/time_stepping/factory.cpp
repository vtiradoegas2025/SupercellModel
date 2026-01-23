#include "factory.hpp"
#include "schemes/rk3/rk3.hpp"
#include "schemes/rk4/rk4.hpp"

/*This function creates the time stepping scheme.
Takes in the scheme name and creates the time stepping scheme.*/
std::unique_ptr<TimeSteppingSchemeBase> create_time_stepping_scheme(const std::string& scheme_name) 
{
    if (scheme_name == "rk3") 
    {
        return std::make_unique<RK3Scheme>();
    }
    else if (scheme_name == "rk4") 
    {
        return std::make_unique<RK4Scheme>();
    }
    else 
    {
        throw std::runtime_error("Unknown time stepping scheme: " + scheme_name);
    }
}

/*This function gets the available time stepping schemes.
Takes in the available time stepping schemes.*/
std::vector<std::string> get_available_time_stepping_schemes() 
{
    return {"rk3", "rk4"};
}
