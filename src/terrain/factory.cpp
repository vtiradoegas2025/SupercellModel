#include "factory.hpp"
#include "schemes/bell/bell.hpp"
#include "schemes/schar/schar.hpp"
#include "schemes/none.hpp"

/*This file contains the implementation of the terrain scheme factory.
It manages the creation of the terrain schemes.*/
std::unique_ptr<TerrainSchemeBase> create_terrain_scheme(const std::string& scheme_name) 
{
    if (scheme_name == "bell") 
    {
        return std::make_unique<BellScheme>();
    }
    else if (scheme_name == "schar") 
    {
        return std::make_unique<ScharScheme>();
    }
    else if (scheme_name == "none") {
        return std::make_unique<NoneScheme>();
    }
    else 
    {
        return std::make_unique<NoneScheme>();  // default
    }
}

/*This function gets the available terrain schemes.
Takes in the available terrain schemes and returns the available terrain schemes.*/
std::vector<std::string> get_available_terrain_schemes() 
{
    return {"none", "bell", "schar"};
}
