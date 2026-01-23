#include "factory.hpp"
#include "schemes/slab/slab.hpp"
#include "schemes/ysu/ysu.hpp"
#include "schemes/mynn/mynn.hpp"

/*This function creates a boundary layer scheme based on the scheme name.
Takes in the scheme name and creates a boundary layer scheme based on the scheme name.
Passes out the boundary layer scheme to the calling function for use in the simulation.*/
std::unique_ptr<BoundaryLayerSchemeBase> create_boundary_layer_scheme(const std::string& scheme_name) 
{
    // If slab
    if (scheme_name == "slab") 
    {
        return std::make_unique<SlabScheme>();
    }

    // If ysu
    else if (scheme_name == "ysu") 
    {
        return std::make_unique<YSUScheme>();
    }

    // If mynn
    else if (scheme_name == "mynn") 
    {
        return std::make_unique<MYNNScheme>();
    }
    else 
    {
        throw std::runtime_error("Unknown boundary layer scheme: " + scheme_name);
    }
}

/*This function returns the available boundary layer schemes.
Takes in no arguments and returns the available boundary layer schemes.
Passes out the available boundary layer schemes to the calling function for use in the simulation.*/
std::vector<std::string> get_available_boundary_layer_schemes() 
{
    return {"slab", "ysu", "mynn"};
}
