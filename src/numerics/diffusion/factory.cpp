#include "factory.hpp"
#include "schemes/explicit/explicit.hpp"
#include "schemes/implicit/implicit.hpp"

/*This function creates the diffusion scheme.
Takes in the scheme name and creates the diffusion scheme.*/
std::unique_ptr<DiffusionSchemeBase> create_diffusion_scheme(const std::string& scheme_name) {
    if (scheme_name == "explicit") 
    {
        return std::make_unique<ExplicitDiffusionScheme>();
    }
    else if (scheme_name == "implicit") 
    {
        return std::make_unique<ImplicitDiffusionScheme>();
    }
    else
    {
        throw std::runtime_error("Unknown diffusion scheme: " + scheme_name);
    }
}

/*This function gets the available diffusion schemes.
Takes in the available diffusion schemes.*/
std::vector<std::string> get_available_diffusion_schemes() 
{
    return {"explicit", "implicit"};
}
