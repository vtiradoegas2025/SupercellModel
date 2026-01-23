#include "factory.hpp"
#include "schemes/kessler/kessler.hpp"
#include "schemes/lin/lin.hpp"
#include "schemes/thompson/thompson.hpp"
#include "schemes/milbrandt/milbrandt.hpp"

// Factory function to create microphysics schemes
std::unique_ptr<MicrophysicsScheme> create_microphysics_scheme(const std::string& scheme_name) 
{
    if (scheme_name == "kessler") 
    {
        return std::make_unique<KesslerScheme>();
    } else if (scheme_name == "lin") 
    {
        return std::make_unique<LinScheme>();
    } 
    else if (scheme_name == "thompson") 
    {
        return std::make_unique<ThompsonScheme>();
    } 
    else if (scheme_name == "milbrandt")
    {
        return std::make_unique<MilbrandtScheme>();
    } 
    else 
    {
        throw std::runtime_error("Unknown microphysics scheme: " + scheme_name);
    }
}

// Registry of available schemes
std::vector<std::string> get_available_schemes() 
{
    return {"kessler", "lin", "thompson", "milbrandt"};
}
