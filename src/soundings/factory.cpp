#include "factory.hpp"
#include "schemes/sharpy/sharpy_sounding.hpp"
#include <iostream>
#include <stdexcept>

/*This function creates the SHARPY sounding scheme.
Takes in the SHARPY sounding scheme and creates the SHARPY sounding scheme.*/
std::unique_ptr<SoundingScheme> create_sharpy_sounding_scheme() 
{
    return std::make_unique<SharpySoundingScheme>();
}

// Implementation of the factory function declared in soundings_base.hpp
std::unique_ptr<SoundingScheme> create_sounding_scheme(const std::string& scheme_id) 
{
    // If the scheme id is "sharpy" or "SHARPY", create the SHARPY sounding scheme.
    if (scheme_id == "sharpy" || scheme_id == "SHARPY") 
    {
        return create_sharpy_sounding_scheme();
    }

    // If the scheme id is "none" or the scheme id is empty, return null pointer.   
    else if (scheme_id == "none" || scheme_id.empty()) {
        // Return null pointer for "none" scheme
        return nullptr;
    } 
    else 
    {
        throw std::runtime_error("Unknown sounding scheme: " + scheme_id +
                               ". Available schemes: 'sharpy', 'none'");
    }
}
