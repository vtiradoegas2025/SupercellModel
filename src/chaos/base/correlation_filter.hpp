#pragma once
#include <vector>
#include <complex>
#include <string>

namespace chaos 
{

//==============================================================================
// Spatial correlation filters for perturbation fields
//==============================================================================

/**
 * @brief Base class for spatial correlation filters
 */
class CorrelationFilter 
{
public:
    virtual ~CorrelationFilter() = default;

    /**
     * @brief Apply correlation filter to 2D field
     * @param field Input/output field (modified in-place)
     * @param dx Grid spacing in x-direction
     * @param dy Grid spacing in y-direction
     */
    virtual     void apply_2d(
        std::vector<std::vector<double>>& field,
        double dx,
        double dy
    ) = 0;

    /**
     * @brief Apply correlation filter to 3D field (horizontal correlation only)
     * @param field Input/output field (modified in-place)
     * @param dx Grid spacing in x-direction
     * @param dy Grid spacing in y-direction
     */
    virtual void apply_3d(
        std::vector<std::vector<std::vector<float>>>& field,
        double dx,
        double dy
    ) 
    {
        // Default implementation: apply 2D filter to each horizontal slice
        if (field.empty() || field[0].empty() || field[0][0].empty()) 
        {
            return;
        }

        size_t NR = field.size();
        size_t NTH = field[0].size();
        size_t NZ = field[0][0].size();

        // Iterate over the vertical levels.
        for (size_t k = 0; k < NZ; ++k) 
        {
            std::vector<std::vector<double>> horizontal_slice(NR, std::vector<double>(NTH));

            // Iterate over the horizontal levels.
            for (size_t i = 0; i < NR; ++i) 
            {
                // Iterate over the horizontal levels.
                for (size_t j = 0; j < NTH; ++j) 
                {
                    horizontal_slice[i][j] = static_cast<double>(field[i][j][k]);
                }
            }

            // Apply 2D filter
            apply_2d(horizontal_slice, dx, dy);

            // Iterate over the horizontal levels to put back the filtered field.
            for (size_t i = 0; i < NR; ++i) 
            {
                // Iterate over the horizontal levels to put back the filtered field.
                for (size_t j = 0; j < NTH; ++j) 
                {
                    field[i][j][k] = static_cast<float>(horizontal_slice[i][j]);
                }
            }
        }
    }

    /**
     * @brief Get filter name
     */
    virtual std::string name() const = 0;
};

/**
 * @brief Spectral Gaussian correlation filter
 *
 * TODO: Re-implement spectral correlation filtering
 * This is a COMEBACK SECTION - SpectralGaussianFilter class was temporarily removed
 * due to vtable compilation issues during initial integration.
 *
 * Full implementation should:
 * 1. Apply Gaussian correlation in spectral space using FFT
 * 2. Correlation function: exp(-0.5 * (k*L)^2) where L is correlation length
 * 3. More accurate than recursive filtering for large correlation scales
 * 4. Handle 2D FFT operations efficiently
 *
 * Technical issue: vtable generation failed during initial compilation
 * Workaround: Using RecursiveGaussianFilter as default
 */

/**
 * @brief Recursive filter approximation (future implementation)
 *
 * Simpler and faster than spectral filter but approximate.
 * Uses successive passes of 1D filters.
 */
class RecursiveGaussianFilter : public CorrelationFilter 
{
public:
    RecursiveGaussianFilter(double Lx, double Ly);
    void apply_2d(std::vector<std::vector<double>>& field, double dx, double dy) override;
    std::string name() const override { return "recursive_gaussian"; }

private:
    double Lx_, Ly_;
};

//==============================================================================
// Factory function for correlation filters
//==============================================================================

/**
 * @brief Create correlation filter instance
 * @param filter_id Filter type ("spectral_gaussian", "recursive_gaussian")
 * @param Lx Correlation length x
 * @param Ly Correlation length y
 * @return Unique pointer to correlation filter
 */
std::unique_ptr<CorrelationFilter> create_correlation_filter(
    const std::string& filter_id,
    double Lx,
    double Ly
);

} // namespace chaos
