#ifndef SIMD_UTILS_HPP
#define SIMD_UTILS_HPP

#include <cstdint>
#include <cstring>

/**
 * @file simd_utils.hpp
 * @brief SIMD capability detection and vector math utility declarations.
 *
 * Declares architecture-aware vector kernels and scalar fallbacks.
 * Compile-time macros expose the available SIMD width for callers.
 * Implementations choose optimized instructions when available.
 */

namespace simd_utils
{

#if defined(__AVX512F__)
#define SIMD_AVX512_AVAILABLE
#define SIMD_WIDTH 16
#elif defined(__AVX__)
#define SIMD_AVX_AVAILABLE
#define SIMD_WIDTH 8
#elif defined(__SSE4_1__) || defined(__SSE__)
#define SIMD_SSE_AVAILABLE
#define SIMD_WIDTH 4
#else
#define SIMD_WIDTH 1
#endif

enum class SIMDType
{
    NONE = 0,
    SSE = 1,
    AVX = 2,
    AVX512 = 3
};

/**
 * @brief Detects the SIMD ISA supported by the current build/runtime path.
 * @return Selected SIMD type.
 */
SIMDType get_available_simd();

/**
 * @brief Processes vectors using AVX-width lanes where possible.
 */
void process_8_floats(const float* a, const float* b, const float* c, float* result, int count);

/**
 * @brief Processes vectors using SSE-width lanes where possible.
 */
void process_4_floats(const float* a, const float* b, const float* c, float* result, int count);

/**
 * @brief Computes element-wise vector addition.
 */
void add_vectors(const float* a, const float* b, float* result, int count);

/**
 * @brief Computes element-wise vector multiplication.
 */
void multiply_vectors(const float* a, const float* b, float* result, int count);

/**
 * @brief Computes element-wise fused multiply-add.
 */
void fma_vectors(const float* a, const float* b, const float* c, float* result, int count);

namespace scalar
{
/**
 * @brief Scalar fallback for element-wise vector addition.
 */
void add_vectors(const float* a, const float* b, float* result, int count);

/**
 * @brief Scalar fallback for element-wise vector multiplication.
 */
void multiply_vectors(const float* a, const float* b, float* result, int count);

/**
 * @brief Scalar fallback for element-wise fused multiply-add.
 */
void fma_vectors(const float* a, const float* b, const float* c, float* result, int count);
} // namespace scalar

} // namespace simd_utils

#endif // SIMD_UTILS_HPP
