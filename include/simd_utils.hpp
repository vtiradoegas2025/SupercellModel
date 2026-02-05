#ifndef SIMD_UTILS_HPP
#define SIMD_UTILS_HPP

#include <cstdint>
#include <cstring>

/* SIMD utility functions for vectorized operations
 * Supports SSE, AVX, and AVX-512 instruction sets
 * Falls back to scalar operations if SIMD unavailable
 */

namespace simd_utils {

// Detect available SIMD instruction sets at compile time
#if defined(__AVX512F__)
    #define SIMD_AVX512_AVAILABLE
    #define SIMD_WIDTH 16  // 16 floats (512 bits)
#elif defined(__AVX__)
    #define SIMD_AVX_AVAILABLE
    #define SIMD_WIDTH 8   // 8 floats (256 bits)
#elif defined(__SSE4_1__) || defined(__SSE__)
    #define SIMD_SSE_AVAILABLE
    #define SIMD_WIDTH 4   // 4 floats (128 bits)
#else
    #define SIMD_WIDTH 1   // Scalar fallback
#endif

// Runtime detection of SIMD capabilities
enum class SIMDType {
    NONE = 0,
    SSE = 1,
    AVX = 2,
    AVX512 = 3
};

// Get available SIMD type at runtime
SIMDType get_available_simd();

// Process 8 floats at once (AVX-256)
// Element-wise operations: result[i] = a[i] * b[i] + c[i]
void process_8_floats(
    const float* a,
    const float* b,
    const float* c,
    float* result,
    int count
);

// Process 4 floats at once (SSE-128)
void process_4_floats(
    const float* a,
    const float* b,
    const float* c,
    float* result,
    int count
);

// Vectorized addition: result = a + b
void add_vectors(
    const float* a,
    const float* b,
    float* result,
    int count
);

// Vectorized multiplication: result = a * b
void multiply_vectors(
    const float* a,
    const float* b,
    float* result,
    int count
);

// Vectorized fused multiply-add: result = a * b + c
void fma_vectors(
    const float* a,
    const float* b,
    const float* c,
    float* result,
    int count
);

// Scalar fallback implementations
namespace scalar {
    void add_vectors(
        const float* a,
        const float* b,
        float* result,
        int count
    );
    
    void multiply_vectors(
        const float* a,
        const float* b,
        float* result,
        int count
    );
    
    void fma_vectors(
        const float* a,
        const float* b,
        const float* c,
        float* result,
        int count
    );
}

} // namespace simd_utils

#endif // SIMD_UTILS_HPP
