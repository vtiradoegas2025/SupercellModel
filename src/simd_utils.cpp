#include "simd_utils.hpp"

// Include SIMD intrinsics based on available instruction sets and architecture
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    // x86/x64 architecture
    #ifdef SIMD_AVX512_AVAILABLE
        #include <immintrin.h>
    #elif defined(SIMD_AVX_AVAILABLE)
        #include <immintrin.h>
    #elif defined(SIMD_SSE_AVAILABLE)
        #include <xmmintrin.h>  // SSE
        #include <emmintrin.h>  // SSE2
        #ifdef __SSE4_1__
            #include <smmintrin.h>  // SSE4.1
        #endif
    #endif
#elif defined(__aarch64__) || defined(__arm__)
    // ARM architecture - use NEON
    #include <arm_neon.h>
    #undef SIMD_WIDTH
    #define SIMD_WIDTH 4  // 4 floats (128 bits) for NEON
#else
    // Unknown architecture - scalar fallback
    #undef SIMD_WIDTH
    #define SIMD_WIDTH 1
#endif

namespace simd_utils {

// Runtime detection (simplified - in production would use CPUID)
SIMDType get_available_simd() {
#ifdef SIMD_AVX512_AVAILABLE
    return SIMDType::AVX512;
#elif defined(SIMD_AVX_AVAILABLE)
    return SIMDType::AVX;
#elif defined(SIMD_SSE_AVAILABLE)
    return SIMDType::SSE;
#else
    return SIMDType::NONE;
#endif
}

// Scalar fallback implementations
namespace scalar {
    void add_vectors(
        const float* a,
        const float* b,
        float* result,
        int count
    ) {
        for (int i = 0; i < count; ++i) {
            result[i] = a[i] + b[i];
        }
    }
    
    void multiply_vectors(
        const float* a,
        const float* b,
        float* result,
        int count
    ) {
        for (int i = 0; i < count; ++i) {
            result[i] = a[i] * b[i];
        }
    }
    
    void fma_vectors(
        const float* a,
        const float* b,
        const float* c,
        float* result,
        int count
    ) {
        for (int i = 0; i < count; ++i) {
            result[i] = a[i] * b[i] + c[i];
        }
    }
}

#ifdef SIMD_AVX_AVAILABLE
// AVX-256 implementations (8 floats)
void add_vectors(
    const float* a,
    const float* b,
    float* result,
    int count
) {
    int i = 0;
    int simd_end = (count / 8) * 8;
    
    // Vectorized section
    for (; i < simd_end; i += 8) {
        __m256 va = _mm256_loadu_ps(&a[i]);
        __m256 vb = _mm256_loadu_ps(&b[i]);
        __m256 vresult = _mm256_add_ps(va, vb);
        _mm256_storeu_ps(&result[i], vresult);
    }
    
    // Scalar remainder
    for (; i < count; ++i) {
        result[i] = a[i] + b[i];
    }
}

void multiply_vectors(
    const float* a,
    const float* b,
    float* result,
    int count
) {
    int i = 0;
    int simd_end = (count / 8) * 8;
    
    // Vectorized section
    for (; i < simd_end; i += 8) {
        __m256 va = _mm256_loadu_ps(&a[i]);
        __m256 vb = _mm256_loadu_ps(&b[i]);
        __m256 vresult = _mm256_mul_ps(va, vb);
        _mm256_storeu_ps(&result[i], vresult);
    }
    
    // Scalar remainder
    for (; i < count; ++i) {
        result[i] = a[i] * b[i];
    }
}

void fma_vectors(
    const float* a,
    const float* b,
    const float* c,
    float* result,
    int count
) {
    int i = 0;
    int simd_end = (count / 8) * 8;
    
    // Vectorized section (using FMA if available)
    #ifdef __FMA__
    for (; i < simd_end; i += 8) {
        __m256 va = _mm256_loadu_ps(&a[i]);
        __m256 vb = _mm256_loadu_ps(&b[i]);
        __m256 vc = _mm256_loadu_ps(&c[i]);
        __m256 vresult = _mm256_fmadd_ps(va, vb, vc);
        _mm256_storeu_ps(&result[i], vresult);
    }
    #else
    // Fallback: multiply then add
    for (; i < simd_end; i += 8) {
        __m256 va = _mm256_loadu_ps(&a[i]);
        __m256 vb = _mm256_loadu_ps(&b[i]);
        __m256 vc = _mm256_loadu_ps(&c[i]);
        __m256 vmul = _mm256_mul_ps(va, vb);
        __m256 vresult = _mm256_add_ps(vmul, vc);
        _mm256_storeu_ps(&result[i], vresult);
    }
    #endif
    
    // Scalar remainder
    for (; i < count; ++i) {
        result[i] = a[i] * b[i] + c[i];
    }
}

#elif defined(SIMD_SSE_AVAILABLE)
// SSE-128 implementations (4 floats)
void add_vectors(
    const float* a,
    const float* b,
    float* result,
    int count
) {
    int i = 0;
    int simd_end = (count / 4) * 4;
    
    // Vectorized section
    for (; i < simd_end; i += 4) {
        __m128 va = _mm_loadu_ps(&a[i]);
        __m128 vb = _mm_loadu_ps(&b[i]);
        __m128 vresult = _mm_add_ps(va, vb);
        _mm_storeu_ps(&result[i], vresult);
    }
    
    // Scalar remainder
    for (; i < count; ++i) {
        result[i] = a[i] + b[i];
    }
}

void multiply_vectors(
    const float* a,
    const float* b,
    float* result,
    int count
) {
    int i = 0;
    int simd_end = (count / 4) * 4;
    
    // Vectorized section
    for (; i < simd_end; i += 4) {
        __m128 va = _mm_loadu_ps(&a[i]);
        __m128 vb = _mm_loadu_ps(&b[i]);
        __m128 vresult = _mm_mul_ps(va, vb);
        _mm_storeu_ps(&result[i], vresult);
    }
    
    // Scalar remainder
    for (; i < count; ++i) {
        result[i] = a[i] * b[i];
    }
}

void fma_vectors(
    const float* a,
    const float* b,
    const float* c,
    float* result,
    int count
) {
    int i = 0;
    int simd_end = (count / 4) * 4;
    
    // Vectorized section
    for (; i < simd_end; i += 4) {
        __m128 va = _mm_loadu_ps(&a[i]);
        __m128 vb = _mm_loadu_ps(&b[i]);
        __m128 vc = _mm_loadu_ps(&c[i]);
        __m128 vmul = _mm_mul_ps(va, vb);
        __m128 vresult = _mm_add_ps(vmul, vc);
        _mm_storeu_ps(&result[i], vresult);
    }
    
    // Scalar remainder
    for (; i < count; ++i) {
        result[i] = a[i] * b[i] + c[i];
    }
}

#else
// Scalar fallback when no SIMD available
void add_vectors(
    const float* a,
    const float* b,
    float* result,
    int count
) {
    scalar::add_vectors(a, b, result, count);
}

void multiply_vectors(
    const float* a,
    const float* b,
    float* result,
    int count
) {
    scalar::multiply_vectors(a, b, result, count);
}

void fma_vectors(
    const float* a,
    const float* b,
    const float* c,
    float* result,
    int count
) {
    scalar::fma_vectors(a, b, c, result, count);
}
#endif

// Process 8 floats (wrapper for AVX)
void process_8_floats(
    const float* a,
    const float* b,
    const float* c,
    float* result,
    int count
) {
    fma_vectors(a, b, c, result, count);
}

// Process 4 floats (wrapper for SSE)
void process_4_floats(
    const float* a,
    const float* b,
    const float* c,
    float* result,
    int count
) {
    fma_vectors(a, b, c, result, count);
}

} // namespace simd_utils
