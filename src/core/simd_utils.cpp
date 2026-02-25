/**
 * @file simd_utils.cpp
 * @brief Core runtime implementation for the tornado model.
 *
 * Provides simulation orchestration and subsystem integration
 * for dynamics, numerics, physics, and runtime execution paths.
 * This file belongs to the primary src/core execution layer.
 */

#include "simd_utils.hpp"

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
    #ifdef SIMD_AVX512_AVAILABLE
        #include <immintrin.h>
    #elif defined(SIMD_AVX_AVAILABLE)
        #include <immintrin.h>
    #elif defined(SIMD_SSE_AVAILABLE)
        #include <xmmintrin.h>
        #include <emmintrin.h>
        #ifdef __SSE4_1__
            #include <smmintrin.h>
        #endif
    #endif
#elif defined(__aarch64__) || defined(__arm__)
    #include <arm_neon.h>
    #undef SIMD_WIDTH
    #define SIMD_WIDTH 4
#else
    #undef SIMD_WIDTH
    #define SIMD_WIDTH 1
#endif

namespace simd_utils {

/**
 * @brief Detects the highest SIMD instruction set enabled at compile time.
 */
SIMDType get_available_simd() 
{
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

namespace scalar 
{
    /**
     * @brief Scalar fallback for elementwise vector addition.
     */
    void add_vectors(const float* a, const float* b, float* result, int count) 
    {
        for (int i = 0; i < count; ++i) {result[i] = a[i] + b[i];}
    }
    
    /**
     * @brief Scalar fallback for elementwise vector multiplication.
     */
    void multiply_vectors(const float* a, const float* b, float* result, int count) 
    {
        for (int i = 0; i < count; ++i) {result[i] = a[i] * b[i]; }
    }
    
    /**
     * @brief Scalar fallback for fused multiply-add style operation.
     */
    void fma_vectors(const float* a, const float* b, const float* c, float* result, int count) 
    {
        for (int i = 0; i < count; ++i) {result[i] = a[i] * b[i] + c[i]; }
    }
}

#ifdef SIMD_AVX_AVAILABLE
/**
 * @brief AVX-accelerated elementwise vector addition.
 */
void add_vectors(const float* a,const float* b, float* result, int count) 
{
    int i = 0;
    int simd_end = (count / 8) * 8;
    
    for (; i < simd_end; i += 8) 
    {
        __m256 va = _mm256_loadu_ps(&a[i]);
        __m256 vb = _mm256_loadu_ps(&b[i]);
        __m256 vresult = _mm256_add_ps(va, vb);
        _mm256_storeu_ps(&result[i], vresult);
    }
    
    for (; i < count; ++i) {result[i] = a[i] + b[i];}
}

/**
 * @brief AVX-accelerated elementwise vector multiplication.
 */
void multiply_vectors( const float* a, const float* b, float* result, int count) 
{
    int i = 0;
    int simd_end = (count / 8) * 8;
    
    for (; i < simd_end; i += 8) 
    {
        __m256 va = _mm256_loadu_ps(&a[i]);
        __m256 vb = _mm256_loadu_ps(&b[i]);
        __m256 vresult = _mm256_mul_ps(va, vb);
        _mm256_storeu_ps(&result[i], vresult);
    }
    
    for (; i < count; ++i) { result[i] = a[i] * b[i];}
}

/**
 * @brief AVX-accelerated fused multiply-add style operation.
 */
void fma_vectors(const float* a, const float* b, const float* c, float* result, int count) 
{
    int i = 0;
    int simd_end = (count / 8) * 8;
    
    #ifdef __FMA__
    for (; i < simd_end; i += 8) 
    {
        __m256 va = _mm256_loadu_ps(&a[i]);
        __m256 vb = _mm256_loadu_ps(&b[i]);
        __m256 vc = _mm256_loadu_ps(&c[i]);
        __m256 vresult = _mm256_fmadd_ps(va, vb, vc);
        _mm256_storeu_ps(&result[i], vresult);
    }
    #else
    for (; i < simd_end; i += 8) 
    {
        __m256 va = _mm256_loadu_ps(&a[i]);
        __m256 vb = _mm256_loadu_ps(&b[i]);
        __m256 vc = _mm256_loadu_ps(&c[i]);
        __m256 vmul = _mm256_mul_ps(va, vb);
        __m256 vresult = _mm256_add_ps(vmul, vc);
        _mm256_storeu_ps(&result[i], vresult);
    }
    #endif
    
    for (; i < count; ++i) {result[i] = a[i] * b[i] + c[i];}
}

#elif defined(SIMD_SSE_AVAILABLE)
/**
 * @brief SSE-accelerated elementwise vector addition.
 */
void add_vectors( const float* a, const float* b, float* result, int count) 
{
    int i = 0;
    int simd_end = (count / 4) * 4;
    
    for (; i < simd_end; i += 4) 
    {
        __m128 va = _mm_loadu_ps(&a[i]);
        __m128 vb = _mm_loadu_ps(&b[i]);
        __m128 vresult = _mm_add_ps(va, vb);
        _mm_storeu_ps(&result[i], vresult);
    }
    
    for (; i < count; ++i) { result[i] = a[i] + b[i];}
}

/**
 * @brief SSE-accelerated elementwise vector multiplication.
 */
void multiply_vectors(const float* a, const float* b, float* result, int count) 
{
    int i = 0;
    int simd_end = (count / 4) * 4;
    
    for (; i < simd_end; i += 4) 
    {
        __m128 va = _mm_loadu_ps(&a[i]);
        __m128 vb = _mm_loadu_ps(&b[i]);
        __m128 vresult = _mm_mul_ps(va, vb);
        _mm_storeu_ps(&result[i], vresult);
    }
    
    for (; i < count; ++i) { result[i] = a[i] * b[i];}
}

/**
 * @brief SSE-accelerated fused multiply-add style operation.
 */
void fma_vectors(const float* a, const float* b, const float* c, float* result, int count) 
{
    int i = 0;
    int simd_end = (count / 4) * 4;
    
    for (; i < simd_end; i += 4) 
    {
        __m128 va = _mm_loadu_ps(&a[i]);
        __m128 vb = _mm_loadu_ps(&b[i]);
        __m128 vc = _mm_loadu_ps(&c[i]);
        __m128 vmul = _mm_mul_ps(va, vb);
        __m128 vresult = _mm_add_ps(vmul, vc);
        _mm_storeu_ps(&result[i], vresult);
    }
    
    for (; i < count; ++i) {result[i] = a[i] * b[i] + c[i];}
}

#else
/**
 * @brief Scalar-dispatch elementwise vector addition.
 */
void add_vectors(const float* a, const float* b, float* result, int count) {scalar::add_vectors(a, b, result, count);}

/**
 * @brief Scalar-dispatch elementwise vector multiplication.
 */
void multiply_vectors(const float* a, const float* b, float* result,int count){scalar::multiply_vectors(a, b, result, count);}

/**
 * @brief Scalar-dispatch fused multiply-add style operation.
 */
void fma_vectors( const float* a, const float* b, const float* c, float* result, int count){scalar::fma_vectors(a, b, c, result, count);}
#endif

/**
 * @brief Processes up to eight floats using the active SIMD path.
 */
void process_8_floats(const float* a, const float* b, const float* c,float* result, int count) {fma_vectors(a, b, c, result, count);}

/**
 * @brief Processes up to four floats using the active SIMD path.
 */
void process_4_floats(const float* a, const float* b, const float* c,float* result, int count) {fma_vectors(a, b, c, result, count);}

}
