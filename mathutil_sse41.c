#include"mathutil_sse41.h"
#include"cpudetect.h"
#include<smmintrin.h>

#if _MSC_VER
#if !__INTELLISENSE__
#define _compiler_barrier _ReadWriteBarrier()
#else
#define _compiler_barrier
#endif
#elif defined(__GNUC__) || defined(clang)
#define _compiler_barrier asm volatile("" ::: "memory")
#endif

#if MATHUTIL_VAR_NOT_ALIGNED && !MATHUTIL_VAR_ASSUME_ALIGNED
#define _mm_load_ps_a _mm_loadu_ps
#define _mm_load_pd_a _mm_loadu_pd
#define _mm_store_ps_a _mm_storeu_ps
#define _mm_store_pd_a _mm_storeu_pd
#define _mm_load_si128_a _mm_lddqu_si128
#define _mm_store_si128_a _mm_storeu_si128
#else // !MATHUTIL_VAR_NOT_ALIGNED || MATHUTIL_VAR_ASSUME_ALIGNED
#define _mm_load_ps_a _mm_load_ps
#define _mm_load_pd_a _mm_load_pd
#define _mm_store_ps_a _mm_store_ps
#define _mm_store_pd_a _mm_store_pd
#define _mm_load_si128_a _mm_load_si128
#define _mm_store_si128_a _mm_store_si128
#endif

#define mathsimd_func(r,n) r n

#if !MATHUTIL_USE_DOUBLE

mathsimd_func(real_t,vec4_dot_sse41)(vec4_t v1, vec4_t v2)
{
	real_t r;
	_mm_store_ss(&r, _mm_dp_ps(_mm_load_ps_a(&v1.x), _mm_load_ps_a(&v2.x), 0xFF));
	_mm_sfence();
	_compiler_barrier;
	return r;
}

mathsimd_func(real_t,vec4_length_sse41)(vec4_t v)
{
	__m128 mv = _mm_load_ps_a(&v.x);
	real_t r;
	_mm_store_ss(&r, _mm_sqrt_ss(_mm_dp_ps(mv, mv, 0xFF)));
	_mm_sfence();
	_compiler_barrier;
	return r;
}

mathsimd_func(vec4_t,vec4_normalize_sse41)(vec4_t v)
{
	vec4_t r;
	__m128 mv = _mm_load_ps_a(&v.x);
	_mm_store_ps_a(&r.x, _mm_mul_ps(mv, _mm_rsqrt_ps(_mm_dp_ps(mv, mv, 0xFF))));
	return r;
}

mathsimd_func(vec4_t,vec4_min_sse41)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	_mm_store_si128_a((__m128i*)&r.x, _mm_min_epi32(_mm_load_si128_a((__m128i*)&v1.x), _mm_load_si128_a((__m128i*)&v2.x)));
	return r;
}

mathsimd_func(vec4_t,vec4_max_sse41)(vec4_t v1, vec4_t v2)
{
	vec4_t r;\
	_mm_store_si128_a((__m128i*)&r.x, _mm_max_epi32(_mm_load_si128_a((__m128i*)&v1.x), _mm_load_si128_a((__m128i*)&v2.x)));
	return r;
}

#if MATHUTIL_DETECT_CPU
int mathutil_sse41_implements()
{
	if(!CPUID_SSE41()) return 0;
	
	vec4_dot = vec4_dot_sse41;
	vec4_min = vec4_min_sse41;
	vec4_max = vec4_max_sse41;
	vec4_length = vec4_length_sse41;
	vec4_normalize = vec4_normalize_sse41;

	return 1;
}
#endif // MATHUTIL_DETECT_CPU

#else // MATHUTIL_USE_DOUBLE

#if MATHUTIL_DETECT_CPU
int mathutil_sse41_implements()
{
	if(!CPUID_SSE41()) return 0;
	

	return 1;
}
#endif // MATHUTIL_DETECT_CPU

#endif // !MATHUTIL_USE_DOUBLE
