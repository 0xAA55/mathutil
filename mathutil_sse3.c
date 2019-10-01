// MIT License
// 
// Copyright (c) 2019 0xaa55
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include"mathutil_sse3.h"
#include"cpudetect.h"
#include<pmmintrin.h>

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

mathsimd_func(vec4_t,vec4_abs_sse3)(vec4_t v)
{
	vec4_t r;
	_mm_store_si128_a((__m128i*)&r.x, _mm_and_si128(_mm_load_si128_a((__m128i*)&v.x), _mm_set1_epi32(0x7fffffff)));
	return r;
}

mathsimd_func(vec4_t,vec4_sgn_sse3)(vec4_t v)
{
	vec4_t r;
	_mm_store_si128_a((__m128i*)&r.x, _mm_or_si128(_mm_and_si128(_mm_load_si128_a((__m128i*)&v.x), _mm_set1_epi32(0x80000000)), _mm_set1_epi32(0x3F800000)));
	return r;
}

mathsimd_func(vec4_t,vec4_invert_sse3)(vec4_t v)
{
	vec4_t r;
	_mm_store_si128_a((__m128i*)&r.x, _mm_xor_si128(_mm_load_si128_a((__m128i*)&v.x), _mm_set1_epi32(0x80000000)));
	return r;
}

mathsimd_func(real_t,vec4_dot_sse3)(vec4_t v1, vec4_t v2)
{
	__m128 r1 = _mm_mul_ps(_mm_load_ps_a(&v1.x), _mm_load_ps_a(&v2.x));
	__m128 r2 = _mm_hadd_ps(r1, r1);
	real_t r;
	_mm_store_ss(&r, _mm_hadd_ps(r2, r2));
	_mm_sfence();
	_compiler_barrier;
	return r;
}

mathsimd_func(real_t,vec4_length_sse3)(vec4_t v)
{
	__m128 mv = _mm_load_ps_a(&v.x);
	__m128 m = _mm_mul_ps(mv, mv);
	__m128 t = _mm_hadd_ps(m, m);
	__m128 length = _mm_sqrt_ss(_mm_hadd_ps(t, t));
	real_t r;
	_mm_store_ss(&r, length);
	_mm_sfence();
	_compiler_barrier;
	return r;
}

mathsimd_func(vec4_t,vec4_normalize_sse3)(vec4_t v)
{
	vec4_t r;
	__m128 mv = _mm_load_ps_a(&v.x);
	__m128 m = _mm_mul_ps(mv, mv);
	__m128 t = _mm_hadd_ps(m, m);
	_mm_store_ps_a(&r.x, _mm_mul_ps(mv, _mm_rsqrt_ps(_mm_hadd_ps(t, t))));
	return r;
}

#if MATHUTIL_DETECT_CPU
int mathutil_sse3_implements()
{
	if(!CPUID_SSE3()) return 0;
	
	vec4_abs = vec4_abs_sse3;
	vec4_sgn = vec4_sgn_sse3;
	vec4_dot = vec4_dot_sse3;
	vec4_invert = vec4_invert_sse3;
	vec4_length = vec4_length_sse3;
	vec4_normalize = vec4_normalize_sse3;

	return 1;
}
#endif // MATHUTIL_DETECT_CPU

#else // MATHUTIL_USE_DOUBLE

mathsimd_func(vec4_t,vec4_abs_sse3)(vec4_t v)
{
	vec4_t r;
	const ALIGNED_(16) int64_t ins[2] = {0x7fffffffffffffffll, 0x7fffffffffffffffll};
	__m128i mns = _mm_load_si128((__m128i*)&ins);
	_mm_store_si128_a((__m128i*)&r.x, _mm_and_si128(_mm_load_si128_a((__m128i*)&v.x), mns));
	_mm_store_si128_a((__m128i*)&r.z, _mm_and_si128(_mm_load_si128_a((__m128i*)&v.z), mns));
	return r;
}

mathsimd_func(vec4_t,vec4_sgn_sse3)(vec4_t v)
{
	vec4_t r;
	const ALIGNED_(16) int64_t isgnb[2] = {0x8000000000000000ll, 0x8000000000000000ll};
	const ALIGNED_(16) int64_t i1d[2] = {0x3ff0000000000000ll, 0x3ff0000000000000ll};
	__m128i msgnb = _mm_load_si128((__m128i*)&isgnb);
	__m128i m1d = _mm_load_si128((__m128i*)&i1d);
	_mm_store_si128_a((__m128i*)&r.x, _mm_or_si128(_mm_and_si128(_mm_load_si128_a((__m128i*)&v.x), msgnb), m1d));
	_mm_store_si128_a((__m128i*)&r.z, _mm_or_si128(_mm_and_si128(_mm_load_si128_a((__m128i*)&v.z), msgnb), m1d));
	return r;
}

mathsimd_func(vec4_t,vec4_invert_sse3)(vec4_t v)
{
	vec4_t r;
	const ALIGNED_(16) int64_t isgnb[2] = {0x8000000000000000ll, 0x8000000000000000ll};
	__m128i msgnb = _mm_load_si128((__m128i*)&isgnb);
	_mm_store_si128_a((__m128i*)&r.x, _mm_xor_si128(_mm_load_si128_a((__m128i*)&v.x), msgnb));
	_mm_store_si128_a((__m128i*)&r.z, _mm_xor_si128(_mm_load_si128_a((__m128i*)&v.z), msgnb));
	return r;
}

mathsimd_func(real_t,vec4_dot_sse3)(vec4_t v1, vec4_t v2)
{
	__m128d mxaddz_yaddw = _mm_add_pd(
		_mm_mul_pd(_mm_load_pd_a(&v1.x), _mm_load_pd_a(&v2.x)),
		_mm_mul_pd(_mm_load_pd_a(&v1.z), _mm_load_pd_a(&v2.z)));
	real_t r;
	_mm_store_sd(&r, _mm_hadd_pd(mxaddz_yaddw, mxaddz_yaddw));
	_mm_sfence();
	_compiler_barrier;
	return r;
}

mathsimd_func(real_t,vec4_length_sse3)(vec4_t v)
{
	__m128d mvxy = _mm_load_pd_a(&v.x);
	__m128d mvzw = _mm_load_pd_a(&v.z);
	__m128d mxaddz_yaddw = _mm_add_pd(_mm_mul_pd(mvxy, mvxy), _mm_mul_pd(mvzw, mvzw));
	__m128d lsq = _mm_hadd_pd(mxaddz_yaddw, mxaddz_yaddw);

	real_t r;
	_mm_store_sd(&r, _mm_sqrt_pd(lsq));
	_mm_sfence();
	_compiler_barrier;
	return r;
}

mathsimd_func(vec4_t,vec4_normalize_sse3)(vec4_t v)
{
	vec4_t r;
	__m128d mvxy = _mm_load_pd_a(&v.x);
	__m128d mvzw = _mm_load_pd_a(&v.z);
	__m128d mxaddz_yaddw = _mm_add_pd(_mm_mul_pd(mvxy, mvxy), _mm_mul_pd(mvzw, mvzw));
	__m128d lsq = _mm_hadd_pd(mxaddz_yaddw, mxaddz_yaddw);
	__m128d length = _mm_sqrt_pd(lsq);
	_mm_store_pd_a(&r.x, _mm_div_pd(mvxy, length));
	_mm_store_pd_a(&r.z, _mm_div_pd(mvzw, length));
	return r;
}

#if MATHUTIL_DETECT_CPU
int mathutil_sse3_implements()
{
	if(!CPUID_SSE3()) return 0;
	
	vec4_abs = vec4_abs_sse3;
	vec4_sgn = vec4_sgn_sse3;
	vec4_dot = vec4_dot_sse3;
	vec4_invert = vec4_invert_sse3;
	vec4_length = vec4_length_sse3;
	vec4_normalize = vec4_normalize_sse3;

	return 1;
}
#endif // MATHUTIL_DETECT_CPU

#endif // !MATHUTIL_USE_DOUBLE
