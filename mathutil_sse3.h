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

// mathutil_sse3
// Author: 0xAA55
// The implemention of mathutil API which uses SSE2 intrinsics.

#include"mathutil_sse_common.h"

#if HAVE_SSE3 || __INTELLISENSE__

#ifndef sse3_func
  #if __INTELLISENSE__
    #define sse3_func(r,n) r n
    #include<pmmintrin.h>
  #endif // __INTELLISENSE__
#endif

#if !MATHUTIL_USE_DOUBLE

#if !vec4_abs_implemented || MATHUTIL_DETECT_CPU
sse3_func(vec4_t, vec4_abs)(vec4_t v)
{
	vec4_t r;
	vec4_set_iresult(&r, _mm_and_si128(vec4_loadi(v), _mm_set1_epi32(0x7fffffff)));
	return r;
}
#define vec4_abs_implemented 1
#endif

#if !vec4_sgn_implemented || MATHUTIL_DETECT_CPU
sse3_func(vec4_t, vec4_sgn)(vec4_t v)
{
	vec4_t r;
	vec4_set_iresult(&r, _mm_or_si128(_mm_and_si128(vec4_loadi(v), _mm_set1_epi32(0x80000000)), _mm_set1_epi32(0x3F800000)));
	return r;
}
#define vec4_sgn_implemented 1
#endif

#if !vec4_invert_implemented || MATHUTIL_DETECT_CPU
sse3_func(vec4_t, vec4_invert)(vec4_t v)
{
	vec4_t r;
	vec4_set_iresult(&r, _mm_xor_si128(vec4_loadi(v), _mm_set1_epi32(0x80000000)));
	return r;
}
#define vec4_invert_implemented 1
#endif

#if !vec4_dot_implemented || MATHUTIL_DETECT_CPU
sse3_func(real_t, vec4_dot)(vec4_t v1, vec4_t v2)
{
	__m128 r1 = _mm_mul_ps(_mm_load_ps_a(&v1.x), _mm_load_ps_a(&v2.x));
	__m128 r2 = _mm_hadd_ps(r1, r1);
	real_t r;
	_mm_store_ss(&r, _mm_hadd_ps(r2, r2));
	_mm_sfence();
	_compiler_barrier;
	return r;
}
#define vec4_dot_implemented 1
#endif

#if !vec4_length_implemented || MATHUTIL_DETECT_CPU
sse3_func(real_t, vec4_length)(vec4_t v)
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
#define vec4_length_implemented 1
#endif

#if !vec4_normalize_implemented || MATHUTIL_DETECT_CPU
sse3_func(vec4_t, vec4_normalize)(vec4_t v)
{
	vec4_t r;
	__m128 mv = _mm_load_ps_a(&v.x);
	__m128 m = _mm_mul_ps(mv, mv);
	__m128 t = _mm_hadd_ps(m, m);
	_mm_store_ps_a(&r.x, _mm_mul_ps(mv, _mm_rsqrt_ps(_mm_hadd_ps(t, t))));
	return r;
}
#define vec4_normalize_implemented 1
#endif

#else // MATHUTIL_USE_DOUBLE

#if !vec4_abs_implemented || MATHUTIL_DETECT_CPU
sse3_func(vec4_t, vec4_abs)(vec4_t v)
{
	vec4_t r;
	const ALIGNED_(16) int64_t ins[2] = { 0x7fffffffffffffffll, 0x7fffffffffffffffll };
	__m128i mns = _mm_load_si128((__m128i*) & ins);
	vec4_set_iresult_xy(&r, _mm_and_si128(vec4_loadi_xy(v), mns));
	vec4_set_iresult_zw(&r, _mm_and_si128(vec4_loadi_zw(v), mns));
	return r;
}
#define vec4_abs_implemented 1
#endif

#if !vec4_sgn_implemented || MATHUTIL_DETECT_CPU
sse3_func(vec4_t, vec4_sgn)(vec4_t v)
{
	vec4_t r;
	const ALIGNED_(16) int64_t isgnb[2] = { 0x8000000000000000ll, 0x8000000000000000ll };
	const ALIGNED_(16) int64_t i1d[2] = { 0x3ff0000000000000ll, 0x3ff0000000000000ll };
	__m128i msgnb = _mm_load_si128((__m128i*) & isgnb);
	__m128i m1d = _mm_load_si128((__m128i*) & i1d);
	vec4_set_iresult_xy(&r, _mm_or_si128(_mm_and_si128(vec4_loadi_xy(v), msgnb), m1d));
	vec4_set_iresult_zw(&r, _mm_or_si128(_mm_and_si128(vec4_loadi_zw(v), msgnb), m1d));
	return r;
}
#define vec4_sgn_implemented 1
#endif

#if !vec4_invert_implemented || MATHUTIL_DETECT_CPU
sse3_func(vec4_t, vec4_invert)(vec4_t v)
{
	vec4_t r;
	const ALIGNED_(16) int64_t isgnb[2] = { 0x8000000000000000ll, 0x8000000000000000ll };
	__m128i msgnb = _mm_load_si128((__m128i*) & isgnb);
	vec4_set_iresult_xy(&r, _mm_xor_si128(vec4_loadi_xy(v), msgnb));
	vec4_set_iresult_zw(&r, _mm_xor_si128(vec4_loadi_zw(v), msgnb));
	return r;
}
#define vec4_invert_implemented 1
#endif

#if !vec4_dot_implemented || MATHUTIL_DETECT_CPU
sse3_func(real_t, vec4_dot)(vec4_t v1, vec4_t v2)
{
	__m128d mxaddz_yaddw = _mm_add_pd(
		_mm_mul_pd(vec4_load_xy(v1), vec4_load_xy(v2)),
		_mm_mul_pd(vec4_load_zw(v1), vec4_load_zw(v2)));
	real_t r;
	_mm_store_sd(&r, _mm_hadd_pd(mxaddz_yaddw, mxaddz_yaddw));
	_mm_sfence();
	_compiler_barrier;
	return r;
}
#define vec4_dot_implemented 1
#endif

#if !vec4_length_implemented || MATHUTIL_DETECT_CPU
sse3_func(real_t, vec4_length)(vec4_t v)
{
	__m128d mvxy = vec4_load_xy(v);
	__m128d mvzw = vec4_load_zw(v);
	__m128d mxaddz_yaddw = _mm_add_pd(_mm_mul_pd(mvxy, mvxy), _mm_mul_pd(mvzw, mvzw));
	__m128d lsq = _mm_hadd_pd(mxaddz_yaddw, mxaddz_yaddw);

	real_t r;
	_mm_store_sd(&r, _mm_sqrt_pd(lsq));
	_mm_sfence();
	_compiler_barrier;
	return r;
}
#define vec4_length_implemented 1
#endif

#if !vec4_normalize_implemented || MATHUTIL_DETECT_CPU
sse3_func(vec4_t, vec4_normalize)(vec4_t v)
{
	vec4_t r;
	__m128d mvxy = vec4_load_xy(v);
	__m128d mvzw = vec4_load_zw(v);
	__m128d mxaddz_yaddw = _mm_add_pd(_mm_mul_pd(mvxy, mvxy), _mm_mul_pd(mvzw, mvzw));
	__m128d lsq = _mm_hadd_pd(mxaddz_yaddw, mxaddz_yaddw);
	__m128d length = _mm_sqrt_pd(lsq);
	vec4_set_result_xy(&r, _mm_div_pd(mvxy, length));
	vec4_set_result_zw(&r, _mm_div_pd(mvzw, length));
	return r;
}
#define vec4_normalize_implemented 1
#endif

#endif // MATHUTIL_USE_DOUBLE

#endif
