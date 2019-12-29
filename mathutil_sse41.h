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

#include"mathutil_sse_common.h"

#if HAVE_SSE41 || __INTELLISENSE__

#ifndef sse41_func
  #if __INTELLISENSE__
    #define sse41_func(r,n) r n
    #include<smmintrin.h>
  #endif // __INTELLISENSE__
#endif

#if !MATHUTIL_USE_DOUBLE
#if !vec4_dot_implemented || MATHUTIL_DETECT_CPU
sse41_func(real_t, vec4_dot)(vec4_t v1, vec4_t v2)
{
	real_t r;
	_mm_store_ss(&r, _mm_dp_ps(vec4_load(v1), vec4_load(v2), 0xFF));
	_mm_sfence();
	_compiler_barrier;
	return r;
}
#define vec4_dot_implemented 1
#endif

#if !vec4_length_implemented || MATHUTIL_DETECT_CPU
sse41_func(real_t, vec4_length)(vec4_t v)
{
	__m128 mv = vec4_load(v);
	real_t r;
	_mm_store_ss(&r, _mm_sqrt_ss(_mm_dp_ps(mv, mv, 0xFF)));
	_mm_sfence();
	_compiler_barrier;
	return r;
}
#define vec4_length_implemented 1
#endif

#if !vec4_normalize_implemented || MATHUTIL_DETECT_CPU
sse41_func(vec4_t, vec4_normalize)(vec4_t v)
{
	vec4_t r;
	__m128 mv = vec4_load(v);
	vec4_set_result(&r, _mm_mul_ps(mv, _mm_rsqrt_ps(_mm_dp_ps(mv, mv, 0xFF))));
	return r;
}
#define vec4_normalize_implemented 1
#endif

#if !vec4_min_implemented || MATHUTIL_DETECT_CPU
sse41_func(vec4_t, vec4_min)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	vec4_set_iresult(&r, _mm_min_epi32(vec4_loadi(v1), vec4_loadi(v2)));
	return r;
}
#define vec4_min_implemented 1
#endif

#if !vec4_max_implemented || MATHUTIL_DETECT_CPU
sse41_func(vec4_t, vec4_max)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	vec4_set_iresult(&r, _mm_max_epi32(vec4_loadi(v1), vec4_loadi(v2)));
	return r;
}
#define vec4_max_implemented 1
#endif
#else // MATHUTIL_USE_DOUBLE

#endif // MATHUTIL_USE_DOUBLE

#endif