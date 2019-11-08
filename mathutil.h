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

#ifndef _MATHUTIL_H_
#define _MATHUTIL_H_

#include<stdint.h>
#include<math.h>
#include<float.h>

#if MATHUTIL_USE_DOUBLE
typedef double real_t, *real_p;
#define REAL_DEFINED 1

#define real_max DBL_MAX
#define real_dig DBL_DIG
#define real_mant_dig DBL_MANT_DIG
#define real_max_exp DBL_MAX_EXP
#define real_min_exp DBL_MIN_EXP
#define real_max_10_exp DBL_MAX_10_EXP
#define real_min_10_exp DBL_MIN_10_EXP

#else // !MATHUTIL_USE_DOUBLE
typedef float real_t, *real_p;
#define REAL_DEFINED 1

#define real_max FLT_MAX
#define real_dig FLT_DIG
#define real_mant_dig FLT_MANT_DIG
#define real_max_exp FLT_MAX_EXP
#define real_min_exp FLT_MIN_EXP
#define real_max_10_exp FLT_MAX_10_EXP
#define real_min_10_exp FLT_MIN_10_EXP

#endif // !MATHUTIL_USE_DOUBLE

#if !MATHUTIL_DETECT_CPU && !MATHUTIL_REFONLY

#if defined(HAVE_SSE)
#include<xmmintrin.h>
#include<immintrin.h>

#if !MATHUTIL_USE_DOUBLE
typedef union vec4_union
{
	struct
	{
		real_t x, y, z, w;
	};
	__m128 m_xyzw;
}vec4_t, *vec4_p;
#define VEC4_DEFINED 1
#define VEC4_WITH_M128_XYZW 1
#endif // !MATHUTIL_USE_DOUBLE
#endif // !defined(HAVE_SSE)

#if defined(HAVE_SSE2)
#include<emmintrin.h>
#if !defined(VEC4_DEFINED) && !defined(HAVE_AVX)
typedef union vec4_union
{
	struct
	{
		real_t x, y, z, w;
	};
	__m128 m_xy;
	__m128 m_zw;
}vec4_t, *vec4_p;
#define VEC4_DEFINED 1
#define VEC4_WITH_M128_XY_ZW 1
#endif // !defined(VEC4_DEFINED) && !defined(HAVE_AVX)
#endif // !defined(HAVE_SSE2)

#if !defined(VEC4_DEFINED) && defined(HAVE_AVX)
typedef union vec4_union
{
	struct
	{
		real_t x, y, z, w;
	};
	__m256 m_xyzw;
}vec4_t, *vec4_p;
#define VEC4_DEFINED 1
#define VEC4_WITH_M256_XYZW 1
#endif // !defined(VEC4_DEFINED) && defined(HAVE_AVX)
#endif // !MATHUTIL_DETECT_CPU && !MATHUTIL_REFONLY

#if !defined(VEC4_DEFINED)
typedef struct vec4_struct
{
	real_t x, y, z, w;
}vec4_t, *vec4_p;
#define VEC4_DEFINED 1
#endif // !defined(VEC4_DEFINED)

typedef vec4_t quat_t, *quat_p;

typedef struct mat4_struct
{
	vec4_t x, y, z, w;
}mat4_t, *mat4_p;

#define r_pi (real_t)3.141592653589793238462643383279502884
#define r_1 (real_t)1
#define r_epsilon (real_t)0.000001

#ifdef __cplusplus
#define math_extern extern"C"
#else
#define math_extern extern
#endif

#if !defined(MATHUTIL_INTERNAL)

#define math_func(r,n,arg,carg) math_extern r n arg
#include"mathutil_funclist.h"
#undef math_func

#endif

#endif // !_MATHUTIL_H_
