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

#include"mathutil_conf.h"

#include<stdint.h>
#include<math.h>
#include<float.h>

#if MATHUTIL_USE_DOUBLE
typedef double real_t, * real_p;
#define REAL_DEFINED 1

#define real_max DBL_MAX
#define real_dig DBL_DIG
#define real_mant_dig DBL_MANT_DIG
#define real_max_exp DBL_MAX_EXP
#define real_min_exp DBL_MIN_EXP
#define real_max_10_exp DBL_MAX_10_EXP
#define real_min_10_exp DBL_MIN_10_EXP

#else // !MATHUTIL_USE_DOUBLE
typedef float real_t, * real_p;
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
}vec4_t, * vec4_p;
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
	__m128d m_xy;
	__m128d m_zw;
}vec4_t, * vec4_p;
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
}vec4_t, * vec4_p;
#define VEC4_DEFINED 1
#define VEC4_WITH_M256_XYZW 1
#endif // !defined(VEC4_DEFINED) && defined(HAVE_AVX)
#endif // !MATHUTIL_DETECT_CPU && !MATHUTIL_REFONLY

#if !defined(VEC4_DEFINED)
typedef struct vec4_struct
{
	real_t x, y, z, w;
}vec4_t, * vec4_p;
#define VEC4_DEFINED 1
#endif // !defined(VEC4_DEFINED)

typedef vec4_t quat_t, * quat_p;

typedef struct mat4_struct
{
	vec4_t x, y, z, w;
}mat4_t, * mat4_p;

#define r_pi (real_t)3.141592653589793238462643383279502884
#define r_1 (real_t)1
#define r_epsilon (real_t)0.000001

#if __cplusplus
#define math_extern extern "C"
#else
#define math_extern extern
#endif

#if !defined(MATHUTIL_INTERNAL)

#define math_func(r, n) math_extern r n

// Functions for scalar numbers
math_func(real_t, r_rnd)(uint32_t* p_seed);
math_func(real_t, r_sin)(real_t x);
math_func(real_t, r_cos)(real_t x);
math_func(real_t, r_tan)(real_t x);
math_func(real_t, r_abs)(real_t x);
math_func(real_t, r_sgn)(real_t x);
math_func(real_t, r_sqr)(real_t x);
math_func(real_t, r_floor)(real_t x);
math_func(real_t, r_ceil)(real_t x);
math_func(real_t, r_atan)(real_t x);
math_func(real_t, r_exp)(real_t x);
math_func(real_t, r_log)(real_t x);
math_func(real_t, r_pow)(real_t x, real_t y);
math_func(real_t, r_mod)(real_t x, real_t y);
math_func(real_t, r_max)(real_t a, real_t b);
math_func(real_t, r_min)(real_t a, real_t b);
math_func(real_t, r_atan2)(real_t y, real_t x);
math_func(real_t, r_clamp)(real_t n, real_t min_, real_t max_);
math_func(real_t, r_lerp)(real_t a, real_t b, real_t s);
math_func(real_t, r_hermite)(real_t s);
math_func(real_t, r_slerp)(real_t a, real_t b, real_t s);

// Functions for vectors
math_func(vec4_t, vec4_t_ctor)(real_t x, real_t y, real_t z, real_t w);
math_func(vec4_t, vec4_flushcomp)(vec4_t v);
math_func(vec4_t, vec4_abs)(vec4_t v);
math_func(vec4_t, vec4_sgn)(vec4_t v);
math_func(vec4_t, vec4_invert)(vec4_t v);
math_func(real_t, vec4_length)(vec4_t v);
math_func(vec4_t, vec4_normalize)(vec4_t v);
math_func(vec4_t, vec4_scale)(vec4_t v, real_t s);
math_func(vec4_t, vec4_pow)(vec4_t v, real_t n);
math_func(real_t, vec4_dot)(vec4_t v1, vec4_t v2);
math_func(vec4_t, vec4_add)(vec4_t v1, vec4_t v2);
math_func(vec4_t, vec4_sub)(vec4_t v1, vec4_t v2);
math_func(vec4_t, vec4_mul)(vec4_t v1, vec4_t v2);
math_func(vec4_t, vec4_div)(vec4_t v1, vec4_t v2);
math_func(vec4_t, vec4_min)(vec4_t v1, vec4_t v2);
math_func(vec4_t, vec4_max)(vec4_t v1, vec4_t v2);
math_func(vec4_t, vec4_cross3)(vec4_t v1, vec4_t v2);
math_func(vec4_t, vec4_clamp)(vec4_t v, real_t min_, real_t max_);
math_func(vec4_t, vec4_rot_quat)(vec4_t v, quat_t q);
math_func(vec4_t, vec4_mul_mat4)(vec4_t v, mat4_t m);
math_func(vec4_t, vec4_mul_mat4_transpose)(vec4_t v, mat4_t m);
math_func(vec4_t, vec4_lerp)(vec4_t v1, vec4_t v2, real_t s);
math_func(vec4_t, vec4_slerp)(vec4_t v1, vec4_t v2, real_t s);

// Functions for quaternions, which can be used to describe rotation
math_func(quat_t, quat_t_ctor)(real_t x, real_t y, real_t z, real_t w);
math_func(quat_t, quat_flushcomp)(quat_t q);
math_func(quat_t, quat_rot_axis)(vec4_t axis, real_t angle);
math_func(quat_t, quat_mul)(quat_t q1, quat_t q2);
math_func(quat_t, quat_add_vec)(quat_t q, vec4_t v, real_t s);

// Functions for matrices
math_func(mat4_t, mat4_t_ctor)(vec4_t mx, vec4_t my, vec4_t mz, vec4_t mw);
math_func(mat4_t, mat4_flushcomp)(mat4_t m);
math_func(mat4_t, mat4_rot_x)(real_t angle);
math_func(mat4_t, mat4_rot_y)(real_t angle);
math_func(mat4_t, mat4_rot_z)(real_t angle);
math_func(mat4_t, mat4_rot_axis)(vec4_t axis, real_t angle);
math_func(mat4_t, mat4_rot_euler)(real_t yaw, real_t pitch, real_t roll);
math_func(mat4_t, mat4_from_quat)(quat_t q);
math_func(mat4_t, mat4_from_quat_transpose)(quat_t q);
math_func(mat4_t, mat4_transpose)(mat4_t m);
math_func(real_t, mat4_det)(mat4_t m);
math_func(int, mat4_inverse)(mat4_t m, real_p det_out, mat4_p mat_out);
math_func(mat4_t, mat4_add)(mat4_t l, mat4_t r);
math_func(mat4_t, mat4_add_s)(mat4_t m, real_t s);
math_func(mat4_t, mat4_add_transpose)(mat4_t l, mat4_t r);
math_func(mat4_t, mat4_sub)(mat4_t l, mat4_t r);
math_func(mat4_t, mat4_sub_s)(mat4_t m, real_t s);
math_func(mat4_t, mat4_sub_transpose)(mat4_t l, mat4_t r);
math_func(mat4_t, mat4_mul)(mat4_t l, mat4_t r);
math_func(mat4_t, mat4_mul_s)(mat4_t m, real_t s);
math_func(mat4_t, mat4_mul_transpose)(mat4_t l, mat4_t r);

#define vec4 vec4_t_ctor
#define quat quat_t_ctor
#define mat4 mat4_t_ctor

#endif

#endif // !_MATHUTIL_H_
