#ifndef _MATHUTIL_H_
#define _MATHUTIL_H_

#include<stdint.h>
#include<math.h>
#include<float.h>

#if MATHUTIL_USE_DOUBLE
typedef double real_t, *real_p;

#define real_max DBL_MAX
#define real_dig DBL_DIG
#define real_mant_dig DBL_MANT_DIG
#define real_max_exp DBL_MAX_EXP
#define real_min_exp DBL_MIN_EXP
#define real_max_10_exp DBL_MAX_10_EXP
#define real_min_10_exp DBL_MIN_10_EXP

#else // !MATHUTIL_USE_DOUBLE
typedef float real_t, *real_p;

#define real_max FLT_MAX
#define real_dig FLT_DIG
#define real_mant_dig FLT_MANT_DIG
#define real_max_exp FLT_MAX_EXP
#define real_min_exp FLT_MIN_EXP
#define real_max_10_exp FLT_MAX_10_EXP
#define real_min_10_exp FLT_MIN_10_EXP

#endif // !MATHUTIL_USE_DOUBLE

#if defined(_MSC_VER)
#define ALIGNED_(x) __declspec(align(x))
#elif defined(__GNUC__)
#define ALIGNED_(x) __attribute__ ((aligned(x)))
#else
#define ALIGNED_(x)
#endif
#define MATHUTIL_NOT_ALIGNED 1
#define MATHUTIL_ASSUME_ALIGNED 0
#define MATHUTIL_ALLOW_TESTED_SLOW_IMPLEMENTS 1

typedef struct vec4_struct
{
	real_t x, y, z, w;
}vec4_t, *vec4_p;

typedef vec4_t quat_t, quat_p;

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

#if MATHUTIL_DETECT_SIMD
#define math_func(r_hint,r,n,arg,carg) math_extern r(*n) arg
#else // !MATHUTIL_DETECT_SIMD
#define math_func(r_hint,r,n,arg,carg) r n arg
#endif

#include"mathutil_funclist.h"

#undef math_func

#endif // !_MATHUTIL_H_
