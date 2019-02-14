#ifndef _MATHUTIL_REF_H_
#define _MATHUTIL_REF_H_

#include"mathutil.h"

#if MATHUTIL_DETECT_SIMD
#define math_func(r_hint,r,n,arg,carg) r n ## _ref arg

#include"mathutil_funclist.h"

#undef math_func
#endif // !MATHUTIL_DETECT_SIMD

#endif // _MATHUTIL_BASE_H_
