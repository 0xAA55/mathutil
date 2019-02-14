#ifndef _MATHUTIL_SSSE3_H_
#define _MATHUTIL_SSSE3_H_

#include"mathutil.h"

#define math_func(r_hint,r,n,arg,carg) r n ## _ssse3 arg

#include"mathutil_funclist.h"

#undef math_func

#if MATHUTIL_DETECT_SIMD
int mathutil_ssse3_implements();
#endif // !MATHUTIL_DETECT_SIMD

#endif // _MATHUTIL_MMX_H_
