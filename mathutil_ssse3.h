#ifndef _MATHUTIL_SSSE3_H_
#define _MATHUTIL_SSSE3_H_

#include"mathutil_conf.h"
#include"mathutil.h"

#define math_func(r_hint,r,n,arg,carg) r n ## _ssse3 arg

#include"mathutil_funclist.h"

#undef math_func

#if MATHUTIL_DETECT_CPU
int mathutil_ssse3_implements();
#endif // !MATHUTIL_DETECT_CPU

#endif // _MATHUTIL_MMX_H_
