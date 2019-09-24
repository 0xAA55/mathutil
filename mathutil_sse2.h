#ifndef _MATHUTIL_SSE2_H_
#define _MATHUTIL_SSE2_H_

#include"mathutil_conf.h"
#include"mathutil.h"

#define math_func(r,n,arg) r n ## _sse2 arg

#include"mathutil_funclist.h"

#undef math_func

#if MATHUTIL_DETECT_CPU
int mathutil_sse2_implements();
#endif // !MATHUTIL_DETECT_CPU

#endif // _MATHUTIL_MMX_H_
