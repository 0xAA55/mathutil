#ifndef _MATHUTIL_SSE3_H_
#define _MATHUTIL_SSE3_H_

#include"mathutil_conf.h"
#include"mathutil.h"

#define math_func(r,n,arg) r n ## _sse3 arg

#include"mathutil_funclist.h"

#undef math_func

#if MATHUTIL_DETECT_CPU
int mathutil_sse3_implements();
#endif // !MATHUTIL_DETECT_CPU

#endif // _MATHUTIL_MMX_H_
