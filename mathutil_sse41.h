#ifndef _MATHUTIL_SSE41_H_
#define _MATHUTIL_SSE41_H_

#include"mathutil_conf.h"
#include"mathutil.h"

#define math_func(r_hint,r,n,arg,carg) r n ## _sse41 arg

#include"mathutil_funclist.h"

#undef math_func

#if MATHUTIL_DETECT_CPU
int mathutil_sse41_implements();
#endif // !MATHUTIL_DETECT_CPU

#endif // _MATHUTIL_MMX_H_
