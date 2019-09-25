#ifndef _MATHUTIL_REF_H_
#define _MATHUTIL_REF_H_

#include"mathutil_conf.h"
#include"mathutil.h"

#if MATHUTIL_REFONLY
#define math_func(r,n,arg) r n arg
#else
#define math_func(r,n,arg) r n ## _ref arg
#endif // !MATHUTIL_DETECT_CPU
#include"mathutil_funclist.h"
#undef math_func

#endif // _MATHUTIL_BASE_H_
