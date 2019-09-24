#ifndef _MATHUTIL_REF_H_
#define _MATHUTIL_REF_H_

#include"mathutil_conf.h"
#include"mathutil.h"

#if MATHUTIL_DETECT_CPU
#define math_func(r,n,arg) r n ## _ref arg

#include"mathutil_funclist.h"

#undef math_func
#endif // !MATHUTIL_DETECT_CPU

#endif // _MATHUTIL_BASE_H_
