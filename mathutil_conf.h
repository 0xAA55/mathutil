#ifndef _MATHUTIL_CONF_H_
#define _MATHUTIL_CONF_H_ 1

#if defined(_MSC_VER)
#  define ALIGNED_(x) __declspec(align(x))
#  ifndef MATHUTIL_VAR_NOT_ALIGNED
#    define MATHUTIL_VAR_NOT_ALIGNED 0
#  endif
#  if _MSC_VER < 1910 || defined(_DEBUG)
#    ifndef MATHUTIL_DETECT_CPU
#      define MATHUTIL_DETECT_CPU 1
#    endif
#    ifndef MATHUTIL_REFONLY
#      define MATHUTIL_REFONLY 0
#    endif
#  endif
#  if _MSC_VER >= 1910
#    ifndef MATHUTIL_DETECT_CPU
#      define MATHUTIL_DETECT_CPU 0
#    endif
#    ifndef MATHUTIL_REFONLY
#      define MATHUTIL_REFONLY 0
#    endif
#  endif
#elif defined(__GNUC__) || defined(__clang__)
#  define ALIGNED_(x) __attribute__ ((aligned(x)))
#  ifndef MATHUTIL_DETECT_CPU
#    define MATHUTIL_DETECT_CPU 0
#  endif
#  ifndef MATHUTIL_REFONLY
#    define MATHUTIL_REFONLY 0
#  endif
#  ifndef MATHUTIL_VAR_NOT_ALIGNED
#    define MATHUTIL_VAR_NOT_ALIGNED 0
#  endif
#else
#  define ALIGNED_(x)
#  ifndef MATHUTIL_DETECT_CPU
#    define MATHUTIL_DETECT_CPU 0
#  endif
#  ifndef MATHUTIL_REFONLY
#    define MATHUTIL_REFONLY 1
#  endif
#  ifndef MATHUTIL_VAR_NOT_ALIGNED
#    define MATHUTIL_VAR_NOT_ALIGNED 1
#  endif
#endif
#ifndef MATHUTIL_VAR_ASSUME_ALIGNED
#  define MATHUTIL_VAR_ASSUME_ALIGNED 0
#endif
#ifndef MATHUTIL_ALLOW_TESTED_SLOW_IMPLEMENTS
#  define MATHUTIL_ALLOW_TESTED_SLOW_IMPLEMENTS 1
#endif

#if MATHUTIL_DETECT_CPU
#define math_func(r_hint,r,n,arg,carg) math_extern r(*n) arg
#else // !MATHUTIL_DETECT_CPU
#define math_func(r_hint,r,n,arg,carg) r n arg
#endif // !MATHUTIL_DETECT_CPU

#include"mathutil_funclist.h"

#undef math_func

#endif _MATHUTIL_CONF_H_
