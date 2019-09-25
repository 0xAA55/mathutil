#ifndef _MATHUTIL_CONF_H_
#define _MATHUTIL_CONF_H_ 1

// Option: MATHUTIL_USE_DOUBLE
// This option is used to determine the basic calculating unit (which is 'real_t') of mathutil.
// Normally, 'real_t' will be defined to 'float', and it can be defined to 'double'.
// This option changes API parameter type, 
#ifndef MATHUTIL_USE_DOUBLE
#  define MATHUTIL_USE_DOUBLE 0
#endif
// Option: MATHUTIL_DETECT_CPU
// This option changes the behavior of runtime CPU detection.
// The runtime CPU detection is designed for compatibility, but it may causes the complexity of data storeing,
// result in low efficiency.
// Define this option to 0 disables the runtime CPU detection, and enables compile time intrinsics selection.
#ifndef MATHUTIL_DETECT_CPU
#  define MATHUTIL_DETECT_CPU 0
#endif

#if defined(_MSC_VER)
#  define ALIGNED_(x) __declspec(align(x))
#  ifndef MATHUTIL_VAR_NOT_ALIGNED
#    define MATHUTIL_VAR_NOT_ALIGNED 0
#  endif
#  if _MSC_VER < 1910 || defined(_DEBUG)
#    ifndef MATHUTIL_REFONLY
#      define MATHUTIL_REFONLY 0
#    endif
#  endif
#  if _MSC_VER >= 1910
#    ifndef MATHUTIL_REFONLY
#      define MATHUTIL_REFONLY 0
#    endif
#  endif
#  if defined(_M_AMD64)
#    define __amd64__ 1
#  endif
#elif defined(__GNUC__) || defined(__clang__)
#  define ALIGNED_(x) __attribute__ ((aligned(x)))
#  ifndef MATHUTIL_REFONLY
#    define MATHUTIL_REFONLY 0
#  endif
#  ifndef MATHUTIL_VAR_NOT_ALIGNED
#    define MATHUTIL_VAR_NOT_ALIGNED 0
#  endif
#else
#  define ALIGNED_(x)
#  if !MATHUTIL_DETECT_CPU
#    ifndef MATHUTIL_REFONLY
#      define MATHUTIL_REFONLY 1
#    endif
#  endif
#  ifndef MATHUTIL_VAR_NOT_ALIGNED
#    define MATHUTIL_VAR_NOT_ALIGNED 1
#  endif
#endif
#ifndef MATHUTIL_VAR_ASSUME_ALIGNED
#  define MATHUTIL_VAR_ASSUME_ALIGNED 0
#endif // !MATHUTIL_VAR_ASSUME_ALIGNED
#if MATHUTIL_DETECT_CPU
#ifndef MATHUTIL_ALLOW_TESTED_SLOW_IMPLEMENTS
#  define MATHUTIL_ALLOW_TESTED_SLOW_IMPLEMENTS 0
#endif
#endif // !MATHUTIL_DETECT_CPU

#if defined(__SSE2__) || defined(__x86_64__) || defined(__amd64__)
#define HAVE_SSE2 1
#endif

#if !MATHUTIL_DETECT_CPU && !MATHUTIL_REFONLY
#  if defined(__SSSE3__)
#    define HAVE_SSSE3 1
#  endif
#  if defined(__SSE4_1__)
#    define HAVE_SSE41 1
#  endif
#  if defined(__AVX__)
#    define HAVE_AVX 1
#  endif
#  if defined(__XOP__)
#    define HAVE_XOP 1
#  endif
#endif // !MATHUTIL_DETECT_CPU

#ifdef HAVE_AVX2
#  ifndef HAVE_AVX
#    define HAVE_AVX 1
#  endif
#endif
#ifdef HAVE_XOP
#  ifndef HAVE_AVX
#    define HAVE_AVX 1
#  endif
#endif
#ifdef HAVE_AVX
#  ifndef HAVE_SSE41
#    define HAVE_SSE41 1
#  endif
#endif
#ifdef HAVE_SSE41
#  ifndef HAVE_SSSE3
#    define HAVE_SSSE3 1
#  endif
#endif
#ifdef HAVE_SSSE3
#  define HAVE_SSE2 1
#endif
#ifdef HAVE_SSE2
#  define HAVE_SSE 1
#endif

#include <stddef.h>
#include <stdint.h>

#if MATHUTIL_DETECT_CPU
#define math_func(r,n,arg) math_extern r(*n)arg
#endif // !MATHUTIL_DETECT_CPU

#include "mathutil_funclist.h"

#if defined(math_func)
#undef math_func
#endif // !defined(math_func)

#endif _MATHUTIL_CONF_H_
