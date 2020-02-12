// MIT License
// 
// Copyright (c) 2019 0xaa55
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

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
#  define COMPILER_FLAVOR 1
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
#  pragma warning(disable:4201)
#  pragma warning(disable:4239)
#elif defined(__GNUC__) || defined(__clang__)
#  define COMPILER_FLAVOR 2
#  define ALIGNED_(x) __attribute__ ((aligned(x)))
#  ifndef MATHUTIL_REFONLY
#    define MATHUTIL_REFONLY 0
#  endif
#  ifndef MATHUTIL_VAR_NOT_ALIGNED
#    define MATHUTIL_VAR_NOT_ALIGNED 0
#  endif
#else
#  define COMPILER_FLAVOR 0
#  define ALIGNED_(x)
#  if !MATHUTIL_DETECT_CPU
#    ifndef MATHUTIL_REFONLY
#      define MATHUTIL_REFONLY 0
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
#endif // !!MATHUTIL_DETECT_CPU && !MATHUTIL_REFONLY

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
#  define HAVE_SSE3 1
#endif
#ifdef HAVE_SSE3
#  define HAVE_SSE2 1
#endif
#ifdef HAVE_SSE2
#  define HAVE_SSE 1
#endif

#if !defined(MATHUTIL_REFONLY)
#define MATHUTIL_REFONLY (!(HAVE_SSE || HAVE_SSE2 || HAVE_SSE3 || HAVE_SSSE3 || HAVE_SSE41 || HAVE_AVX || HAVE_XOP || HAVE_AVX2))
#endif // !defined(MATHUTIL_REFONLY)

#ifdef __cplusplus
  #include <cstddef>
  #include <cstdint>
#else
  #include <stddef.h>
  #include <stdint.h>
#endif

#endif // _MATHUTIL_CONF_H_
