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

#include"mathutil.h"
#include"mathutil_ref.h"
#include"mathutil_sse.h"
#include"mathutil_sse2.h"
#include"mathutil_sse3.h"
#include"mathutil_ssse3.h"
#include"mathutil_sse41.h"

#if !MATHUTIL_REFONLY
#if MATHUTIL_DETECT_CPU
void mathutil_ref_implements()
{
#	define math_func(r,n,arg) n = n ## _ref
#	include"mathutil_funclist.h"
#	undef math_func
}

static int g_mathutil_initialized = 0;

static void mathutil_init()
{
	if(g_mathutil_initialized) return;

	mathutil_ref_implements();
	mathutil_sse_implements();
	mathutil_sse2_implements();
	mathutil_sse3_implements();
	mathutil_ssse3_implements();
	mathutil_sse41_implements();

	g_mathutil_initialized = 1;
}

#	define math_func(r,n,arg) r n ## _first arg{mathutil_init(); return n carg;} r(*n)arg = n ## _first
#	include"mathutil_funclist.h"
#	undef math_func

#else // !MATHUTIL_DETECT_CPU

#if HAVE_AVX2

#if COMPILER_FLAVOR == 2

#define math_func(r,n,arg) r n ##  arg __attribute__ ((weak, alias ("" # n # "_avx2")));
#include"mathutil_funclist.h"
#undef math_func

#else // COMPILER_FLAVOR == else

#define math_func(r,n,arg, carg) r n ##  arg { return n ## _avx2 carg;}
#include"mathutil_funclist_withcarg.h"
#undef math_func

#endif // COMPILER_FLAVOR

#elif HAVE_XOP

#if COMPILER_FLAVOR == 2

#define math_func(r,n,arg) r n ##  arg __attribute__ ((weak, alias ("" # n # "_xop")));
#include"mathutil_funclist.h"
#undef math_func

#else // COMPILER_FLAVOR == else

#define math_func(r,n,arg, carg) r n ##  arg { return n ## _xop carg;}
#include"mathutil_funclist_withcarg.h"
#undef math_func

#endif // COMPILER_FLAVOR

#elif HAVE_AVX

#if COMPILER_FLAVOR == 2

#define math_func(r,n,arg) r n ##  arg __attribute__ ((weak, alias ("" # n # "_avx")));
#include"mathutil_funclist.h"
#undef math_func

#else // COMPILER_FLAVOR == else

#define math_func(r,n,arg, carg) r n ##  arg { return n ## _avx carg;}
#include"mathutil_funclist_withcarg.h"
#undef math_func

#endif // COMPILER_FLAVOR

#elif HAVE_SSE41

#if COMPILER_FLAVOR == 2

#define math_func(r,n,arg) r n ##  arg __attribute__ ((weak, alias ("" # n # "_sse41")));
#include"mathutil_funclist.h"
#undef math_func

#else // COMPILER_FLAVOR == else

#define math_func(r,n,arg, carg) r n ##  arg { return n ## _sse41 carg;}
#include"mathutil_funclist_withcarg.h"
#undef math_func

#endif // COMPILER_FLAVOR

#elif HAVE_SSSE3

#if COMPILER_FLAVOR == 2

#define math_func(r,n,arg) r n ##  arg __attribute__ ((weak, alias ("" # n # "_ssse3")));
#include"mathutil_funclist.h"
#undef math_func

#else // COMPILER_FLAVOR == else

#define math_func(r,n,arg, carg) r n ##  arg { return n ## _ssse3 carg;}
#include"mathutil_funclist_withcarg.h"
#undef math_func

#endif // COMPILER_FLAVOR

#elif HAVE_SSE3

#if COMPILER_FLAVOR == 2

#define math_func(r,n,arg) r n ##  arg __attribute__ ((weak, alias ("" # n # "_sse3")));
#include"mathutil_funclist.h"
#undef math_func

#else // COMPILER_FLAVOR == else

#define math_func(r,n,arg, carg) r n ##  arg { return n ## _sse3 carg;}
#include"mathutil_funclist_withcarg.h"
#undef math_func

#endif // COMPILER_FLAVOR

#elif HAVE_SSE2

#if COMPILER_FLAVOR == 2

#define math_func(r,n,arg) r n ##  arg __attribute__ ((weak, alias ("" # n # "_sse2")));
#include"mathutil_funclist.h"
#undef math_func

#else // COMPILER_FLAVOR == else

#define math_func(r,n,arg, carg) r n ##  arg { return n ## _sse2 carg;}
#include"mathutil_funclist_withcarg.h"
#undef math_func

#endif // COMPILER_FLAVOR

#elif HAVE_SSE

#if COMPILER_FLAVOR == 2

#define math_func(r,n,arg) r n ##  arg __attribute__ ((weak, alias ("" # n # "_sse")));
#include"mathutil_funclist.h"
#undef math_func

#else // COMPILER_FLAVOR == else

#define math_func(r,n,arg, carg) r n ##  arg { return n ## _sse carg;}
#include"mathutil_funclist_withcarg.h"
#undef math_func

#endif // COMPILER_FLAVOR

#endif

#endif // !MATHUTIL_DETECT_CPU
#endif // !MATHUTIL_REFONLY