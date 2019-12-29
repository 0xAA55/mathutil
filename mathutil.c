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

#if !MATHUTIL_REFONLY
#if MATHUTIL_DETECT_CPU
void mathutil_ref_implements()
{
#	define math_func(r,n,arg,carg) n = n ## _ref
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

#	define math_func(r,n,arg,carg) r n ## _first arg{mathutil_init(); return n carg;} r(*n)arg = n ## _first
#	include"mathutil_funclist.h"
#	undef math_func

#else // !MATHUTIL_DETECT_CPU

#define sse41_func(r, n) r n
#include"mathutil_sse41.h"

#define ssse3_func(r, n) r n
#include"mathutil_ssse3.h"

#define sse3_func(r, n) r n
#include"mathutil_sse3.h"

#define sse2_func(r, n) r n
#include"mathutil_sse2.h"

#define sse_func(r, n) r n
#include"mathutil_sse.h"

#define ref_func(r, n) r n
#include"mathutil_ref.h"

#endif // !MATHUTIL_DETECT_CPU
#endif // MATHUTIL_REFONLY