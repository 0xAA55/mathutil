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

#include"mathutil_ssse3.h"
#include"cpudetect.h"
#include<tmmintrin.h>

#if _MSC_VER
#if !__INTELLISENSE__
#define _compiler_barrier _ReadWriteBarrier()
#else
#define _compiler_barrier
#endif
#elif defined(__GNUC__) || defined(clang)
#define _compiler_barrier asm volatile("" ::: "memory")
#endif

#if MATHUTIL_VAR_NOT_ALIGNED && !MATHUTIL_VAR_ASSUME_ALIGNED
#define _mm_load_ps_a _mm_loadu_ps
#define _mm_load_pd_a _mm_loadu_pd
#define _mm_store_ps_a _mm_storeu_ps
#define _mm_store_pd_a _mm_storeu_pd
#define _mm_load_si128_a _mm_lddqu_si128
#define _mm_store_si128_a _mm_storeu_si128
#else // !MATHUTIL_VAR_NOT_ALIGNED || MATHUTIL_VAR_ASSUME_ALIGNED
#define _mm_load_ps_a _mm_load_ps
#define _mm_load_pd_a _mm_load_pd
#define _mm_store_ps_a _mm_store_ps
#define _mm_store_pd_a _mm_store_pd
#define _mm_load_si128_a _mm_load_si128
#define _mm_store_si128_a _mm_store_si128
#endif

#define mathsimd_func(r,n) r n

#if !MATHUTIL_USE_DOUBLE


#if MATHUTIL_DETECT_CPU
int mathutil_ssse3_implements()
{
	if(!CPUID_SSSE3()) return 0;

	// Boo

	return 1;
}
#endif // MATHUTIL_DETECT_CPU

#else // MATHUTIL_USE_DOUBLE

#if MATHUTIL_DETECT_CPU
int mathutil_ssse3_implements()
{
	if(!CPUID_SSSE3()) return 0;
	
	// Boo

	return 1;
}
#endif // MATHUTIL_DETECT_CPU

#endif // !MATHUTIL_USE_DOUBLE
