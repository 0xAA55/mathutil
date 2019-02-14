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

#if MATHUTIL_NOT_ALIGNED && !MATHUTIL_ASSUME_ALIGNED
#define _mm_load_ps_a _mm_loadu_ps
#define _mm_load_pd_a _mm_loadu_pd
#define _mm_store_ps_a _mm_storeu_ps
#define _mm_store_pd_a _mm_storeu_pd
#define _mm_load_si128_a _mm_lddqu_si128
#define _mm_store_si128_a _mm_storeu_si128
#else // !MATHUTIL_NOT_ALIGNED || MATHUTIL_ASSUME_ALIGNED
#define _mm_load_ps_a _mm_load_ps
#define _mm_load_pd_a _mm_load_pd
#define _mm_store_ps_a _mm_store_ps
#define _mm_store_pd_a _mm_store_pd
#define _mm_load_si128_a _mm_load_si128
#define _mm_store_si128_a _mm_store_si128
#endif

#define mathsimd_func(r,n) r n

#if !MATHUTIL_USE_DOUBLE


#if MATHUTIL_DETECT_SIMD
int mathutil_ssse3_implements()
{
	if(!CPUID_SSSE3()) return 0;

	// Boo

	return 1;
}
#endif // MATHUTIL_DETECT_SIMD

#else // MATHUTIL_USE_DOUBLE

#if MATHUTIL_DETECT_SIMD
int mathutil_ssse3_implements()
{
	if(!CPUID_SSSE3()) return 0;
	
	// Boo

	return 1;
}
#endif // MATHUTIL_DETECT_SIMD

#endif // !MATHUTIL_USE_DOUBLE
