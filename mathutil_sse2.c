#include"mathutil_sse2.h"
#include"cpudetect.h"
#include<emmintrin.h>

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
#define _mm_load_si128_a _mm_loadu_si128
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

mathsimd_func(vec4_t,vec4_abs_sse2)(vec4_t v)
{
	vec4_t r;
	_mm_store_si128_a((__m128i*)&r.x, _mm_and_si128(_mm_load_si128_a((__m128i*)&v.x), _mm_set1_epi32(0x7fffffff)));
	return r;
}

mathsimd_func(vec4_t,vec4_sgn_sse2)(vec4_t v)
{
	vec4_t r;
	_mm_store_si128_a((__m128i*)&r.x, _mm_or_si128(_mm_and_si128(_mm_load_si128_a((__m128i*)&v.x), _mm_set1_epi32(0x80000000)), _mm_set1_epi32(0x3F800000)));
	return r;
}

mathsimd_func(vec4_t,vec4_invert_sse2)(vec4_t v)
{
	vec4_t r;
	_mm_store_si128_a((__m128i*)&r.x, _mm_xor_si128(_mm_load_si128_a((__m128i*)&v.x), _mm_set1_epi32(0x80000000)));
	return r;
}

#if MATHUTIL_DETECT_SIMD
int mathutil_sse2_implements()
{
	if(!CPUID_SSE2()) return 0;
	
	vec4_abs = vec4_abs_sse2;
	vec4_sgn = vec4_sgn_sse2;
	vec4_invert = vec4_invert_sse2;

	// Boo

	return 1;
}
#endif // MATHUTIL_DETECT_SIMD

#else // MATHUTIL_USE_DOUBLE

mathsimd_func(vec4_t,vec4_abs_sse2)(vec4_t v)
{
	vec4_t r;
	static const ALIGNED_(16) int64_t ins[2] = {0x7fffffffffffffffll, 0x7fffffffffffffffll};
	__m128i mns = _mm_load_si128((__m128i*)&ins);
	_mm_store_si128_a((__m128i*)&r.x, _mm_and_si128(_mm_load_si128_a((__m128i*)&v.x), mns));
	_mm_store_si128_a((__m128i*)&r.z, _mm_and_si128(_mm_load_si128_a((__m128i*)&v.z), mns));
	return r;
}

mathsimd_func(vec4_t,vec4_sgn_sse2)(vec4_t v)
{
	vec4_t r;
	static const ALIGNED_(16) int64_t isgnb[2] = {0x8000000000000000ll, 0x8000000000000000ll};
	static const ALIGNED_(16) int64_t i1d[2] = {0x3ff0000000000000ll, 0x3ff0000000000000ll};
	__m128i msgnb = _mm_load_si128((__m128i*)&isgnb);
	__m128i m1d = _mm_load_si128((__m128i*)&i1d);
	_mm_store_si128_a((__m128i*)&r.x, _mm_or_si128(_mm_and_si128(_mm_load_si128_a((__m128i*)&v.x), msgnb), m1d));
	_mm_store_si128_a((__m128i*)&r.z, _mm_or_si128(_mm_and_si128(_mm_load_si128_a((__m128i*)&v.z), msgnb), m1d));
	return r;
}

mathsimd_func(vec4_t,vec4_invert_sse2)(vec4_t v)
{
	vec4_t r;
	static const ALIGNED_(16) int64_t isgnb[2] = {0x8000000000000000ll, 0x8000000000000000ll};
	__m128i msgnb = _mm_load_si128((__m128i*)&isgnb);
	_mm_store_si128_a((__m128i*)&r.x, _mm_xor_si128(_mm_load_si128_a((__m128i*)&v.x), msgnb));
	_mm_store_si128_a((__m128i*)&r.z, _mm_xor_si128(_mm_load_si128_a((__m128i*)&v.z), msgnb));
	return r;
}

mathsimd_func(real_t,vec4_length_sse2)(vec4_t v)
{
	__m128d mvxy = _mm_load_pd_a(&v.x);
	__m128d mvzw = _mm_load_pd_a(&v.z);
	__m128d mxaddz_yaddw = _mm_add_pd(_mm_mul_pd(mvxy, mvxy), _mm_mul_pd(mvzw, mvzw));
	__m128d lsq = _mm_add_sd(mxaddz_yaddw, _mm_unpackhi_pd(mxaddz_yaddw, mxaddz_yaddw));

	real_t r;
	_mm_store_sd(&r, _mm_sqrt_pd(lsq));
	_mm_sfence();
	_compiler_barrier;
	return r;
}

mathsimd_func(vec4_t,vec4_normalize_sse2)(vec4_t v)
{
	vec4_t r;
	__m128d mvxy = _mm_load_pd_a(&v.x);
	__m128d mvzw = _mm_load_pd_a(&v.z);
	__m128d mxaddz_yaddw = _mm_add_pd(_mm_mul_pd(mvxy, mvxy), _mm_mul_pd(mvzw, mvzw));
	__m128d lsq = _mm_add_sd(mxaddz_yaddw, _mm_unpackhi_pd(mxaddz_yaddw, mxaddz_yaddw));
	__m128d length = _mm_sqrt_pd(lsq);
	_mm_store_pd_a(&r.x, _mm_div_pd(mvxy, length));
	_mm_store_pd_a(&r.z, _mm_div_pd(mvzw, length));
	return r;
}

mathsimd_func(vec4_t,vec4_scale_sse2)(vec4_t v, real_t s)
{
	vec4_t r;
	__m128d ms = _mm_load1_pd(&s);
	_mm_store_pd_a(&r.x, _mm_mul_pd(_mm_load_pd_a(&v.x), ms));
	_mm_store_pd_a(&r.z, _mm_mul_pd(_mm_load_pd_a(&v.z), ms));
	return r;
}

mathsimd_func(real_t,vec4_dot_sse2)(vec4_t v1, vec4_t v2)
{
	__m128d mxaddz_yaddw = _mm_add_pd(
		_mm_mul_pd(_mm_load_pd_a(&v1.x), _mm_load_pd_a(&v2.x)),
		_mm_mul_pd(_mm_load_pd_a(&v1.z), _mm_load_pd_a(&v2.z)));
	real_t r;
	_mm_store_sd(&r, _mm_add_sd(mxaddz_yaddw, _mm_unpackhi_pd(mxaddz_yaddw, mxaddz_yaddw)));
	_mm_sfence();
	_compiler_barrier;
	return r;
}

mathsimd_func(vec4_t,vec4_add_sse2)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	__m128d mxy = _mm_add_pd(_mm_load_pd_a(&v1.x), _mm_load_pd_a(&v2.x));
	__m128d mzw = _mm_add_pd(_mm_load_pd_a(&v1.z), _mm_load_pd_a(&v2.z));
	_mm_store_pd_a(&r.x, mxy);
	_mm_store_pd_a(&r.z, mzw);
	return r;
}

mathsimd_func(vec4_t,vec4_sub_sse2)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	__m128d mxy = _mm_sub_pd(_mm_load_pd_a(&v1.x), _mm_load_pd_a(&v2.x));
	__m128d mzw = _mm_sub_pd(_mm_load_pd_a(&v1.z), _mm_load_pd_a(&v2.z));
	_mm_store_pd_a(&r.x, mxy);
	_mm_store_pd_a(&r.z, mzw);
	return r;
}

mathsimd_func(vec4_t,vec4_mul_sse2)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	__m128d mxy = _mm_mul_pd(_mm_load_pd_a(&v1.x), _mm_load_pd_a(&v2.x));
	__m128d mzw = _mm_mul_pd(_mm_load_pd_a(&v1.z), _mm_load_pd_a(&v2.z));
	_mm_store_pd_a(&r.x, mxy);
	_mm_store_pd_a(&r.z, mzw);
	return r;
}

mathsimd_func(vec4_t,vec4_div_sse2)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	__m128d mxy = _mm_div_pd(_mm_load_pd_a(&v1.x), _mm_load_pd_a(&v2.x));
	__m128d mzw = _mm_div_pd(_mm_load_pd_a(&v1.z), _mm_load_pd_a(&v2.z));
	_mm_store_pd_a(&r.x, mxy);
	_mm_store_pd_a(&r.z, mzw);
	return r;
}

mathsimd_func(vec4_t,vec4_min_sse2)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	_mm_store_pd_a(&r.x, _mm_min_pd(_mm_load_pd_a(&v1.x), _mm_load_pd_a(&v2.x)));
	_mm_store_pd_a(&r.z, _mm_min_pd(_mm_load_pd_a(&v1.z), _mm_load_pd_a(&v2.z)));
	return r;
}

mathsimd_func(vec4_t,vec4_max_sse2)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	_mm_store_pd_a(&r.x, _mm_max_pd(_mm_load_pd_a(&v1.x), _mm_load_pd_a(&v2.x)));
	_mm_store_pd_a(&r.z, _mm_max_pd(_mm_load_pd_a(&v1.z), _mm_load_pd_a(&v2.z)));
	return r;
}

mathsimd_func(vec4_t,vec4_mul_mat4_sse2)(vec4_t v, mat4_t m)
{
	vec4_t r;
	
	__m128d mvxy = _mm_load_pd_a(&v.x);
	__m128d mvzw = _mm_load_pd_a(&v.z);
	
	_mm_store_pd_a(&r.x, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_load_pd_a(&m.x.x)),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_load_pd_a(&m.y.x))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_load_pd_a(&m.z.x)),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_load_pd_a(&m.w.x)))));
	_mm_store_pd_a(&r.z, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_load_pd_a(&m.x.z)),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_load_pd_a(&m.y.z))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_load_pd_a(&m.z.z)),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_load_pd_a(&m.w.z)))));

	return r;
}

mathsimd_func(vec4_t,vec4_mul_mat4_transpose_sse2)(vec4_t v, mat4_t m)
{
	vec4_t r;
	
	__m128d mvxy = _mm_load_pd_a(&v.x);
	__m128d mvzw = _mm_load_pd_a(&v.z);
	
	__m128d mx_xy = _mm_load_pd_a(&m.x.x);
	__m128d my_xy = _mm_load_pd_a(&m.y.x);
	__m128d mz_xy = _mm_load_pd_a(&m.z.x);
	__m128d mw_xy = _mm_load_pd_a(&m.w.x);

	__m128d mx_zw = _mm_load_pd_a(&m.x.z);
	__m128d my_zw = _mm_load_pd_a(&m.y.z);
	__m128d mz_zw = _mm_load_pd_a(&m.z.z);
	__m128d mw_zw = _mm_load_pd_a(&m.w.z);
	
	_mm_store_pd_a(&r.x, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mx_xy, my_xy)),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mx_xy, my_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mx_zw, my_zw)),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mx_zw, my_zw)))));
	_mm_store_pd_a(&r.z, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mz_xy, mw_xy)),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mz_xy, mw_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mz_zw, mw_zw)),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mz_zw, mw_zw)))));

	return r;
}

mathsimd_func(vec4_t,vec4_lerp_sse2)(vec4_t v1, vec4_t v2, real_t s)
{
	__m128d maxy = _mm_load_pd_a(&v1.x);
	__m128d mazw = _mm_load_pd_a(&v1.z);
	__m128d mbxy = _mm_load_pd_a(&v2.x);
	__m128d mbzw = _mm_load_pd_a(&v2.z);
	__m128d ms = _mm_load1_pd(&s);
	vec4_t r;
	
	_mm_store_pd_a(&r.x, _mm_add_pd(maxy, _mm_mul_pd(_mm_sub_pd(mbxy, maxy), ms)));
	_mm_store_pd_a(&r.z, _mm_add_pd(mazw, _mm_mul_pd(_mm_sub_pd(mbzw, mazw), ms)));
	return r;
}

mathsimd_func(vec4_t,vec4_slerp_sse2)(vec4_t v1, vec4_t v2, real_t s)
{
	__m128d maxy = _mm_load_pd_a(&v1.x);
	__m128d mazw = _mm_load_pd_a(&v1.z);
	__m128d mbxy = _mm_load_pd_a(&v2.x);
	__m128d mbzw = _mm_load_pd_a(&v2.z);
	__m128d ms = _mm_load1_pd(&s);
	vec4_t r;
	
	ms = _mm_mul_pd(_mm_mul_pd(ms, ms), _mm_sub_pd(_mm_set1_pd(3), _mm_mul_pd(_mm_set1_pd(2), ms)));
	ms = _mm_max_pd(_mm_setzero_pd(), _mm_min_pd(_mm_set1_pd(1), ms));
	_mm_store_pd_a(&r.x, _mm_add_pd(maxy, _mm_mul_pd(_mm_sub_pd(mbxy, maxy), ms)));
	_mm_store_pd_a(&r.z, _mm_add_pd(mazw, _mm_mul_pd(_mm_sub_pd(mbzw, mazw), ms)));
	return r;
}

mathsimd_func(quat_t,quat_mul_sse2)(quat_t q1, quat_t q2)
{
	__m128d mq1xy = _mm_load_pd_a(&q1.x);
	__m128d mq1zw = _mm_load_pd_a(&q1.z);
	__m128d mq2xy = _mm_load_pd_a(&q2.x);
	__m128d mq2zw = _mm_load_pd_a(&q2.z);
	__m128d mpn = _mm_set_pd(-1, 1);
	__m128d mnp = _mm_set_pd( 1,-1);
	__m128d mqxy;
	__m128d mqzw;
	quat_t r;
	
	mqxy = _mm_mul_pd(_mm_unpackhi_pd(mq1zw, mq1zw), mq2xy);
	mqxy = _mm_add_pd(mqxy, _mm_mul_pd(_mm_unpackhi_pd(mq1xy, mq1xy), mq2zw));
	mqxy = _mm_add_pd(mqxy, _mm_mul_pd(_mm_mul_pd(mpn, _mm_unpacklo_pd(mq1xy, mq1xy)), _mm_shuffle_pd(mq2zw, mq2zw, _MM_SHUFFLE2(0,1))));
	mqxy = _mm_add_pd(mqxy, _mm_mul_pd(_mm_mul_pd(mnp, _mm_unpacklo_pd(mq1zw, mq1zw)), _mm_shuffle_pd(mq2xy, mq2xy, _MM_SHUFFLE2(0,1))));
	
	mqzw = _mm_mul_pd(_mm_unpackhi_pd(mq1zw, mq1zw), mq2zw);
	mqzw = _mm_sub_pd(mqzw, _mm_mul_pd(_mm_unpackhi_pd(mq1xy, mq1xy), mq2xy));
	mqzw = _mm_add_pd(mqzw, _mm_mul_pd(_mm_mul_pd(mpn, _mm_unpacklo_pd(mq1xy, mq1xy)), _mm_shuffle_pd(mq2xy, mq2xy, _MM_SHUFFLE2(0,1))));
	mqzw = _mm_add_pd(mqzw, _mm_mul_pd(_mm_mul_pd(mpn, _mm_unpacklo_pd(mq1zw, mq1zw)), _mm_shuffle_pd(mq2zw, mq2zw, _MM_SHUFFLE2(0,1))));

	_mm_store_pd_a(&r.x, mqxy);
	_mm_store_pd_a(&r.z, mqzw);
	return r;
}

mathsimd_func(quat_t,quat_add_vec_sse2)(quat_t q, vec4_t v, real_t s)
{
	__m128d mq1xy = _mm_load_pd_a(&v.x);
	__m128d mq1zw = _mm_load_sd(&v.z);
	__m128d mq2xy = _mm_load_pd_a(&q.x);
	__m128d mq2zw = _mm_load_pd_a(&q.z);
	__m128d mpn = _mm_set_pd(-1, 1);
	__m128d mnp = _mm_set_pd( 1,-1);
	__m128d mqxy;
	__m128d mqzw;
	__m128d ms = _mm_load1_pd(&s);
	__m128d mhf = _mm_set1_pd(0.5);
	quat_t r;

	mq1xy = _mm_mul_pd(mq1xy, ms);
	mq1zw = _mm_mul_sd(mq1zw, ms);

	mqxy = _mm_mul_pd(_mm_unpackhi_pd(mq1xy, mq1xy), mq2zw);
	mqxy = _mm_add_pd(mqxy, _mm_mul_pd(_mm_mul_pd(mpn, _mm_unpacklo_pd(mq1xy, mq1xy)), _mm_shuffle_pd(mq2zw, mq2zw, _MM_SHUFFLE2(0,1))));
	mqxy = _mm_add_pd(mqxy, _mm_mul_pd(_mm_mul_pd(mnp, _mm_unpacklo_pd(mq1zw, mq1zw)), _mm_shuffle_pd(mq2xy, mq2xy, _MM_SHUFFLE2(0,1))));

	mqzw = _mm_mul_pd(_mm_mul_pd(mpn, _mm_unpacklo_pd(mq1xy, mq1xy)), _mm_shuffle_pd(mq2xy, mq2xy, _MM_SHUFFLE2(0,1)));
	mqzw = _mm_add_pd(mqzw, _mm_mul_pd(_mm_mul_pd(mpn, _mm_unpacklo_pd(mq1zw, mq1zw)), _mm_shuffle_pd(mq2zw, mq2zw, _MM_SHUFFLE2(0,1))));
	mqzw = _mm_sub_pd(mqzw, _mm_mul_pd(_mm_unpackhi_pd(mq1xy, mq1xy), mq2xy));

	_mm_store_pd_a(&r.x, _mm_add_pd(_mm_mul_pd(mqxy, mhf), mq2xy));
	_mm_store_pd_a(&r.z, _mm_add_pd(_mm_mul_pd(mqzw, mhf), mq2zw));
	return r;
}

mathsimd_func(mat4_t,mat4_transpose_sse2)(mat4_t m)
{
	__m128d mx_xy;
	__m128d my_xy;
	__m128d mz_xy;
	__m128d mw_xy;

	__m128d mx_zw;
	__m128d my_zw;
	__m128d mz_zw;
	__m128d mw_zw;

	__m128d t;
	mat4_t r;
	
	mx_xy = _mm_load_pd_a(&m.x.x);
	my_xy = _mm_load_pd_a(&m.y.x);
	mz_xy = _mm_load_pd_a(&m.z.x);
	mw_xy = _mm_load_pd_a(&m.w.x);

	t = _mm_unpacklo_pd(mx_xy, my_xy);
	_mm_store_pd_a(&r.x.x, t);
	
	t = _mm_unpacklo_pd(mz_xy, mw_xy);
	_mm_store_pd_a(&r.x.z, t);
	
	t = _mm_unpackhi_pd(mx_xy, my_xy);
	_mm_store_pd_a(&r.y.x, t);
	
	t = _mm_unpackhi_pd(mz_xy, mw_xy);
	_mm_store_pd_a(&r.y.z, t);
	
	mx_zw = _mm_load_pd_a(&m.x.z);
	my_zw = _mm_load_pd_a(&m.y.z);
	mz_zw = _mm_load_pd_a(&m.z.z);
	mw_zw = _mm_load_pd_a(&m.w.z);

	t = _mm_unpacklo_pd(mx_zw, my_zw);
	_mm_store_pd_a(&r.z.x, t);
	
	t = _mm_unpacklo_pd(mz_zw, mw_zw);
	_mm_store_pd_a(&r.z.z, t);
	
	t = _mm_unpackhi_pd(mx_zw, my_zw);
	_mm_store_pd_a(&r.w.x, t);
	
	t = _mm_unpackhi_pd(mz_zw, mw_zw);
	_mm_store_pd_a(&r.w.z, t);

	return r;
}

mathsimd_func(mat4_t,mat4_add_sse2)(mat4_t l, mat4_t r)
{
	mat4_t o;
	_mm_store_pd_a(&o.x.x, _mm_add_pd(_mm_load_pd_a(&l.x.x), _mm_load_pd_a(&r.x.x)));
	_mm_store_pd_a(&o.y.x, _mm_add_pd(_mm_load_pd_a(&l.y.x), _mm_load_pd_a(&r.y.x)));
	_mm_store_pd_a(&o.z.x, _mm_add_pd(_mm_load_pd_a(&l.z.x), _mm_load_pd_a(&r.z.x)));
	_mm_store_pd_a(&o.w.x, _mm_add_pd(_mm_load_pd_a(&l.w.x), _mm_load_pd_a(&r.w.x)));
	_mm_store_pd_a(&o.x.z, _mm_add_pd(_mm_load_pd_a(&l.x.z), _mm_load_pd_a(&r.x.z)));
	_mm_store_pd_a(&o.y.z, _mm_add_pd(_mm_load_pd_a(&l.y.z), _mm_load_pd_a(&r.y.z)));
	_mm_store_pd_a(&o.z.z, _mm_add_pd(_mm_load_pd_a(&l.z.z), _mm_load_pd_a(&r.z.z)));
	_mm_store_pd_a(&o.w.z, _mm_add_pd(_mm_load_pd_a(&l.w.z), _mm_load_pd_a(&r.w.z)));
	return o;
}

mathsimd_func(mat4_t,mat4_add_s_sse2)(mat4_t m, real_t s)
{
	__m128d ms = _mm_load1_pd(&s);
	mat4_t r;
	_mm_store_pd_a(&r.x.x, _mm_add_pd(_mm_load_pd_a(&m.x.x), ms));
	_mm_store_pd_a(&r.y.x, _mm_add_pd(_mm_load_pd_a(&m.y.x), ms));
	_mm_store_pd_a(&r.z.x, _mm_add_pd(_mm_load_pd_a(&m.z.x), ms));
	_mm_store_pd_a(&r.w.x, _mm_add_pd(_mm_load_pd_a(&m.w.x), ms));
	_mm_store_pd_a(&r.x.z, _mm_add_pd(_mm_load_pd_a(&m.x.z), ms));
	_mm_store_pd_a(&r.y.z, _mm_add_pd(_mm_load_pd_a(&m.y.z), ms));
	_mm_store_pd_a(&r.z.z, _mm_add_pd(_mm_load_pd_a(&m.z.z), ms));
	_mm_store_pd_a(&r.w.z, _mm_add_pd(_mm_load_pd_a(&m.w.z), ms));
	return r;
}

mathsimd_func(mat4_t,mat4_add_transpose_sse2)(mat4_t l, mat4_t r)
{
	__m128d mx_xy;
	__m128d my_xy;
	__m128d mz_xy;
	__m128d mw_xy;

	__m128d mx_zw;
	__m128d my_zw;
	__m128d mz_zw;
	__m128d mw_zw;

	__m128d t;
	mat4_t o;
	
	mx_xy = _mm_load_pd_a(&r.x.x);
	my_xy = _mm_load_pd_a(&r.y.x);
	mz_xy = _mm_load_pd_a(&r.z.x);
	mw_xy = _mm_load_pd_a(&r.w.x);

	t = _mm_unpacklo_pd(mx_xy, my_xy);
	_mm_store_pd_a(&o.x.x, _mm_add_pd(_mm_load_pd_a(&l.x.x), t));
	
	t = _mm_unpacklo_pd(mz_xy, mw_xy);
	_mm_store_pd_a(&o.x.z, _mm_add_pd(_mm_load_pd_a(&l.x.z), t));
	
	t = _mm_unpackhi_pd(mx_xy, my_xy);
	_mm_store_pd_a(&o.y.x, _mm_add_pd(_mm_load_pd_a(&l.y.x), t));
	
	t = _mm_unpackhi_pd(mz_xy, mw_xy);
	_mm_store_pd_a(&o.y.z, _mm_add_pd(_mm_load_pd_a(&l.y.z), t));
	
	mx_zw = _mm_load_pd_a(&r.x.z);
	my_zw = _mm_load_pd_a(&r.y.z);
	mz_zw = _mm_load_pd_a(&r.z.z);
	mw_zw = _mm_load_pd_a(&r.w.z);

	t = _mm_unpacklo_pd(mx_zw, my_zw);
	_mm_store_pd_a(&o.z.x, _mm_add_pd(_mm_load_pd_a(&l.z.x), t));
	
	t = _mm_unpacklo_pd(mz_zw, mw_zw);
	_mm_store_pd_a(&o.z.z, _mm_add_pd(_mm_load_pd_a(&l.z.z), t));
	
	t = _mm_unpackhi_pd(mx_zw, my_zw);
	_mm_store_pd_a(&o.w.x, _mm_add_pd(_mm_load_pd_a(&l.w.x), t));
	
	t = _mm_unpackhi_pd(mz_zw, mw_zw);
	_mm_store_pd_a(&o.w.z, _mm_add_pd(_mm_load_pd_a(&l.w.z), t));

	return o;
}

mathsimd_func(mat4_t,mat4_sub_sse2)(mat4_t l, mat4_t r)
{
	mat4_t o;
	_mm_store_pd_a(&o.x.x, _mm_sub_pd(_mm_load_pd_a(&l.x.x), _mm_load_pd_a(&r.x.x)));
	_mm_store_pd_a(&o.y.x, _mm_sub_pd(_mm_load_pd_a(&l.y.x), _mm_load_pd_a(&r.y.x)));
	_mm_store_pd_a(&o.z.x, _mm_sub_pd(_mm_load_pd_a(&l.z.x), _mm_load_pd_a(&r.z.x)));
	_mm_store_pd_a(&o.w.x, _mm_sub_pd(_mm_load_pd_a(&l.w.x), _mm_load_pd_a(&r.w.x)));
	_mm_store_pd_a(&o.x.z, _mm_sub_pd(_mm_load_pd_a(&l.x.z), _mm_load_pd_a(&r.x.z)));
	_mm_store_pd_a(&o.y.z, _mm_sub_pd(_mm_load_pd_a(&l.y.z), _mm_load_pd_a(&r.y.z)));
	_mm_store_pd_a(&o.z.z, _mm_sub_pd(_mm_load_pd_a(&l.z.z), _mm_load_pd_a(&r.z.z)));
	_mm_store_pd_a(&o.w.z, _mm_sub_pd(_mm_load_pd_a(&l.w.z), _mm_load_pd_a(&r.w.z)));
	return o;
}

mathsimd_func(mat4_t,mat4_sub_s_sse2)(mat4_t m, real_t s)
{
	__m128d ms = _mm_load1_pd(&s);
	mat4_t r;
	_mm_store_pd_a(&r.x.x, _mm_sub_pd(_mm_load_pd_a(&m.x.x), ms));
	_mm_store_pd_a(&r.y.x, _mm_sub_pd(_mm_load_pd_a(&m.y.x), ms));
	_mm_store_pd_a(&r.z.x, _mm_sub_pd(_mm_load_pd_a(&m.z.x), ms));
	_mm_store_pd_a(&r.w.x, _mm_sub_pd(_mm_load_pd_a(&m.w.x), ms));
	_mm_store_pd_a(&r.x.z, _mm_sub_pd(_mm_load_pd_a(&m.x.z), ms));
	_mm_store_pd_a(&r.y.z, _mm_sub_pd(_mm_load_pd_a(&m.y.z), ms));
	_mm_store_pd_a(&r.z.z, _mm_sub_pd(_mm_load_pd_a(&m.z.z), ms));
	_mm_store_pd_a(&r.w.z, _mm_sub_pd(_mm_load_pd_a(&m.w.z), ms));
	return r;
}

mathsimd_func(mat4_t,mat4_sub_transpose_sse2)(mat4_t l, mat4_t r)
{
	__m128d mx_xy;
	__m128d my_xy;
	__m128d mz_xy;
	__m128d mw_xy;

	__m128d mx_zw;
	__m128d my_zw;
	__m128d mz_zw;
	__m128d mw_zw;

	__m128d t;
	mat4_t o;
	
	mx_xy = _mm_load_pd_a(&r.x.x);
	my_xy = _mm_load_pd_a(&r.y.x);
	mz_xy = _mm_load_pd_a(&r.z.x);
	mw_xy = _mm_load_pd_a(&r.w.x);

	t = _mm_unpacklo_pd(mx_xy, my_xy);
	_mm_store_pd_a(&o.x.x, _mm_sub_pd(_mm_load_pd_a(&l.x.x), t));
	
	t = _mm_unpacklo_pd(mz_xy, mw_xy);
	_mm_store_pd_a(&o.x.z, _mm_sub_pd(_mm_load_pd_a(&l.x.z), t));
	
	t = _mm_unpackhi_pd(mx_xy, my_xy);
	_mm_store_pd_a(&o.y.x, _mm_sub_pd(_mm_load_pd_a(&l.y.x), t));
	
	t = _mm_unpackhi_pd(mz_xy, mw_xy);
	_mm_store_pd_a(&o.y.z, _mm_sub_pd(_mm_load_pd_a(&l.y.z), t));
	
	mx_zw = _mm_load_pd_a(&r.x.z);
	my_zw = _mm_load_pd_a(&r.y.z);
	mz_zw = _mm_load_pd_a(&r.z.z);
	mw_zw = _mm_load_pd_a(&r.w.z);

	t = _mm_unpacklo_pd(mx_zw, my_zw);
	_mm_store_pd_a(&o.z.x, _mm_sub_pd(_mm_load_pd_a(&l.z.x), t));
	
	t = _mm_unpacklo_pd(mz_zw, mw_zw);
	_mm_store_pd_a(&o.z.z, _mm_sub_pd(_mm_load_pd_a(&l.z.z), t));
	
	t = _mm_unpackhi_pd(mx_zw, my_zw);
	_mm_store_pd_a(&o.w.x, _mm_sub_pd(_mm_load_pd_a(&l.w.x), t));
	
	t = _mm_unpackhi_pd(mz_zw, mw_zw);
	_mm_store_pd_a(&o.w.z, _mm_sub_pd(_mm_load_pd_a(&l.w.z), t));
	return o;
}

mathsimd_func(mat4_t,mat4_mul_sse2)(mat4_t l, mat4_t r)
{
	__m128d mrx_xy = _mm_load_pd_a(&r.x.x);
	__m128d mry_xy = _mm_load_pd_a(&r.y.x);
	__m128d mrz_xy = _mm_load_pd_a(&r.z.x);
	__m128d mrw_xy = _mm_load_pd_a(&r.w.x);
	__m128d mrx_zw = _mm_load_pd_a(&r.x.x);
	__m128d mry_zw = _mm_load_pd_a(&r.y.x);
	__m128d mrz_zw = _mm_load_pd_a(&r.z.x);
	__m128d mrw_zw = _mm_load_pd_a(&r.w.x);
	__m128d mvxy;
	__m128d mvzw;
	
	mat4_t o;
	mvxy = _mm_load_pd_a(&l.x.x);
	mvzw = _mm_load_pd_a(&l.x.z);
	_mm_store_pd_a(&o.x.x, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), mrx_xy),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), mry_xy)),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), mrz_xy),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), mrw_xy))));
	_mm_store_pd_a(&o.x.z, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), mrx_zw),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), mry_zw)),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), mrz_zw),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), mrw_zw))));
	
	mvxy = _mm_load_pd_a(&l.y.x);
	mvzw = _mm_load_pd_a(&l.y.z);
	_mm_store_pd_a(&o.y.x, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), mrx_xy),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), mry_xy)),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), mrz_xy),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), mrw_xy))));
	_mm_store_pd_a(&o.y.z, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), mrx_zw),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), mry_zw)),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), mrz_zw),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), mrw_zw))));
	
	mvxy = _mm_load_pd_a(&l.z.x);
	mvzw = _mm_load_pd_a(&l.z.z);
	_mm_store_pd_a(&o.z.x, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), mrx_xy),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), mry_xy)),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), mrz_xy),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), mrw_xy))));
	_mm_store_pd_a(&o.z.z, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), mrx_zw),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), mry_zw)),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), mrz_zw),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), mrw_zw))));
	
	mvxy = _mm_load_pd_a(&l.w.x);
	mvzw = _mm_load_pd_a(&l.w.z);
	_mm_store_pd_a(&o.w.x, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), mrx_xy),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), mry_xy)),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), mrz_xy),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), mrw_xy))));
	_mm_store_pd_a(&o.w.z, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), mrx_zw),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), mry_zw)),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), mrz_zw),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), mrw_zw))));

	return o;
}

mathsimd_func(mat4_t,mat4_mul_transpose_sse2)(mat4_t l, mat4_t r)
{
	__m128d mrx_xy = _mm_load_pd_a(&r.x.x);
	__m128d mry_xy = _mm_load_pd_a(&r.y.x);
	__m128d mrz_xy = _mm_load_pd_a(&r.z.x);
	__m128d mrw_xy = _mm_load_pd_a(&r.w.x);
	__m128d mrx_zw = _mm_load_pd_a(&r.x.x);
	__m128d mry_zw = _mm_load_pd_a(&r.y.x);
	__m128d mrz_zw = _mm_load_pd_a(&r.z.x);
	__m128d mrw_zw = _mm_load_pd_a(&r.w.x);
	__m128d mvxy;
	__m128d mvzw;
	
	mat4_t o;
	mvxy = _mm_load_pd_a(&l.x.x);
	mvzw = _mm_load_pd_a(&l.x.z);
	_mm_store_pd_a(&o.x.x, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mrx_xy, mry_xy)),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mrx_xy, mry_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mrx_zw, mry_zw)),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mrx_zw, mry_zw)))));
	_mm_store_pd_a(&o.x.z, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mrz_xy, mrw_xy)),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mrz_xy, mrw_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mrz_zw, mrw_zw)),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mrz_zw, mrw_zw)))));
	
	mvxy = _mm_load_pd_a(&l.y.x);
	mvzw = _mm_load_pd_a(&l.y.z);
	_mm_store_pd_a(&o.y.x, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mrx_xy, mry_xy)),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mrx_xy, mry_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mrx_zw, mry_zw)),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mrx_zw, mry_zw)))));
	_mm_store_pd_a(&o.y.z, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mrz_xy, mrw_xy)),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mrz_xy, mrw_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mrz_zw, mrw_zw)),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mrz_zw, mrw_zw)))));
	
	mvxy = _mm_load_pd_a(&l.z.x);
	mvzw = _mm_load_pd_a(&l.z.z);
	_mm_store_pd_a(&o.z.x, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mrx_xy, mry_xy)),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mrx_xy, mry_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mrx_zw, mry_zw)),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mrx_zw, mry_zw)))));
	_mm_store_pd_a(&o.z.z, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mrz_xy, mrw_xy)),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mrz_xy, mrw_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mrz_zw, mrw_zw)),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mrz_zw, mrw_zw)))));
	
	mvxy = _mm_load_pd_a(&l.w.x);
	mvzw = _mm_load_pd_a(&l.w.z);
	_mm_store_pd_a(&o.w.x, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mrx_xy, mry_xy)),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mrx_xy, mry_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mrx_zw, mry_zw)),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mrx_zw, mry_zw)))));
	_mm_store_pd_a(&o.w.z, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mrz_xy, mrw_xy)),_mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mrz_xy, mrw_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mrz_zw, mrw_zw)),_mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mrz_zw, mrw_zw)))));

	return o;
}

#if MATHUTIL_DETECT_SIMD
int mathutil_sse2_implements()
{
	if(!CPUID_SSE2()) return 0;

	vec4_abs = vec4_abs_sse2;
	vec4_sgn = vec4_sgn_sse2;
	vec4_invert = vec4_invert_sse2;
	vec4_length = vec4_length_sse2;
	vec4_scale = vec4_scale_sse2;
	vec4_dot = vec4_dot_sse2;
	vec4_add = vec4_add_sse2;
	vec4_sub = vec4_sub_sse2;
	vec4_mul = vec4_mul_sse2;
	vec4_div = vec4_div_sse2;
	vec4_mul_mat4_transpose = vec4_mul_mat4_transpose_sse2;
	vec4_normalize = vec4_normalize_sse2;
	
	vec4_mul_mat4 = vec4_mul_mat4_sse2;

	vec4_min = vec4_min_sse2;
	vec4_max = vec4_max_sse2;
	vec4_lerp = vec4_lerp_sse2;
	vec4_slerp = vec4_slerp_sse2;

	quat_mul = quat_mul_sse2;
	quat_add_vec = quat_add_vec_sse2;

	mat4_transpose = mat4_transpose_sse2;
	mat4_add = mat4_add_sse2;
	mat4_add_s = mat4_add_s_sse2;
	mat4_add_transpose = mat4_add_transpose_sse2;
	mat4_sub = mat4_sub_sse2;
	mat4_sub_s = mat4_sub_s_sse2;
	mat4_sub_transpose = mat4_sub_transpose_sse2;
	mat4_mul = mat4_mul_sse2;
	mat4_mul_transpose = mat4_mul_transpose_sse2;
	return 1;
}
#endif // MATHUTIL_DETECT_SIMD

#endif // !MATHUTIL_USE_DOUBLE