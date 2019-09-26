#include"mathutil_sse.h"
#include"mathutil_ref.h"
#include"cpudetect.h"
#include<xmmintrin.h>
#include<immintrin.h>

#if defined(_MSC_VER)
#  if !__INTELLISENSE__
#    define _compiler_barrier _ReadWriteBarrier()
#  else
#    define _compiler_barrier
#  endif
#elif defined(__GNUC__) || defined(__clang__)
#  define _compiler_barrier asm volatile("" ::: "memory")
#else
#  define _compiler_barrier (void)1
#endif

#ifndef _MSC_VER
#  define _mm_castps_si128(a) (__m128i)(a)
#  define _mm_castsi128_ps(a) (__m128)(a)
#endif

#if MATHUTIL_VAR_NOT_ALIGNED && !MATHUTIL_VAR_ASSUME_ALIGNED
#define _mm_load_ps_a _mm_loadu_ps
#define _mm_load_pd_a _mm_loadu_pd
#define _mm_store_ps_a _mm_storeu_ps
#define _mm_store_pd_a _mm_storeu_pd
#else // !MATHUTIL_VAR_NOT_ALIGNED || MATHUTIL_VAR_ASSUME_ALIGNED
#define _mm_load_ps_a _mm_load_ps
#define _mm_load_pd_a _mm_load_pd
#define _mm_store_ps_a _mm_store_ps
#define _mm_store_pd_a _mm_store_pd
#endif

#define mathsimd_func(r,n) r n ## _sse

static __m128 vec4_load(vec4_t v)
{
#if VEC4_WITH_M128_XYZW
	return v.m_xyzw;
#else
	return _mm_load_ps_a(&v.x);
#endif
}

static void vec4_set_result(vec4_p v, __m128 m)
{
#if VEC4_WITH_M128_XYZW
	v->m_xyzw = m;
#else
	_mm_store_ps_a(&v->x, m);
#endif
}

#if !MATHUTIL_USE_DOUBLE

mathsimd_func(real_t, r_rnd)(uint32_t *p_seed)
{
	return r_rnd_ref(p_seed);
}

mathsimd_func(real_t,r_sin)(real_t x)
{
	__m128 mfac3_5_7_9 = _mm_set_ps(
		-6.0f,
		 120.0f,
		-5040.0f,
		 362880.0f);
	__m128 mfac11_13_15_17 = _mm_set_ps(
		-39916800.0f,
		 6227020800.0f,
		-1307674368000.0f,
		 355687428096000.0f);
	__m128 m1111 = _mm_set1_ps(1);
	__m128 m111x;
	__m128 m11xx;
	__m128 m1xxx;
	__m128 mxxxx;
	__m128 mxp2_4_6_8;
	__m128 mxp3_5_7_9;
	__m128 mxp11_13_15_17;
	__m128 mxp8_8_8_8;
	__m128 m1234;
	__m128 mr;
	real_t r;

	x = (r_floor_sse(x / (r_pi * 2.0f) + 0.5f) - 0.5f) * (r_pi * 2) - x;
	if(x > r_pi * 0.5f) x = x >= 0 ? r_pi - x : -r_pi - x;

	mxxxx = _mm_load1_ps(&x);
	m111x = _mm_move_ss(m1111,mxxxx);
	m11xx = _mm_movelh_ps(mxxxx, m1111);
	m1xxx = _mm_movelh_ps(mxxxx, m111x);
	
	mxp2_4_6_8 = _mm_mul_ps(mxxxx, m111x);
	mxp2_4_6_8 = _mm_mul_ps(mxp2_4_6_8, m11xx);
	mxp2_4_6_8 = _mm_mul_ps(mxp2_4_6_8, m1xxx);
	mxp2_4_6_8 = _mm_mul_ps(mxp2_4_6_8, mxp2_4_6_8);

	mxp8_8_8_8 = _mm_shuffle_ps(mxp2_4_6_8, mxp2_4_6_8, _MM_SHUFFLE(0,0,0,0));

	mxp3_5_7_9 = _mm_mul_ps(mxp2_4_6_8, mxxxx);
	mxp11_13_15_17 = _mm_mul_ps(mxp3_5_7_9, mxp8_8_8_8);

	m1234 = _mm_add_ps(_mm_div_ps(mxp3_5_7_9, mfac3_5_7_9),
		_mm_div_ps(mxp11_13_15_17, mfac11_13_15_17));

	mr = _mm_add_ps(m1234, _mm_shuffle_ps(m1234, m1234, _MM_SHUFFLE(2,3,0,1)));
	mr = _mm_add_ss(mxxxx, _mm_add_ss(mr, _mm_shuffle_ps(mr, mr, _MM_SHUFFLE(1,0,3,2))));

	_mm_store_ss(&r, mr);
	_mm_sfence();
	_compiler_barrier;
	return r;
}

mathsimd_func(real_t,r_cos)(real_t x)
{
	__m128 mfac2_4_6_8 = _mm_set_ps(
		-2.0f,
		 24.0f,
		-720.0f,
		 40320.0f);
	__m128 mfac10_12_14_16 = _mm_set_ps(
		-3628800.0f,
		 479001600.0f,
		-87178291200.0f,
		 20922789888000.0f);
	__m128 m1111 = _mm_set1_ps(1);
	__m128 m111x;
	__m128 m11xx;
	__m128 m1xxx;
	__m128 mxxxx;
	__m128 mxp2_4_6_8;
	__m128 mxp10_12_14_16;
	__m128 mxp8_8_8_8;
	__m128 m1234;
	__m128 mr;
	real_t r;

	x = x - (float)floor(x / (r_pi * 2)) * r_pi * 2;
	if(x > r_pi) x = r_pi * 2 - x;
	
	mxxxx = _mm_load1_ps(&x);
	m111x = _mm_move_ss(m1111,mxxxx);
	m11xx = _mm_movelh_ps(mxxxx, m1111);
	m1xxx = _mm_movelh_ps(mxxxx, m111x);
	
	mxp2_4_6_8 = _mm_mul_ps(mxxxx, m111x);
	mxp2_4_6_8 = _mm_mul_ps(mxp2_4_6_8, m11xx);
	mxp2_4_6_8 = _mm_mul_ps(mxp2_4_6_8, m1xxx);
	mxp2_4_6_8 = _mm_mul_ps(mxp2_4_6_8, mxp2_4_6_8);
	
	mxp8_8_8_8 = _mm_shuffle_ps(mxp2_4_6_8, mxp2_4_6_8, _MM_SHUFFLE(0,0,0,0));

	mxp10_12_14_16 = _mm_mul_ps(mxp2_4_6_8, mxp8_8_8_8);

	m1234 = _mm_add_ps(_mm_div_ps(mxp2_4_6_8, mfac2_4_6_8),
		_mm_div_ps(mxp10_12_14_16, mfac10_12_14_16));

	mr = _mm_add_ps(m1234, _mm_shuffle_ps(m1234, m1234, _MM_SHUFFLE(2,3,0,1)));
	mr = _mm_add_ss(m1111, _mm_add_ss(mr, _mm_shuffle_ps(mr, mr, _MM_SHUFFLE(1,0,3,2))));

	_mm_store_ss(&r, mr);
	_mm_sfence();
	_compiler_barrier;
	return r;
}

mathsimd_func(real_t, r_tan)(real_t x)
{
	return r_tan_ref(x);
}

mathsimd_func(real_t, r_abs)(real_t x)
{
	return r_abs_ref(x);
}

mathsimd_func(real_t, r_sgn)(real_t x)
{
	return r_sgn_ref(x);
}

mathsimd_func(real_t, r_sqr)(real_t x)
{
	return r_sqr_ref(x);
}

mathsimd_func(real_t, r_floor)(real_t x)
{
	return r_floor_ref(x);
}

mathsimd_func(real_t, r_ceil)(real_t x)
{
	return r_ceil_ref(x);
}

mathsimd_func(real_t, r_atan)(real_t x)
{
	return r_atan_ref(x);
}

mathsimd_func(real_t, r_exp)(real_t x)
{
	return r_exp_ref(x);
}

mathsimd_func(real_t,r_log)(real_t x)
{
	__m128 mpowpart;
	__m128 mdivpart = _mm_set_ps(-2, 3,-4, 5);
	__m128 m1111 = _mm_set1_ps(1);
	__m128 mnpnp = _mm_set_ps(-1, 1,-1, 1);
	__m128 m111x;
	__m128 m11xx;
	__m128 m1xxx;
	__m128 mxxxx;
	__m128 mx4x4x4x4;
	__m128 mr = _mm_set_ss(-1);
	__m128 ma;
	__m128 me = _mm_set1_ps(r_epsilon);
	real_t r;
	const unsigned repeat_min = 2;
	const unsigned repeat_max = 8;
	unsigned i;
	
	mxxxx = _mm_load1_ps(&x);
	m111x = _mm_move_ss(m1111,mxxxx);
	m11xx = _mm_movelh_ps(mxxxx, m1111);
	m1xxx = _mm_movelh_ps(mxxxx, m111x);
	mr = _mm_add_ss(mr, mxxxx); // x+1, 0, 0, 0

	mpowpart = _mm_mul_ps(m111x, _mm_mul_ps(m11xx, _mm_mul_ps(m1xxx, mxxxx)));
	mx4x4x4x4 = _mm_shuffle_ps(mpowpart, mpowpart, _MM_SHUFFLE(0,0,0,0));
	mpowpart = _mm_mul_ps(mpowpart, mxxxx); // x^2, x^3, x^4, x^5
	ma = _mm_div_ps(mpowpart, mdivpart);
	mr = _mm_add_ps(mr, ma);

	for(i = 0; i < repeat_min; i++)
	{
		mpowpart = _mm_mul_ps(mpowpart, mx4x4x4x4);
		mdivpart = _mm_add_ps(mdivpart, mnpnp);
		ma = _mm_div_ps(mpowpart, mdivpart);
		mr = _mm_add_ps(mr, ma);
	}

	for(i = 0; i < repeat_max - repeat_min; i++)
	{
		mpowpart = _mm_mul_ps(mpowpart, mx4x4x4x4);
		mdivpart = _mm_add_ps(mdivpart, mnpnp);
		ma = _mm_div_ps(mpowpart, mdivpart);
		mr = _mm_add_ps(mr, ma);

		if(_mm_ucomilt_ss(ma, me)) break;
	}

	mr = _mm_add_ps(mr, _mm_shuffle_ps(mr, mr, _MM_SHUFFLE(2,3,0,1)));
	mr = _mm_add_ss(mr, _mm_shuffle_ps(mr, mr, _MM_SHUFFLE(1,0,3,2)));
	
	_mm_store_ss(&r, mr);
	_mm_sfence();
	_compiler_barrier;
	return r;
}

mathsimd_func(real_t,r_pow)(real_t x, real_t y)
{
	return r_exp_sse(y * r_log_sse(x));
}

mathsimd_func(real_t, r_mod)(real_t x, real_t y)
{
	return r_mod_ref(x, y);
}

mathsimd_func(real_t, r_max)(real_t x, real_t y)
{
	return r_max_ref(x, y);
}

mathsimd_func(real_t, r_min)(real_t x, real_t y)
{
	return r_min_ref(x, y);
}

mathsimd_func(real_t, r_atan2)(real_t y, real_t x)
{
	return r_atan2_ref(y, x);
}

mathsimd_func(real_t, r_clamp)(real_t n, real_t min_, real_t max_)
{
	return r_clamp_ref(n, min_, max_);
}

mathsimd_func(real_t, r_lerp)(real_t a, real_t b, real_t s)
{
	return r_lerp_ref(a, b, s);
}

mathsimd_func(real_t, r_hermite)(real_t s)
{
	return r_hermite_ref(s);
}

mathsimd_func(real_t, r_slerp)(real_t a, real_t b, real_t s)
{
	return r_slerp_ref(a, b, s);
}

mathsimd_func(vec4_t, vec4)(real_t x, real_t y, real_t z, real_t w)
{
	vec4_t v;
	v.x = x;
	v.y = y;
	v.z = z;
	v.w = w;
#if VEC4_WITH_M128_XYZW
	v.m_xyzw = _mm_set_ps(x, y, z, w);
#endif
	return v;
}

mathsimd_func(vec4_t, vec4_flushcomp)(vec4_t v)
{
#if VEC4_WITH_M128_XYZW
	_mm_store_ps_a(&v.x, v.m_xyzw);
#endif
	return v;
}

static union _abs_mask_union
{
	uint32_t i;
	float f;
}_abs_mask = { 0x7fffffff };

mathsimd_func(vec4_t,vec4_abs)(vec4_t v)
{
	vec4_t r;
	__m128 _abs_mask_reg = _mm_set1_ps(_abs_mask.f);
	vec4_set_result(&r, _mm_and_ps(vec4_load(v), _abs_mask_reg));
	return r;
}

mathsimd_func(vec4_t,vec4_sgn)(vec4_t v)
{
	vec4_t r;
	__m128 mz = _mm_setzero_ps();
	__m128 mv = vec4_load(v);
    __m128 mp = _mm_and_ps(_mm_cmpgt_ps(mv, mz), _mm_set1_ps(1.0f));
    __m128 mn = _mm_and_ps(_mm_cmplt_ps(mv, mz), _mm_set1_ps(-1.0f));
	vec4_set_result(&r, _mm_or_ps(mp, mn));
	return r;
}

mathsimd_func(vec4_t,vec4_invert)(vec4_t v)
{
	return vec4_scale_sse(v, -1);
}

mathsimd_func(real_t,vec4_length)(vec4_t v)
{
	__m128 mv = vec4_load(v);
	__m128 m = _mm_mul_ps(mv, mv);
	__m128 t = _mm_add_ps(m, _mm_shuffle_ps(m, m, _MM_SHUFFLE(2, 3, 0, 1)));
	__m128 length = _mm_sqrt_ss(_mm_add_ps(t, _mm_shuffle_ps(t, t, _MM_SHUFFLE(1, 0, 3, 2))));
	real_t r;
	_mm_store_ss(&r, length);
	_mm_sfence();
	_compiler_barrier;
	return r;
}

mathsimd_func(vec4_t,vec4_normalize)(vec4_t v)
{
	vec4_t r;
	__m128 mv = vec4_load(v);
	__m128 m = _mm_mul_ps(mv, mv);
	__m128 t = _mm_add_ps(m, _mm_shuffle_ps(m, m, _MM_SHUFFLE(2, 3, 0, 1)));
	vec4_set_result(&r, _mm_mul_ps(mv, _mm_rsqrt_ps(_mm_add_ps(t, _mm_shuffle_ps(t, t, _MM_SHUFFLE(1, 0, 3, 2))))));
	return r;
}

mathsimd_func(vec4_t,vec4_scale)(vec4_t v, real_t s)
{
	vec4_t r;
	vec4_set_result(&r, _mm_mul_ps(vec4_load(v), _mm_load1_ps(&s)));
	return r;
}

mathsimd_func(vec4_t, vec4_pow) (vec4_t v, real_t n)
{
	vec4_t r;
	vec4_flushcomp(v);
	r.x = r_pow(v.x, n);
	r.y = r_pow(v.y, n);
	r.z = r_pow(v.z, n);
	r.w = r_pow(v.w, n);
#if VEC4_WITH_M128_XYZW
	r.m_xyzw = _mm_load_ps_a(&r.x);
#endif
	return r;
}

mathsimd_func(real_t,vec4_dot)(vec4_t v1, vec4_t v2)
{
	/*
	__m128 m = _mm_mul_ps(vec4_load(v1), vec4_load(v2));
	__m128 t = _mm_add_ps(_mm_unpacklo_ps(m, m), _mm_unpackhi_ps(m, m));
	real_t r;
	_mm_store_ss(&r, _mm_add_ss(t, _mm_movehl_ps(t, t)));
	_mm_sfence();
	_compiler_barrier;
	return r;
	*/

	__m128 m = _mm_mul_ps(vec4_load(v1), vec4_load(v2));
	__m128 t = _mm_add_ps(m, _mm_shuffle_ps(m, m, _MM_SHUFFLE(2, 3, 0, 1)));
	real_t r;
	_mm_store_ss(&r, _mm_add_ss(t, _mm_shuffle_ps(t, t, _MM_SHUFFLE(1, 0, 3, 2))));
	_mm_sfence();
	_compiler_barrier;
	return r;
}

mathsimd_func(vec4_t,vec4_add)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	vec4_set_result(&r, _mm_add_ps(vec4_load(v1), vec4_load(v2)));
	return r;
}

mathsimd_func(vec4_t,vec4_sub)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	vec4_set_result(&r, _mm_sub_ps(vec4_load(v1), vec4_load(v2)));
	return r;
}

mathsimd_func(vec4_t,vec4_mul)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	vec4_set_result(&r, _mm_mul_ps(vec4_load(v1), vec4_load(v2)));
	return r;
}

mathsimd_func(vec4_t,vec4_div)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	vec4_set_result(&r, _mm_div_ps(vec4_load(v1), vec4_load(v2)));
	return r;
}

mathsimd_func(vec4_t,vec4_min)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	vec4_set_result(&r, _mm_min_ps(vec4_load(v1), vec4_load(v2)));
	return r;
}

mathsimd_func(vec4_t,vec4_max)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	vec4_set_result(&r, _mm_max_ps(vec4_load(v1), vec4_load(v2)));
	return r;
}

mathsimd_func(vec4_t,vec4_clamp)(vec4_t v, real_t min_, real_t max_)
{
	vec4_t r;
	__m128 ml, mh;
	
	ml = _mm_load1_ps(&min_);
	mh = _mm_load1_ps(&max_);
	
	vec4_set_result(&r, _mm_min_ps(_mm_max_ps(vec4_load(v), ml), mh));
	return r;
}

mathsimd_func(vec4_t,vec4_cross3)(vec4_t v1, vec4_t v2)
{
	__m128 mv1 = vec4_load(v1);
	__m128 mv2 = vec4_load(v2);
	__m128 ml;
	__m128 mr;
	vec4_t r;
	
	ml = _mm_mul_ps(_mm_shuffle_ps(mv1,mv1,_MM_SHUFFLE(3,0,2,1)), _mm_shuffle_ps(mv2,mv2,_MM_SHUFFLE(3,1,0,2)));
	mr = _mm_mul_ps(_mm_shuffle_ps(mv1,mv1,_MM_SHUFFLE(3,1,0,2)), _mm_shuffle_ps(mv2,mv2,_MM_SHUFFLE(3,0,2,1)));

	vec4_set_result(&r, _mm_sub_ps(ml, mr));
	return r;
}

mathsimd_func(vec4_t,vec4_rot_quat)(vec4_t v, quat_t q)
{
	vec4_t r;
	__m128 mv = _mm_mul_ps(vec4_load(v), _mm_set_ps(0,1,1,1));
	__m128 mq = vec4_load(q);
	__m128 mqconj = _mm_mul_ps(mq, _mm_set_ps(1,-1,-1,-1));
	__m128 ms2 = _mm_set_ps(-1, 1,-1, 1);
	__m128 ms3 = _mm_set_ps(-1,-1, 1, 1);
	__m128 ms4 = _mm_set_ps(-1, 1, 1,-1);
	__m128 mr1;
	__m128 mr2;

	mr1 = _mm_mul_ps(_mm_shuffle_ps(mq, mq, _MM_SHUFFLE(3,3,3,3)), mv);
	mr1 = _mm_add_ps(mr1, _mm_mul_ps(ms2, _mm_mul_ps(_mm_shuffle_ps(mq, mq, _MM_SHUFFLE(0,0,0,0)), _mm_shuffle_ps(mv, mv, _MM_SHUFFLE(0,1,2,3)))));
	mr1 = _mm_add_ps(mr1, _mm_mul_ps(ms3, _mm_mul_ps(_mm_shuffle_ps(mq, mq, _MM_SHUFFLE(1,1,1,1)), _mm_shuffle_ps(mv, mv, _MM_SHUFFLE(1,0,3,2)))));
	mr1 = _mm_add_ps(mr1, _mm_mul_ps(ms4, _mm_mul_ps(_mm_shuffle_ps(mq, mq, _MM_SHUFFLE(2,2,2,2)), _mm_shuffle_ps(mv, mv, _MM_SHUFFLE(2,3,0,1)))));
	
	mr2 = _mm_mul_ps(_mm_shuffle_ps(mr1, mr1, _MM_SHUFFLE(3,3,3,3)), mqconj);
	mr2 = _mm_add_ps(mr2, _mm_mul_ps(ms2, _mm_mul_ps(_mm_shuffle_ps(mr1, mr1, _MM_SHUFFLE(0,0,0,0)), _mm_shuffle_ps(mqconj, mqconj, _MM_SHUFFLE(0,1,2,3)))));
	mr2 = _mm_add_ps(mr2, _mm_mul_ps(ms3, _mm_mul_ps(_mm_shuffle_ps(mr1, mr1, _MM_SHUFFLE(1,1,1,1)), _mm_shuffle_ps(mqconj, mqconj, _MM_SHUFFLE(1,0,3,2)))));
	mr2 = _mm_add_ps(mr2, _mm_mul_ps(ms4, _mm_mul_ps(_mm_shuffle_ps(mr1, mr1, _MM_SHUFFLE(2,2,2,2)), _mm_shuffle_ps(mqconj, mqconj, _MM_SHUFFLE(2,3,0,1)))));
	
	vec4_set_result(&r, mr2);
	return r;
}

mathsimd_func(vec4_t,vec4_mul_mat4)(vec4_t v, mat4_t m)
{
	vec4_t r;

	__m128 mv = vec4_load(v);
	__m128 mr;
	
	mr =                _mm_mul_ps(_mm_shuffle_ps(mv, mv, _MM_SHUFFLE(0,0,0,0)), vec4_load(m.x));
	mr = _mm_add_ps(mr, _mm_mul_ps(_mm_shuffle_ps(mv, mv, _MM_SHUFFLE(1,1,1,1)), vec4_load(m.y)));
	mr = _mm_add_ps(mr, _mm_mul_ps(_mm_shuffle_ps(mv, mv, _MM_SHUFFLE(2,2,2,2)), vec4_load(m.z)));
	mr = _mm_add_ps(mr, _mm_mul_ps(_mm_shuffle_ps(mv, mv, _MM_SHUFFLE(3,3,3,3)), vec4_load(m.w)));

	vec4_set_result(&r, mr);
	return r;
}

mathsimd_func(vec4_t,vec4_mul_mat4_transpose)(vec4_t v, mat4_t m)
{
	vec4_t r;

	__m128 mv = vec4_load(v);
	
	__m128 mx = _mm_mul_ps(mv, vec4_load(m.x));
	__m128 my = _mm_mul_ps(mv, vec4_load(m.y));
	__m128 mz = _mm_mul_ps(mv, vec4_load(m.z));
	__m128 mw = _mm_mul_ps(mv, vec4_load(m.w));
	
	_MM_TRANSPOSE4_PS(mx, my, mz, mw);
	
	vec4_set_result(&r, _mm_add_ps(_mm_add_ps(mx, my), _mm_add_ps(mz, mw)));
	return r;
}

mathsimd_func(vec4_t,vec4_lerp)(vec4_t v1, vec4_t v2, real_t s)
{
	__m128 ma = vec4_load(v1);
	__m128 mb = vec4_load(v2);

	vec4_t r;

	vec4_set_result(&r, _mm_add_ps(ma, _mm_mul_ps(_mm_sub_ps(mb, ma), _mm_load1_ps(&s))));
	return r;
}

mathsimd_func(vec4_t,vec4_slerp)(vec4_t v1, vec4_t v2, real_t s)
{
	__m128 ma = vec4_load(v1);
	__m128 mb = vec4_load(v2);
	__m128 ms = _mm_load1_ps(&s);
	vec4_t r;
	
	ms = _mm_mul_ps(_mm_mul_ps(ms, ms), _mm_sub_ps(_mm_set1_ps(3), _mm_mul_ps(_mm_set1_ps(2), ms)));
	ms = _mm_max_ps(_mm_setzero_ps(), _mm_min_ps(_mm_set1_ps(1), ms));
	vec4_set_result(&r, _mm_add_ps(ma, _mm_mul_ps(_mm_sub_ps(mb, ma), ms)));

	return r;
}
mathsimd_func(quat_t, quat)(real_t x, real_t y, real_t z, real_t w)
{
	quat_t q;
	q.x = x;
	q.y = y;
	q.z = z;
	q.w = w;
#if VEC4_WITH_M128_XYZW
	q.m_xyzw = _mm_set_ps(x, y, z, w);
#endif
	return q;
}

static __m128 quat_load(quat_t q)
{
#if VEC4_WITH_M128_XYZW
	return q.m_xyzw;
#else
	return _mm_load_ps_a(&q.x);
#endif
}

static void quat_set_result(quat_p q, __m128 m)
{
#if VEC4_WITH_M128_XYZW
	q->m_xyzw = m;
#else
	_mm_store_ps_a(&q->x, m);
#endif
}

mathsimd_func(quat_t, quat_flushcomp)(quat_t q)
{
#if VEC4_WITH_M128_XYZW
	_mm_store_ps_a(&q.x, q.m_xyzw);
#endif
	return q;
}

mathsimd_func(quat_t, quat_rot_axis)(vec4_t axis, real_t angle)
{
	real_t ha = angle / 2;
	real_t sin_ha = r_sin_sse(ha);
	real_t cos_ha = r_cos_sse(ha);
	__m128 mq = vec4_load(axis);
	quat_t q;

	mq = _mm_mul_ps(_mm_mul_ps(mq, _mm_set_ps(1, 1, 1, 0)), _mm_load1_ps(&sin_ha));
	mq = _mm_add_ps(_mm_mul_ps(_mm_set_ps(0, 0, 0, 1), _mm_load1_ps(&cos_ha)), mq);
	
	quat_set_result(&q, mq);
	return q;
}

mathsimd_func(quat_t,quat_mul)(quat_t q1, quat_t q2)
{
	__m128 mq1 = vec4_load(q1);
	__m128 mq2 = vec4_load(q2);
	__m128 mq;
	quat_t r;

	mq = _mm_mul_ps(_mm_shuffle_ps(mq1, mq1, _MM_SHUFFLE(3,3,3,3)), mq2);
	mq = _mm_add_ps(mq, _mm_mul_ps(_mm_set_ps(-1, 1,-1, 1), _mm_mul_ps(_mm_shuffle_ps(mq1, mq1, _MM_SHUFFLE(0,0,0,0)), _mm_shuffle_ps(mq2, mq2, _MM_SHUFFLE(0,1,2,3)))));
	mq = _mm_add_ps(mq, _mm_mul_ps(_mm_set_ps(-1,-1, 1, 1), _mm_mul_ps(_mm_shuffle_ps(mq1, mq1, _MM_SHUFFLE(1,1,1,1)), _mm_shuffle_ps(mq2, mq2, _MM_SHUFFLE(1,0,3,2)))));
	mq = _mm_add_ps(mq, _mm_mul_ps(_mm_set_ps(-1, 1, 1,-1), _mm_mul_ps(_mm_shuffle_ps(mq1, mq1, _MM_SHUFFLE(2,2,2,2)), _mm_shuffle_ps(mq2, mq2, _MM_SHUFFLE(2,3,0,1)))));

	vec4_set_result(&r, mq);
	return r;
}

mathsimd_func(quat_t,quat_add_vec)(quat_t q, vec4_t v, real_t s)
{
	__m128 mq1 = _mm_mul_ps(_mm_mul_ps(vec4_load(v), _mm_set_ps(0,1,1,1)), _mm_load1_ps(&s));
	__m128 mq2 = vec4_load(q);
	__m128 mq;
	quat_t r;

	// mq = _mm_mul_ps(_mm_shuffle_ps(mq1, mq1, _MM_SHUFFLE(3,3,3,3)), mq2);
	mq =                _mm_mul_ps(_mm_set_ps(-1, 1,-1, 1), _mm_mul_ps(_mm_shuffle_ps(mq1, mq1, _MM_SHUFFLE(0,0,0,0)), _mm_shuffle_ps(mq2, mq2, _MM_SHUFFLE(0,1,2,3))));
	mq = _mm_add_ps(mq, _mm_mul_ps(_mm_set_ps(-1,-1, 1, 1), _mm_mul_ps(_mm_shuffle_ps(mq1, mq1, _MM_SHUFFLE(1,1,1,1)), _mm_shuffle_ps(mq2, mq2, _MM_SHUFFLE(1,0,3,2)))));
	mq = _mm_add_ps(mq, _mm_mul_ps(_mm_set_ps(-1, 1, 1,-1), _mm_mul_ps(_mm_shuffle_ps(mq1, mq1, _MM_SHUFFLE(2,2,2,2)), _mm_shuffle_ps(mq2, mq2, _MM_SHUFFLE(2,3,0,1)))));

	vec4_set_result(&r, _mm_add_ps(_mm_mul_ps(mq, _mm_set1_ps(0.5)), mq2));
	return r;
}

mathsimd_func(mat4_t, mat4)(vec4_t x, vec4_t y, vec4_t z, vec4_t w)
{
	mat4_t r;
	r.x = x;
	r.y = y;
	r.z = z;
	r.w = w;
	return r;
}

mathsimd_func(mat4_t, mat4_flushcomp)(mat4_t m)
{
	mat4_t r;
	r.x = vec4_flushcomp(m.x);
	r.y = vec4_flushcomp(m.y);
	r.z = vec4_flushcomp(m.z);
	r.w = vec4_flushcomp(m.w);
	return r;
}

mathsimd_func(mat4_t, mat4_rot_x) (real_t angle)
{
	real_t sa = r_sin(angle);
	real_t ca = r_cos(angle);
	return mat4
	(
		vec4(1, 0, 0, 0),
		vec4(0, ca, sa, 0),
		vec4(0, -sa, ca, 0),
		vec4(0, 0, 0, 1)
	);
}

mathsimd_func(mat4_t, mat4_rot_y) (real_t angle)
{
	real_t sa = r_sin(angle);
	real_t ca = r_cos(angle);
	return mat4
	(
		vec4(ca, 0, -sa, 0),
		vec4(0, 1, 0, 0),
		vec4(sa, 0, ca, 0),
		vec4(0, 0, 0, 1)
	);
}

mathsimd_func(mat4_t, mat4_rot_z) (real_t angle)
{
	real_t sa = r_sin(angle);
	real_t ca = r_cos(angle);
	return mat4
	(
		vec4(ca, sa, 0, 0),
		vec4(-sa, ca, 0, 0),
		vec4(0, 0, 1, 0),
		vec4(0, 0, 0, 1)
	);
}

mathsimd_func(mat4_t, mat4_rot_axis)(vec4_t axis, real_t angle)
{
	// Cxx+c  Cxy+sz Cxz-sy 0
	// Cxy-sz Cyy+c  Cyz+sx 0
	// Cxz+sy Cyz-sx Czz+c  0
	// 0      0      0      1

	real_t sa = r_sin(angle);
	real_t ca = r_cos(angle);
	vec4_t v = vec4_normalize(axis);
	real_t C = 1 - ca;
	__m128 mC = _mm_load_ss(&C);
	__m128 ms = _mm_load_ss(&sa);
	__m128 mc = _mm_load_ss(&ca);
	__m128 mv = vec4_load(v);
	__m128 mxxx = _mm_shuffle_ps(mv, mv, _MM_SHUFFLE(3,0,0,0));
	__m128 myyy = _mm_shuffle_ps(mv, mv, _MM_SHUFFLE(3,1,1,1));
	__m128 mzzz = _mm_shuffle_ps(mv, mv, _MM_SHUFFLE(3,2,2,2));
	__m128 mCxyz, msxyzc, msxyzm;
	__m128 mrx, mry, mrz, mrw;
	mat4_t m;
	
	mC = _mm_shuffle_ps(mC, mC, _MM_SHUFFLE(3,0,0,0));
	mCxyz = _mm_mul_ps(mC, mv);
	msxyzc = _mm_mul_ps(_mm_shuffle_ps(ms, ms, _MM_SHUFFLE(3,0,0,0)), mv);
	msxyzc = _mm_add_ps(msxyzc, _mm_shuffle_ps(mc, mc, _MM_SHUFFLE(0,3,3,3)));
	msxyzm = _mm_mul_ps(msxyzc, _mm_set_ps(-1,-1,-1,1));
	mrx = _mm_add_ps(_mm_mul_ps(mCxyz, mxxx), _mm_shuffle_ps(msxyzc, msxyzm, _MM_SHUFFLE(3,1,2,3)));
	mry = _mm_add_ps(_mm_mul_ps(mCxyz, myyy), _mm_shuffle_ps(msxyzm, msxyzc, _MM_SHUFFLE(3,0,3,2)));
	mrz = _mm_add_ps(_mm_mul_ps(mCxyz, mzzz), _mm_mul_ps(_mm_shuffle_ps(msxyzc, msxyzc, _MM_SHUFFLE(3,3,0,1)), _mm_set_ps(1,-1,1,0)));
	mrw = _mm_set_ps(0, 0, 0, 1);

	vec4_set_result(&m.x, mrx);
	vec4_set_result(&m.y, mry);
	vec4_set_result(&m.z, mrz);
	vec4_set_result(&m.w, mrw);

	return m;
}

mathsimd_func(mat4_t, mat4_rot_euler)(real_t yaw, real_t pitch, real_t roll)
{
	real_t
		sr = r_sin(yaw),
		cr = r_cos(yaw),
		sp = r_sin(pitch),
		cp = r_cos(pitch),
		sy = r_sin(roll),
		cy = r_cos(roll);

	__m128 msrcrspcp = _mm_set_ps(sr, cr, sp, cp);
	__m128 msycysycy = _mm_set_ps(sy, cy, sy, cy);
	__m128 msrcp_srsp_crcp_crsp = _mm_mul_ps(
		_mm_shuffle_ps(msrcrspcp, msrcrspcp, _MM_SHUFFLE(1,1,0,0)),
		_mm_shuffle_ps(msrcrspcp, msrcrspcp, _MM_SHUFFLE(2,3,2,3)));

	//  crcy 0 -sycr 0
	__m128 alayerx = _mm_mul_ps(_mm_set_ps(1,0,-1,0), _mm_mul_ps(
		_mm_shuffle_ps(msrcrspcp, msrcrspcp, _MM_SHUFFLE(1,1,1,1)),
		_mm_shuffle_ps(msycysycy, msycysycy, _MM_SHUFFLE(0,1,0,1))));

	// -srcy 0  sysr 0
	__m128 alayery = _mm_mul_ps(_mm_set_ps(-1,1,0,0), _mm_mul_ps(
		_mm_shuffle_ps(msrcrspcp, msrcrspcp, _MM_SHUFFLE(0,0,0,0)),
		_mm_shuffle_ps(msycysycy, msycysycy, _MM_SHUFFLE(0,0,1,1))));

	//  sycp 0 cpcy 0
	__m128 alayerz = _mm_mul_ps(_mm_set_ps(1,0,1,0), _mm_mul_ps(
		_mm_shuffle_ps(msrcrspcp, msrcrspcp, _MM_SHUFFLE(2,3,2,3)),
		_mm_shuffle_ps(msycysycy, msycysycy, _MM_SHUFFLE(0,0,1,1))));

	__m128 sy_1_cy_0 = _mm_set_ps(sy, 1, cy, 0);
	__m128 z_msp_z_z = _mm_set_ps(0, -sp, 0, 0);

	__m128 resultx = _mm_add_ps(alayerx, _mm_mul_ps(_mm_shuffle_ps(msrcp_srsp_crcp_crsp, msrcp_srsp_crcp_crsp, _MM_SHUFFLE(3, 1, 0, 1)), sy_1_cy_0));
	__m128 resulty = _mm_add_ps(alayery, _mm_mul_ps(_mm_shuffle_ps(msrcp_srsp_crcp_crsp, msrcp_srsp_crcp_crsp, _MM_SHUFFLE(3, 3, 2, 3)), sy_1_cy_0));
	__m128 resultz = _mm_add_ps(alayerz, z_msp_z_z);
	__m128 resultw = _mm_set_ps(0, 0, 0, 1);

	mat4_t m;

	vec4_set_result(&m.x, resultx);
	vec4_set_result(&m.y, resulty);
	vec4_set_result(&m.z, resultz);
	vec4_set_result(&m.w, resultw);

	return m;
}

mathsimd_func(mat4_t,mat4_from_quat)(quat_t q)
{
	__m128 mq = vec4_load(q);
	__m128 mxxxx = _mm_shuffle_ps(mq, mq, _MM_SHUFFLE(0,0,0,0));
	__m128 myyyy = _mm_shuffle_ps(mq, mq, _MM_SHUFFLE(1,1,1,1));
	__m128 mzzzz = _mm_shuffle_ps(mq, mq, _MM_SHUFFLE(2,2,2,2));
	__m128 myxwz = _mm_shuffle_ps(mq, mq, _MM_SHUFFLE(2,3,0,1));
	__m128 mzwxy = _mm_shuffle_ps(mq, mq, _MM_SHUFFLE(1,0,3,2));
	__m128 mwzyx = _mm_shuffle_ps(mq, mq, _MM_SHUFFLE(0,1,2,3));
	__m128 mr1;
	__m128 mr2;
	mat4_t r;

	// 1 - 2yy - 2zz, 2xy + 2zw, 2xz - 2yw, 0
	// 2xy - 2zw, 1 - 2xx - 2zz, 2yz + 2xw, 0
	// 2xz + 2yw, 2yz - 2xw, 1 - 2xx - 2yy, 0
	// 0, 0, 0, 1

	// -yy  xy -yw  yyy * yxw
	// -zz  zw  xz  zzz * zwx

	//  xy -xx  xw  xxx * yxw
	// -zw -zz  yz  zzz * wzy
	
	//  xz -xw -xx  xxx * zwx
	//  yw  yz -yy  yyy * wzy

	mr1 = _mm_mul_ps(myyyy, myxwz);
	mr1 = _mm_mul_ps(mr1, _mm_set_ps(0,-2, 2,-2));
	mr2 = _mm_mul_ps(mzzzz, mzwxy);
	mr2 = _mm_mul_ps(mr2, _mm_set_ps(0, 2, 2,-2));
	vec4_set_result(&r.x, _mm_add_ps(_mm_set_ps(0,0,0,1), _mm_add_ps(mr1, mr2)));

	mr1 = _mm_mul_ps(mxxxx, myxwz);
	mr1 = _mm_mul_ps(mr1, _mm_set_ps(0, 2,-2, 2));
	mr2 = _mm_mul_ps(mzzzz, mwzyx);
	mr2 = _mm_mul_ps(mr2, _mm_set_ps(0, 2,-2,-2));
	vec4_set_result(&r.y, _mm_add_ps(_mm_set_ps(0,0,1,0), _mm_add_ps(mr1, mr2)));

	mr1 = _mm_mul_ps(mxxxx, mzwxy);
	mr1 = _mm_mul_ps(mr1, _mm_set_ps(0,-2,-2, 2));
	mr2 = _mm_mul_ps(myyyy, mwzyx);
	mr2 = _mm_mul_ps(mr2, _mm_set_ps(0,-2, 2, 2));
	vec4_set_result(&r.z, _mm_add_ps(_mm_set_ps(0,1,0,0), _mm_add_ps(mr1, mr2)));

	vec4_set_result(&r.w, _mm_set_ps(1,0,0,0));
	return r;
}

mathsimd_func(mat4_t,mat4_from_quat_transpose)(quat_t q)
{
	__m128 mq = vec4_load(q);
	__m128 mxxxx = _mm_shuffle_ps(mq, mq, _MM_SHUFFLE(0,0,0,0));
	__m128 myyyy = _mm_shuffle_ps(mq, mq, _MM_SHUFFLE(1,1,1,1));
	__m128 mzzzz = _mm_shuffle_ps(mq, mq, _MM_SHUFFLE(2,2,2,2));
	__m128 myxwz = _mm_shuffle_ps(mq, mq, _MM_SHUFFLE(2,3,0,1));
	__m128 mzwxy = _mm_shuffle_ps(mq, mq, _MM_SHUFFLE(1,0,3,2));
	__m128 mwzyx = _mm_shuffle_ps(mq, mq, _MM_SHUFFLE(0,1,2,3));
	__m128 mr1;
	__m128 mr2;
	mat4_t r;

	// 1 - 2yy - 2zz, 2xy - 2zw, 2xz + 2yw, 0
	// 2xy + 2zw, 1 - 2xx - 2zz, 2yz - 2xw, 0
	// 2xz - 2yw, 2yz + 2xw, 1 - 2xx - 2yy, 0
	// 0, 0, 0, 1

	// -yy  xy  yw  yyy * yxw
	// -zz -zw  xz  zzz * zwx

	//  xy -xx -xw  xxx * yxw
	//  zw -zz  yz  zzz * wzy
	
	//  xz  xw -xx  xxx * zwx
	// -yw  yz -yy  yyy * wzy

	mr1 = _mm_mul_ps(myyyy, myxwz);
	mr1 = _mm_mul_ps(mr1, _mm_set_ps(0, 2, 2,-2));
	mr2 = _mm_mul_ps(mzzzz, mzwxy);
	mr2 = _mm_mul_ps(mr2, _mm_set_ps(0, 2,-2,-2));
	vec4_set_result(&r.x, _mm_add_ps(_mm_set_ps(0,0,0,1), _mm_add_ps(mr1, mr2)));

	mr1 = _mm_mul_ps(mxxxx, myxwz);
	mr1 = _mm_mul_ps(mr1, _mm_set_ps(0,-2,-2, 2));
	mr2 = _mm_mul_ps(mzzzz, mwzyx);
	mr2 = _mm_mul_ps(mr2, _mm_set_ps(0, 2,-2, 2));
	vec4_set_result(&r.y, _mm_add_ps(_mm_set_ps(0,0,1,0), _mm_add_ps(mr1, mr2)));

	mr1 = _mm_mul_ps(mxxxx, mzwxy);
	mr1 = _mm_mul_ps(mr1, _mm_set_ps(0,-2, 2, 2));
	mr2 = _mm_mul_ps(myyyy, mwzyx);
	mr2 = _mm_mul_ps(mr2, _mm_set_ps(0,-2, 2,-2));
	vec4_set_result(&r.z, _mm_add_ps(_mm_set_ps(0,1,0,0), _mm_add_ps(mr1, mr2)));

	vec4_set_result(&r.w, _mm_set_ps(1,0,0,0));
	return r;
}

mathsimd_func(mat4_t,mat4_transpose)(mat4_t m)
{
	__m128 mx = vec4_load(m.x);
	__m128 my = vec4_load(m.y);
	__m128 mz = vec4_load(m.z);
	__m128 mw = vec4_load(m.w);
	mat4_t r;
	
	_MM_TRANSPOSE4_PS(mx, my, mz, mw);
	
	vec4_set_result(&r.x, mx);
	vec4_set_result(&r.y, my);
	vec4_set_result(&r.z, mz);
	vec4_set_result(&r.w, mw);
	return r;
}

mathsimd_func(mat4_t,mat4_add)(mat4_t l, mat4_t r)
{
	mat4_t o;
	vec4_set_result(&o.x, _mm_add_ps(vec4_load(l.x), vec4_load(r.x)));
	vec4_set_result(&o.y, _mm_add_ps(vec4_load(l.y), vec4_load(r.y)));
	vec4_set_result(&o.z, _mm_add_ps(vec4_load(l.z), vec4_load(r.z)));
	vec4_set_result(&o.w, _mm_add_ps(vec4_load(l.w), vec4_load(r.w)));
	return o;
}

mathsimd_func(mat4_t,mat4_add_s)(mat4_t m, real_t s)
{
	__m128 ms = _mm_load1_ps(&s);
	mat4_t r;
	vec4_set_result(&r.x, _mm_add_ps(vec4_load(m.x), ms));
	vec4_set_result(&r.y, _mm_add_ps(vec4_load(m.y), ms));
	vec4_set_result(&r.z, _mm_add_ps(vec4_load(m.z), ms));
	vec4_set_result(&r.w, _mm_add_ps(vec4_load(m.w), ms));
	return r;
}

mathsimd_func(mat4_t,mat4_add_transpose)(mat4_t l, mat4_t r)
{
	mat4_t o;
	__m128 mrx = vec4_load(r.x);
	__m128 mry = vec4_load(r.y);
	__m128 mrz = vec4_load(r.z);
	__m128 mrw = vec4_load(r.w);
	_MM_TRANSPOSE4_PS(mrx, mry, mrz, mrw);
	vec4_set_result(&o.x, _mm_add_ps(vec4_load(l.x), mrx));
	vec4_set_result(&o.y, _mm_add_ps(vec4_load(l.y), mry));
	vec4_set_result(&o.z, _mm_add_ps(vec4_load(l.z), mrz));
	vec4_set_result(&o.w, _mm_add_ps(vec4_load(l.w), mrw));
	return o;
}

mathsimd_func(mat4_t,mat4_sub)(mat4_t l, mat4_t r)
{
	mat4_t o;
	vec4_set_result(&o.x, _mm_sub_ps(vec4_load(l.x), vec4_load(r.x)));
	vec4_set_result(&o.y, _mm_sub_ps(vec4_load(l.y), vec4_load(r.y)));
	vec4_set_result(&o.z, _mm_sub_ps(vec4_load(l.z), vec4_load(r.z)));
	vec4_set_result(&o.w, _mm_sub_ps(vec4_load(l.w), vec4_load(r.w)));
	return o;
}

mathsimd_func(mat4_t,mat4_sub_s)(mat4_t m, real_t s)
{
	__m128 ms = _mm_load1_ps(&s);
	mat4_t r;
	vec4_set_result(&r.x, _mm_sub_ps(vec4_load(m.x), ms));
	vec4_set_result(&r.y, _mm_sub_ps(vec4_load(m.y), ms));
	vec4_set_result(&r.z, _mm_sub_ps(vec4_load(m.z), ms));
	vec4_set_result(&r.w, _mm_sub_ps(vec4_load(m.w), ms));
	return r;
}

mathsimd_func(mat4_t,mat4_sub_transpose)(mat4_t l, mat4_t r)
{
	mat4_t o;
	__m128 mrx = vec4_load(r.x);
	__m128 mry = vec4_load(r.y);
	__m128 mrz = vec4_load(r.z);
	__m128 mrw = vec4_load(r.w);
	_MM_TRANSPOSE4_PS(mrx, mry, mrz, mrw);
	vec4_set_result(&o.x, _mm_sub_ps(vec4_load(l.x), mrx));
	vec4_set_result(&o.y, _mm_sub_ps(vec4_load(l.y), mry));
	vec4_set_result(&o.z, _mm_sub_ps(vec4_load(l.z), mrz));
	vec4_set_result(&o.w, _mm_sub_ps(vec4_load(l.w), mrw));
	return o;
}

mathsimd_func(mat4_t,mat4_mul)(mat4_t l, mat4_t r)
{
	__m128 mrx = vec4_load(r.x);
	__m128 mry = vec4_load(r.y);
	__m128 mrz = vec4_load(r.z);
	__m128 mrw = vec4_load(r.w);
	__m128 mlv;
	__m128 t1, t2, t3, t4;
	
	mat4_t o;
	
	mlv = vec4_load(l.x);
	t1 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(0,0,0,0));
	t2 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(1,1,1,1));
	t3 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(2,2,2,2));
	t4 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(3,3,3,3));
	vec4_set_result(&o.x, _mm_add_ps(
		_mm_add_ps(_mm_mul_ps(t1, mrx), _mm_mul_ps(t2, mry)),
		_mm_add_ps(_mm_mul_ps(t3, mrz), _mm_mul_ps(t4, mrw))));
	
	mlv = vec4_load(l.y);
	t1 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(0,0,0,0));
	t2 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(1,1,1,1));
	t3 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(2,2,2,2));
	t4 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(3,3,3,3));
	vec4_set_result(&o.y, _mm_add_ps(
		_mm_add_ps(_mm_mul_ps(t1, mrx), _mm_mul_ps(t2, mry)),
		_mm_add_ps(_mm_mul_ps(t3, mrz), _mm_mul_ps(t4, mrw))));
	
	mlv = vec4_load(l.z);
	t1 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(0,0,0,0));
	t2 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(1,1,1,1));
	t3 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(2,2,2,2));
	t4 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(3,3,3,3));
	vec4_set_result(&o.z, _mm_add_ps(
		_mm_add_ps(_mm_mul_ps(t1, mrx), _mm_mul_ps(t2, mry)),
		_mm_add_ps(_mm_mul_ps(t3, mrz), _mm_mul_ps(t4, mrw))));
	
	mlv = vec4_load(l.w);
	t1 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(0,0,0,0));
	t2 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(1,1,1,1));
	t3 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(2,2,2,2));
	t4 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(3,3,3,3));
	vec4_set_result(&o.w, _mm_add_ps(
		_mm_add_ps(_mm_mul_ps(t1, mrx), _mm_mul_ps(t2, mry)),
		_mm_add_ps(_mm_mul_ps(t3, mrz), _mm_mul_ps(t4, mrw))));
	return o;
}

mathsimd_func(mat4_t, mat4_mul_s)(mat4_t m, real_t s)
{
	__m128 ms = _mm_load1_ps(&s);
	mat4_t r;
	vec4_set_result(&r.x, _mm_mul_ps(vec4_load(m.x), ms));
	vec4_set_result(&r.y, _mm_mul_ps(vec4_load(m.y), ms));
	vec4_set_result(&r.z, _mm_mul_ps(vec4_load(m.z), ms));
	vec4_set_result(&r.w, _mm_mul_ps(vec4_load(m.w), ms));
	return r;
}

mathsimd_func(mat4_t,mat4_mul_transpose)(mat4_t l, mat4_t r)
{
	__m128 mrx = vec4_load(r.x);
	__m128 mry = vec4_load(r.y);
	__m128 mrz = vec4_load(r.z);
	__m128 mrw = vec4_load(r.w);
	__m128 mlv;
	__m128 t1, t2, t3, t4;
	
	mat4_t o;
	_MM_TRANSPOSE4_PS(mrx, mry, mrz, mrw);
	
	mlv = vec4_load(l.x);
	t1 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(0,0,0,0));
	t2 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(1,1,1,1));
	t3 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(2,2,2,2));
	t4 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(3,3,3,3));
	vec4_set_result(&o.x, _mm_add_ps(
		_mm_add_ps(_mm_mul_ps(t1, mrx), _mm_mul_ps(t2, mry)),
		_mm_add_ps(_mm_mul_ps(t3, mrz), _mm_mul_ps(t4, mrw))));
	
	mlv = vec4_load(l.y);
	t1 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(0,0,0,0));
	t2 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(1,1,1,1));
	t3 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(2,2,2,2));
	t4 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(3,3,3,3));
	vec4_set_result(&o.y, _mm_add_ps(
		_mm_add_ps(_mm_mul_ps(t1, mrx), _mm_mul_ps(t2, mry)),
		_mm_add_ps(_mm_mul_ps(t3, mrz), _mm_mul_ps(t4, mrw))));
	
	mlv = vec4_load(l.z);
	t1 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(0,0,0,0));
	t2 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(1,1,1,1));
	t3 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(2,2,2,2));
	t4 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(3,3,3,3));
	vec4_set_result(&o.z, _mm_add_ps(
		_mm_add_ps(_mm_mul_ps(t1, mrx), _mm_mul_ps(t2, mry)),
		_mm_add_ps(_mm_mul_ps(t3, mrz), _mm_mul_ps(t4, mrw))));
	
	mlv = vec4_load(l.w);
	t1 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(0,0,0,0));
	t2 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(1,1,1,1));
	t3 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(2,2,2,2));
	t4 = _mm_shuffle_ps(mlv, mlv, _MM_SHUFFLE(3,3,3,3));
	vec4_set_result(&o.w, _mm_add_ps(
		_mm_add_ps(_mm_mul_ps(t1, mrx), _mm_mul_ps(t2, mry)),
		_mm_add_ps(_mm_mul_ps(t3, mrz), _mm_mul_ps(t4, mrw))));
	return o;
}

#if MATHUTIL_DETECT_CPU
int mathutil_sse_implements()
{
	if(!CPUID_SSE()) return 0;

	r_sin = r_sin_sse;
	r_cos = r_cos_sse;
	// r_exp = r_exp_sse;
	r_log = r_log_sse;
	// r_pow = r_pow_sse;

	vec4_abs = vec4_abs_sse;
	vec4_sgn = vec4_sgn_sse;
	vec4_invert = vec4_invert_sse;
	vec4_length = vec4_length_sse;
	vec4_scale = vec4_scale_sse;
	vec4_dot = vec4_dot_sse;
	vec4_add = vec4_add_sse;
	vec4_sub = vec4_sub_sse;
	vec4_mul = vec4_mul_sse;
	vec4_div = vec4_div_sse;
	vec4_min = vec4_min_sse;
	vec4_max = vec4_max_sse;
	vec4_clamp = vec4_clamp_sse;
	vec4_mul_mat4_transpose = vec4_mul_mat4_transpose_sse;
	vec4_lerp = vec4_lerp_sse;
	vec4_slerp = vec4_slerp_sse;
	vec4_cross3 = vec4_cross3_sse;
	vec4_rot_quat = vec4_rot_quat_sse;

	vec4_normalize = vec4_normalize_sse;
	vec4_mul_mat4 = vec4_mul_mat4_sse;

	quat_mul = quat_mul_sse;
	quat_add_vec = quat_add_vec_sse;

	mat4_from_quat = mat4_from_quat_sse;
	mat4_transpose = mat4_transpose_sse;
	mat4_from_quat_transpose = mat4_from_quat_transpose_sse;
	mat4_add = mat4_add_sse;
	mat4_add_s = mat4_add_s_sse;
	mat4_add_transpose = mat4_add_transpose_sse;
	mat4_sub = mat4_sub_sse;
	mat4_sub_s = mat4_sub_s_sse;
	mat4_sub_transpose = mat4_sub_transpose_sse;
	mat4_mul = mat4_mul_sse;
	mat4_mul_transpose = mat4_mul_transpose_sse;
	return 1;
}
#endif // MATHUTIL_DETECT_CPU

#else // MATHUTIL_USE_DOUBLE

#if defined(_MSC_VER)

#pragma comment(linker,"/alternatename:r_rnd_sse=r_rnd_ref")
#pragma comment(linker,"/alternatename:r_sin_sse=r_sin_ref")
#pragma comment(linker,"/alternatename:r_cos_sse=r_cos_ref")
#pragma comment(linker,"/alternatename:r_tan_sse=r_tan_ref")
#pragma comment(linker,"/alternatename:r_abs_sse=r_abs_ref")
#pragma comment(linker,"/alternatename:r_sgn_sse=r_sgn_ref")
#pragma comment(linker,"/alternatename:r_sqr_sse=r_sqr_ref")
#pragma comment(linker,"/alternatename:r_floor_sse=r_floor_ref")
#pragma comment(linker,"/alternatename:r_ceil_sse=r_ceil_ref")
#pragma comment(linker,"/alternatename:r_atan_sse=r_atan_ref")
#pragma comment(linker,"/alternatename:r_exp_sse=r_exp_ref")
#pragma comment(linker,"/alternatename:r_log_sse=r_log_ref")
#pragma comment(linker,"/alternatename:r_pow_sse=r_pow_ref")
#pragma comment(linker,"/alternatename:r_mod_sse=r_mod_ref")
#pragma comment(linker,"/alternatename:r_max_sse=r_max_ref")
#pragma comment(linker,"/alternatename:r_min_sse=r_min_ref")
#pragma comment(linker,"/alternatename:r_atan2_sse=r_atan2_ref")
#pragma comment(linker,"/alternatename:r_clamp_sse=r_clamp_ref")
#pragma comment(linker,"/alternatename:r_lerp_sse=r_lerp_ref")
#pragma comment(linker,"/alternatename:r_hermite_sse=r_hermite_ref")
#pragma comment(linker,"/alternatename:r_slerp_sse=r_slerp_ref")
#pragma comment(linker,"/alternatename:vec4_sse=vec4_ref")
#pragma comment(linker,"/alternatename:vec4_flushcomp_sse=vec4_flushcomp_ref")
#pragma comment(linker,"/alternatename:vec4_abs_sse=vec4_abs_ref")
#pragma comment(linker,"/alternatename:vec4_sgn_sse=vec4_sgn_ref")
#pragma comment(linker,"/alternatename:vec4_invert_sse=vec4_invert_ref")
#pragma comment(linker,"/alternatename:vec4_length_sse=vec4_length_ref")
#pragma comment(linker,"/alternatename:vec4_normalize_sse=vec4_normalize_ref")
#pragma comment(linker,"/alternatename:vec4_scale_sse=vec4_scale_ref")
#pragma comment(linker,"/alternatename:vec4_pow_sse=vec4_pow_ref")
#pragma comment(linker,"/alternatename:vec4_dot_sse=vec4_dot_ref")
#pragma comment(linker,"/alternatename:vec4_add_sse=vec4_add_ref")
#pragma comment(linker,"/alternatename:vec4_sub_sse=vec4_sub_ref")
#pragma comment(linker,"/alternatename:vec4_mul_sse=vec4_mul_ref")
#pragma comment(linker,"/alternatename:vec4_div_sse=vec4_div_ref")
#pragma comment(linker,"/alternatename:vec4_min_sse=vec4_min_ref")
#pragma comment(linker,"/alternatename:vec4_max_sse=vec4_max_ref")
#pragma comment(linker,"/alternatename:vec4_cross3_sse=vec4_cross3_ref")
#pragma comment(linker,"/alternatename:vec4_clamp_sse=vec4_clamp_ref")
#pragma comment(linker,"/alternatename:vec4_rot_quat_sse=vec4_rot_quat_ref")
#pragma comment(linker,"/alternatename:vec4_mul_mat4_sse=vec4_mul_mat4_ref")
#pragma comment(linker,"/alternatename:vec4_mul_mat4_transpose_sse=vec4_mul_mat4_transpose_ref")
#pragma comment(linker,"/alternatename:vec4_lerp_sse=vec4_lerp_ref")
#pragma comment(linker,"/alternatename:vec4_slerp_sse=vec4_slerp_ref")
#pragma comment(linker,"/alternatename:quat_sse=quat_ref")
#pragma comment(linker,"/alternatename:quat_flushcomp_sse=quat_flushcomp_ref")
#pragma comment(linker,"/alternatename:quat_rot_axis_sse=quat_rot_axis_ref")
#pragma comment(linker,"/alternatename:quat_mul_sse=quat_mul_ref")
#pragma comment(linker,"/alternatename:quat_add_vec_sse=quat_add_vec_ref")
#pragma comment(linker,"/alternatename:mat4_sse=mat4_ref")
#pragma comment(linker,"/alternatename:mat4_flushcomp_sse=mat4_flushcomp_ref")
#pragma comment(linker,"/alternatename:mat4_rot_x_sse=mat4_rot_x_ref")
#pragma comment(linker,"/alternatename:mat4_rot_y_sse=mat4_rot_y_ref")
#pragma comment(linker,"/alternatename:mat4_rot_z_sse=mat4_rot_z_ref")
#pragma comment(linker,"/alternatename:mat4_rot_axis_sse=mat4_rot_axis_ref")
#pragma comment(linker,"/alternatename:mat4_rot_euler_sse=mat4_rot_euler_ref")
#pragma comment(linker,"/alternatename:mat4_from_quat_sse=mat4_from_quat_ref")
#pragma comment(linker,"/alternatename:mat4_from_quat_transpose_sse=mat4_from_quat_transpose_ref")
#pragma comment(linker,"/alternatename:mat4_transpose_sse=mat4_transpose_ref")
#pragma comment(linker,"/alternatename:mat4_add_sse=mat4_add_ref")
#pragma comment(linker,"/alternatename:mat4_add_s_sse=mat4_add_s_ref")
#pragma comment(linker,"/alternatename:mat4_add_transpose_sse=mat4_add_transpose_ref")
#pragma comment(linker,"/alternatename:mat4_sub_sse=mat4_sub_ref")
#pragma comment(linker,"/alternatename:mat4_sub_s_sse=mat4_sub_s_ref")
#pragma comment(linker,"/alternatename:mat4_sub_transpose_sse=mat4_sub_transpose_ref")
#pragma comment(linker,"/alternatename:mat4_mul_sse=mat4_mul_ref")
#pragma comment(linker,"/alternatename:mat4_mul_s_sse=mat4_mul_s_ref")
#pragma comment(linker,"/alternatename:mat4_mul_transpose_sse=mat4_mul_transpose_ref")

#else

#define math_func(r,n,arg) r n arg __attribute__ ((weak, alias ("" # n # "_ref")));
#include"mathutil_funclist.h"
#undef math_func

#endif

#if MATHUTIL_DETECT_CPU
int mathutil_sse_implements()
{
	if(!CPUID_SSE()) return 0;

	// Yeah

	return 1;
}
#endif // MATHUTIL_DETECT_CPU

#endif // MATHUTIL_USE_DOUBLE
