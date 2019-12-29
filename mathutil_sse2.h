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

// mathutil_sse2
// Author: 0xAA55
// The implemention of mathutil API which uses SSE2 intrinsics.

#include"mathutil_conf.h"

#if HAVE_SSE2 || __INTELLISENSE__

#ifndef sse2_func
  #if __INTELLISENSE__
    #define sse2_func(r,n) r n
    #include<xmmintrin.h>
    #include<immintrin.h>
    #include<emmintrin.h>
  #endif // __INTELLISENSE__
#endif

#include"mathutil_sse_common.h"

#if !MATHUTIL_USE_DOUBLE

#if !vec4_abs_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_abs)(vec4_t v)
{
	vec4_t r;
	vec4_set_iresult(&r, _mm_and_si128(vec4_loadi(v), _mm_set1_epi32(0x7fffffff)));
	return r;
}
#define vec4_abs_implemented 1
#endif

#if !vec4_sgn_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_sgn)(vec4_t v)
{
	vec4_t r;
	vec4_set_iresult(&r, _mm_or_si128(_mm_and_si128(vec4_loadi(v), _mm_set1_epi32(0x80000000)), _mm_set1_epi32(0x3F800000)));
	return r;
}
#define vec4_sgn_implemented 1
#endif

#if !vec4_invert_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_invert)(vec4_t v)
{
	vec4_t r;
	vec4_set_iresult(&r, _mm_xor_si128(vec4_loadi(v), _mm_set1_epi32(0x80000000)));
	return r;
}
#define vec4_invert_implemented 1
#endif

#else // MATHUTIL_USE_DOUBLE

#if !vec4_abs_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_abs)(vec4_t v)
{
	vec4_t r;
	static const ALIGNED_(16) int64_t ins[2] = { 0x7fffffffffffffffll, 0x7fffffffffffffffll };
	__m128i mns = _mm_load_si128((__m128i*) & ins);
	_mm_store_si128_a((__m128i*) & r.x, _mm_and_si128(vec4_loadi_xy(v), mns));
	_mm_store_si128_a((__m128i*) & r.z, _mm_and_si128(vec4_loadi_zw(v), mns));
	return r;
}
#define vec4_abs_implemented 1
#endif

#if !vec4_sgn_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_sgn)(vec4_t v)
{
	vec4_t r;
	static const ALIGNED_(16) int64_t isgnb[2] = { 0x8000000000000000ll, 0x8000000000000000ll };
	static const ALIGNED_(16) int64_t i1d[2] = { 0x3ff0000000000000ll, 0x3ff0000000000000ll };
	__m128i msgnb = _mm_load_si128((void*) & isgnb);
	__m128i m1d = _mm_load_si128((void*) & i1d);
	_mm_store_si128_a((__m128i*) & r.x, _mm_or_si128(_mm_and_si128(vec4_loadi_xy(v), msgnb), m1d));
	_mm_store_si128_a((__m128i*) & r.z, _mm_or_si128(_mm_and_si128(vec4_loadi_zw(v), msgnb), m1d));
	return r;
}
#define vec4_sgn_implemented 1
#endif

#if !vec4_invert_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_invert)(vec4_t v)
{
	vec4_t r;
	static const ALIGNED_(16) int64_t isgnb[2] = { 0x8000000000000000ll, 0x8000000000000000ll };
	__m128i msgnb = _mm_load_si128((__m128i*) & isgnb);
	_mm_store_si128_a((__m128i*) & r.x, _mm_xor_si128(vec4_loadi_xy(v), msgnb));
	_mm_store_si128_a((__m128i*) & r.z, _mm_xor_si128(vec4_loadi_zw(v), msgnb));
	return r;
}
#define vec4_invert_implemented 1
#endif

#if !vec4_length_implemented || MATHUTIL_DETECT_CPU
sse2_func(real_t, vec4_length)(vec4_t v)
{
	__m128d mvxy = vec4_load_xy(v);
	__m128d mvzw = vec4_load_zw(v);
	__m128d mxaddz_yaddw = _mm_add_pd(_mm_mul_pd(mvxy, mvxy), _mm_mul_pd(mvzw, mvzw));
	__m128d lsq = _mm_add_sd(mxaddz_yaddw, _mm_unpackhi_pd(mxaddz_yaddw, mxaddz_yaddw));

	real_t r;
	_mm_store_sd(&r, _mm_sqrt_pd(lsq));
	_mm_sfence();
	_compiler_barrier;
	return r;
}
#define vec4_length_implemented 1
#endif

#if !vec4_normalize_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_normalize)(vec4_t v)
{
	vec4_t r;
	__m128d mvxy = vec4_load_xy(v);
	__m128d mvzw = vec4_load_zw(v);
	__m128d mxaddz_yaddw = _mm_add_pd(_mm_mul_pd(mvxy, mvxy), _mm_mul_pd(mvzw, mvzw));
	__m128d lsq = _mm_add_sd(mxaddz_yaddw, _mm_unpackhi_pd(mxaddz_yaddw, mxaddz_yaddw));
	__m128d length = _mm_sqrt_pd(lsq);
	vec4_set_result_xy(&r, _mm_div_pd(mvxy, length));
	vec4_set_result_zw(&r, _mm_div_pd(mvzw, length));
	return r;
}
#define vec4_normalize_implemented 1
#endif

#if !vec4_scale_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_scale)(vec4_t v, real_t s)
{
	vec4_t r;
	__m128d ms = _mm_load1_pd(&s);
	vec4_set_result_xy(&r, _mm_mul_pd(vec4_load_xy(v), ms));
	vec4_set_result_zw(&r, _mm_mul_pd(vec4_load_zw(v), ms));
	return r;
}
#define vec4_scale_implemented 1
#endif

#if !vec4_clamp_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_clamp)(vec4_t v, real_t min_, real_t max_)
{
	vec4_t r;
	__m128d ml, mh;

	ml = _mm_load1_pd(&min_);
	mh = _mm_load1_pd(&max_);

	vec4_set_result_xy(&r, _mm_min_pd(_mm_max_pd(vec4_load_xy(v), ml), mh));
	vec4_set_result_zw(&r, _mm_min_pd(_mm_max_pd(vec4_load_zw(v), ml), mh));
	return r;
}
#define vec4_clamp_implemented 1
#endif

#if !vec4_dot_implemented || MATHUTIL_DETECT_CPU
sse2_func(real_t, vec4_dot)(vec4_t v1, vec4_t v2)
{
	__m128d mxaddz_yaddw = _mm_add_pd(
		_mm_mul_pd(vec4_load_xy(v1), vec4_load_xy(v2)),
		_mm_mul_pd(vec4_load_zw(v1), vec4_load_zw(v2)));
	real_t r;
	_mm_store_sd(&r, _mm_add_sd(mxaddz_yaddw, _mm_unpackhi_pd(mxaddz_yaddw, mxaddz_yaddw)));
	_mm_sfence();
	_compiler_barrier;
	return r;
}
#define vec4_dot_implemented 1
#endif

#if !vec4_add_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_add)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	__m128d mxy = _mm_add_pd(vec4_load_xy(v1), vec4_load_xy(v2));
	__m128d mzw = _mm_add_pd(vec4_load_zw(v1), vec4_load_zw(v2));
	vec4_set_result_xy(&r, mxy);
	vec4_set_result_zw(&r, mzw);
	return r;
}
#define vec4_add_implemented 1
#endif

#if !vec4_sub_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_sub)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	__m128d mxy = _mm_sub_pd(vec4_load_xy(v1), vec4_load_xy(v2));
	__m128d mzw = _mm_sub_pd(vec4_load_zw(v1), vec4_load_zw(v2));
	vec4_set_result_xy(&r, mxy);
	vec4_set_result_zw(&r, mzw);
	return r;
}
#define vec4_sub_implemented 1
#endif

#if !vec4_mul_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_mul)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	__m128d mxy = _mm_mul_pd(vec4_load_xy(v1), vec4_load_xy(v2));
	__m128d mzw = _mm_mul_pd(vec4_load_zw(v1), vec4_load_zw(v2));
	vec4_set_result_xy(&r, mxy);
	vec4_set_result_zw(&r, mzw);
	return r;
}
#define vec4_mul_implemented 1
#endif

#if !vec4_div_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_div)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	__m128d mxy = _mm_div_pd(vec4_load_xy(v1), vec4_load_xy(v2));
	__m128d mzw = _mm_div_pd(vec4_load_zw(v1), vec4_load_zw(v2));
	vec4_set_result_xy(&r, mxy);
	vec4_set_result_zw(&r, mzw);
	return r;
}
#define vec4_div_implemented 1
#endif

#if !vec4_min_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_min)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	vec4_set_result_xy(&r, _mm_min_pd(vec4_load_xy(v1), vec4_load_xy(v2)));
	vec4_set_result_zw(&r, _mm_min_pd(vec4_load_zw(v1), vec4_load_zw(v2)));
	return r;
}
#define vec4_min_implemented 1
#endif

#if !vec4_max_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_max)(vec4_t v1, vec4_t v2)
{
	vec4_t r;
	vec4_set_result_xy(&r, _mm_max_pd(vec4_load_xy(v1), vec4_load_xy(v2)));
	vec4_set_result_zw(&r, _mm_max_pd(vec4_load_zw(v1), vec4_load_zw(v2)));
	return r;
}
#define vec4_max_implemented 1
#endif

#if !vec4_mul_mat4_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_mul_mat4)(vec4_t v, mat4_t m)
{
	vec4_t r;

	__m128d mvxy = vec4_load_xy(v);
	__m128d mvzw = vec4_load_zw(v);

	vec4_set_result_xy(&r, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), vec4_load_xy(m.x)), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), vec4_load_xy(m.y))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), vec4_load_xy(m.z)), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), vec4_load_xy(m.w)))));
	vec4_set_result_zw(&r, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), vec4_load_zw(m.x)), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), vec4_load_zw(m.y))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), vec4_load_zw(m.z)), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), vec4_load_zw(m.w)))));

	return r;
}
#define vec4_mul_mat4_implemented 1
#endif

#if !vec4_mul_mat4_transpose_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_mul_mat4_transpose)(vec4_t v, mat4_t m)
{
	vec4_t r;

	__m128d mvxy = vec4_load_xy(v);
	__m128d mvzw = vec4_load_zw(v);

	__m128d mx_xy = vec4_load_xy(m.x);
	__m128d my_xy = vec4_load_xy(m.y);
	__m128d mz_xy = vec4_load_xy(m.z);
	__m128d mw_xy = vec4_load_xy(m.w);

	__m128d mx_zw = vec4_load_zw(m.x);
	__m128d my_zw = vec4_load_zw(m.y);
	__m128d mz_zw = vec4_load_zw(m.z);
	__m128d mw_zw = vec4_load_zw(m.w);

	vec4_set_result_xy(&r, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mx_xy, my_xy)), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mx_xy, my_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mx_zw, my_zw)), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mx_zw, my_zw)))));
	vec4_set_result_zw(&r, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mz_xy, mw_xy)), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mz_xy, mw_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mz_zw, mw_zw)), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mz_zw, mw_zw)))));

	return r;
}
#define vec4_mul_mat4_transpose_implemented 1
#endif

#if !vec4_lerp_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_lerp)(vec4_t v1, vec4_t v2, real_t s)
{
	__m128d maxy = vec4_load_xy(v1);
	__m128d mazw = vec4_load_zw(v1);
	__m128d mbxy = vec4_load_xy(v2);
	__m128d mbzw = vec4_load_zw(v2);
	__m128d ms = _mm_load1_pd(&s);
	vec4_t r;

	vec4_set_result_xy(&r, _mm_add_pd(maxy, _mm_mul_pd(_mm_sub_pd(mbxy, maxy), ms)));
	vec4_set_result_zw(&r, _mm_add_pd(mazw, _mm_mul_pd(_mm_sub_pd(mbzw, mazw), ms)));
	return r;
}
#define vec4_lerp_implemented 1
#endif

#if !vec4_slerp_implemented || MATHUTIL_DETECT_CPU
sse2_func(vec4_t, vec4_slerp)(vec4_t v1, vec4_t v2, real_t s)
{
	__m128d maxy = vec4_load_xy(v1);
	__m128d mazw = vec4_load_zw(v1);
	__m128d mbxy = vec4_load_xy(v2);
	__m128d mbzw = vec4_load_zw(v2);
	__m128d ms = _mm_load1_pd(&s);
	vec4_t r;

	ms = _mm_mul_pd(_mm_mul_pd(ms, ms), _mm_sub_pd(_mm_set1_pd(3), _mm_mul_pd(_mm_set1_pd(2), ms)));
	ms = _mm_max_pd(_mm_setzero_pd(), _mm_min_pd(_mm_set1_pd(1), ms));
	vec4_set_result_xy(&r, _mm_add_pd(maxy, _mm_mul_pd(_mm_sub_pd(mbxy, maxy), ms)));
	vec4_set_result_zw(&r, _mm_add_pd(mazw, _mm_mul_pd(_mm_sub_pd(mbzw, mazw), ms)));
	return r;
}
#define vec4_slerp_implemented 1
#endif

#if !quat_mul_implemented || MATHUTIL_DETECT_CPU
sse2_func(quat_t, quat_mul)(quat_t q1, quat_t q2)
{
	__m128d mq1xy = quat_load_xy(q1);
	__m128d mq1zw = quat_load_zw(q1);
	__m128d mq2xy = quat_load_xy(q2);
	__m128d mq2zw = quat_load_zw(q2);
	__m128d mpn = _mm_set_pd(-1, 1);
	__m128d mnp = _mm_set_pd(1, -1);
	__m128d mqxy;
	__m128d mqzw;
	quat_t r;

	mqxy = _mm_mul_pd(_mm_unpackhi_pd(mq1zw, mq1zw), mq2xy);
	mqxy = _mm_add_pd(mqxy, _mm_mul_pd(_mm_unpackhi_pd(mq1xy, mq1xy), mq2zw));
	mqxy = _mm_add_pd(mqxy, _mm_mul_pd(_mm_mul_pd(mpn, _mm_unpacklo_pd(mq1xy, mq1xy)), _mm_shuffle_pd(mq2zw, mq2zw, _MM_SHUFFLE2(0, 1))));
	mqxy = _mm_add_pd(mqxy, _mm_mul_pd(_mm_mul_pd(mnp, _mm_unpacklo_pd(mq1zw, mq1zw)), _mm_shuffle_pd(mq2xy, mq2xy, _MM_SHUFFLE2(0, 1))));

	mqzw = _mm_mul_pd(_mm_unpackhi_pd(mq1zw, mq1zw), mq2zw);
	mqzw = _mm_sub_pd(mqzw, _mm_mul_pd(_mm_unpackhi_pd(mq1xy, mq1xy), mq2xy));
	mqzw = _mm_add_pd(mqzw, _mm_mul_pd(_mm_mul_pd(mpn, _mm_unpacklo_pd(mq1xy, mq1xy)), _mm_shuffle_pd(mq2xy, mq2xy, _MM_SHUFFLE2(0, 1))));
	mqzw = _mm_add_pd(mqzw, _mm_mul_pd(_mm_mul_pd(mpn, _mm_unpacklo_pd(mq1zw, mq1zw)), _mm_shuffle_pd(mq2zw, mq2zw, _MM_SHUFFLE2(0, 1))));

	vec4_set_result_xy(&r, mqxy);
	vec4_set_result_zw(&r, mqzw);
	return r;
}
#define quat_mul_implemented 1
#endif

#if !quat_add_vec_implemented || MATHUTIL_DETECT_CPU
sse2_func(quat_t, quat_add_vec)(quat_t q, vec4_t v, real_t s)
{
	__m128d mq1xy = vec4_load_xy(v);
	__m128d mq1zw = vec4_load_zw(v);
	__m128d mq2xy = quat_load_xy(q);
	__m128d mq2zw = quat_load_zw(q);
	__m128d mpn = _mm_set_pd(-1, 1);
	__m128d mnp = _mm_set_pd(1, -1);
	__m128d mqxy;
	__m128d mqzw;
	__m128d ms = _mm_load1_pd(&s);
	__m128d mhf = _mm_set1_pd(0.5);
	quat_t r;

	mq1xy = _mm_mul_pd(mq1xy, ms);
	mq1zw = _mm_mul_sd(mq1zw, ms);

	mqxy = _mm_mul_pd(_mm_unpackhi_pd(mq1xy, mq1xy), mq2zw);
	mqxy = _mm_add_pd(mqxy, _mm_mul_pd(_mm_mul_pd(mpn, _mm_unpacklo_pd(mq1xy, mq1xy)), _mm_shuffle_pd(mq2zw, mq2zw, _MM_SHUFFLE2(0, 1))));
	mqxy = _mm_add_pd(mqxy, _mm_mul_pd(_mm_mul_pd(mnp, _mm_unpacklo_pd(mq1zw, mq1zw)), _mm_shuffle_pd(mq2xy, mq2xy, _MM_SHUFFLE2(0, 1))));

	mqzw = _mm_mul_pd(_mm_mul_pd(mpn, _mm_unpacklo_pd(mq1xy, mq1xy)), _mm_shuffle_pd(mq2xy, mq2xy, _MM_SHUFFLE2(0, 1)));
	mqzw = _mm_add_pd(mqzw, _mm_mul_pd(_mm_mul_pd(mpn, _mm_unpacklo_pd(mq1zw, mq1zw)), _mm_shuffle_pd(mq2zw, mq2zw, _MM_SHUFFLE2(0, 1))));
	mqzw = _mm_sub_pd(mqzw, _mm_mul_pd(_mm_unpackhi_pd(mq1xy, mq1xy), mq2xy));

	vec4_set_result_xy(&r, _mm_add_pd(_mm_mul_pd(mqxy, mhf), mq2xy));
	vec4_set_result_zw(&r, _mm_add_pd(_mm_mul_pd(mqzw, mhf), mq2zw));
	return r;
}
#define quat_add_vec_implemented 1
#endif

#if !mat4_transpose_implemented || MATHUTIL_DETECT_CPU
sse2_func(mat4_t, mat4_transpose)(mat4_t m)
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

	mx_xy = vec4_load_xy(m.x);
	my_xy = vec4_load_xy(m.y);
	mz_xy = vec4_load_xy(m.z);
	mw_xy = vec4_load_xy(m.w);

	t = _mm_unpacklo_pd(mx_xy, my_xy);
	vec4_set_result_xy(&r.x, t);

	t = _mm_unpacklo_pd(mz_xy, mw_xy);
	vec4_set_result_zw(&r.x, t);

	t = _mm_unpackhi_pd(mx_xy, my_xy);
	vec4_set_result_xy(&r.y, t);

	t = _mm_unpackhi_pd(mz_xy, mw_xy);
	vec4_set_result_zw(&r.y, t);

	mx_zw = vec4_load_zw(m.x);
	my_zw = vec4_load_zw(m.y);
	mz_zw = vec4_load_zw(m.z);
	mw_zw = vec4_load_zw(m.w);

	t = _mm_unpacklo_pd(mx_zw, my_zw);
	vec4_set_result_xy(&r.z, t);

	t = _mm_unpacklo_pd(mz_zw, mw_zw);
	vec4_set_result_zw(&r.z, t);

	t = _mm_unpackhi_pd(mx_zw, my_zw);
	vec4_set_result_xy(&r.w, t);

	t = _mm_unpackhi_pd(mz_zw, mw_zw);
	vec4_set_result_zw(&r.w, t);

	return r;
}
#define mat4_transpose_implemented 1
#endif

#if !mat4_add_implemented || MATHUTIL_DETECT_CPU
sse2_func(mat4_t, mat4_add)(mat4_t l, mat4_t r)
{
	mat4_t o;
	vec4_set_result_xy(&o.x, _mm_add_pd(vec4_load_xy(l.x), vec4_load_xy(r.x)));
	vec4_set_result_xy(&o.y, _mm_add_pd(vec4_load_xy(l.y), vec4_load_xy(r.y)));
	vec4_set_result_xy(&o.z, _mm_add_pd(vec4_load_xy(l.z), vec4_load_xy(r.z)));
	vec4_set_result_xy(&o.w, _mm_add_pd(vec4_load_xy(l.w), vec4_load_xy(r.w)));
	vec4_set_result_zw(&o.x, _mm_add_pd(vec4_load_zw(l.x), vec4_load_zw(r.x)));
	vec4_set_result_zw(&o.y, _mm_add_pd(vec4_load_zw(l.y), vec4_load_zw(r.y)));
	vec4_set_result_zw(&o.z, _mm_add_pd(vec4_load_zw(l.z), vec4_load_zw(r.z)));
	vec4_set_result_zw(&o.w, _mm_add_pd(vec4_load_zw(l.w), vec4_load_zw(r.w)));
	return o;
}
#define mat4_add_implemented 1
#endif

#if !mat4_add_s_implemented || MATHUTIL_DETECT_CPU
sse2_func(mat4_t, mat4_add_s)(mat4_t m, real_t s)
{
	__m128d ms = _mm_load1_pd(&s);
	mat4_t r;
	vec4_set_result_xy(&r.x, _mm_add_pd(vec4_load_xy(m.x), ms));
	vec4_set_result_xy(&r.y, _mm_add_pd(vec4_load_xy(m.y), ms));
	vec4_set_result_xy(&r.z, _mm_add_pd(vec4_load_xy(m.z), ms));
	vec4_set_result_xy(&r.w, _mm_add_pd(vec4_load_xy(m.w), ms));
	vec4_set_result_zw(&r.x, _mm_add_pd(vec4_load_zw(m.x), ms));
	vec4_set_result_zw(&r.y, _mm_add_pd(vec4_load_zw(m.y), ms));
	vec4_set_result_zw(&r.z, _mm_add_pd(vec4_load_zw(m.z), ms));
	vec4_set_result_zw(&r.w, _mm_add_pd(vec4_load_zw(m.w), ms));
	return r;
}
#define mat4_add_s_implemented 1
#endif

#if !mat4_add_transpose_implemented || MATHUTIL_DETECT_CPU
sse2_func(mat4_t, mat4_add_transpose)(mat4_t l, mat4_t r)
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

	mx_xy = vec4_load_xy(r.x);
	my_xy = vec4_load_xy(r.y);
	mz_xy = vec4_load_xy(r.z);
	mw_xy = vec4_load_xy(r.w);

	t = _mm_unpacklo_pd(mx_xy, my_xy);
	vec4_set_result_xy(&o.x, _mm_add_pd(vec4_load_xy(l.x), t));

	t = _mm_unpacklo_pd(mz_xy, mw_xy);
	vec4_set_result_zw(&o.x, _mm_add_pd(vec4_load_zw(l.x), t));

	t = _mm_unpackhi_pd(mx_xy, my_xy);
	vec4_set_result_xy(&o.y, _mm_add_pd(vec4_load_xy(l.y), t));

	t = _mm_unpackhi_pd(mz_xy, mw_xy);
	vec4_set_result_zw(&o.y, _mm_add_pd(vec4_load_zw(l.y), t));

	mx_zw = vec4_load_zw(r.x);
	my_zw = vec4_load_zw(r.y);
	mz_zw = vec4_load_zw(r.z);
	mw_zw = vec4_load_zw(r.w);

	t = _mm_unpacklo_pd(mx_zw, my_zw);
	vec4_set_result_xy(&o.z, _mm_add_pd(vec4_load_xy(l.z), t));

	t = _mm_unpacklo_pd(mz_zw, mw_zw);
	vec4_set_result_zw(&o.z, _mm_add_pd(vec4_load_zw(l.z), t));

	t = _mm_unpackhi_pd(mx_zw, my_zw);
	vec4_set_result_xy(&o.w, _mm_add_pd(vec4_load_xy(l.w), t));

	t = _mm_unpackhi_pd(mz_zw, mw_zw);
	vec4_set_result_zw(&o.w, _mm_add_pd(vec4_load_zw(l.w), t));

	return o;
}
#define mat4_add_transpose_implemented 1
#endif

#if !mat4_sub_implemented || MATHUTIL_DETECT_CPU
sse2_func(mat4_t, mat4_sub)(mat4_t l, mat4_t r)
{
	mat4_t o;
	vec4_set_result_xy(&o.x, _mm_sub_pd(vec4_load_xy(l.x), vec4_load_xy(r.x)));
	vec4_set_result_xy(&o.y, _mm_sub_pd(vec4_load_xy(l.y), vec4_load_xy(r.y)));
	vec4_set_result_xy(&o.z, _mm_sub_pd(vec4_load_xy(l.z), vec4_load_xy(r.z)));
	vec4_set_result_xy(&o.w, _mm_sub_pd(vec4_load_xy(l.w), vec4_load_xy(r.w)));
	vec4_set_result_zw(&o.x, _mm_sub_pd(vec4_load_zw(l.x), vec4_load_zw(r.x)));
	vec4_set_result_zw(&o.y, _mm_sub_pd(vec4_load_zw(l.y), vec4_load_zw(r.y)));
	vec4_set_result_zw(&o.z, _mm_sub_pd(vec4_load_zw(l.z), vec4_load_zw(r.z)));
	vec4_set_result_zw(&o.w, _mm_sub_pd(vec4_load_zw(l.w), vec4_load_zw(r.w)));
	return o;
}
#define mat4_sub_implemented 1
#endif

#if !mat4_sub_s_implemented || MATHUTIL_DETECT_CPU
sse2_func(mat4_t, mat4_sub_s)(mat4_t m, real_t s)
{
	__m128d ms = _mm_load1_pd(&s);
	mat4_t r;
	vec4_set_result_xy(&r.x, _mm_sub_pd(vec4_load_xy(m.x), ms));
	vec4_set_result_xy(&r.y, _mm_sub_pd(vec4_load_xy(m.y), ms));
	vec4_set_result_xy(&r.z, _mm_sub_pd(vec4_load_xy(m.z), ms));
	vec4_set_result_xy(&r.w, _mm_sub_pd(vec4_load_xy(m.w), ms));
	vec4_set_result_zw(&r.x, _mm_sub_pd(vec4_load_zw(m.x), ms));
	vec4_set_result_zw(&r.y, _mm_sub_pd(vec4_load_zw(m.y), ms));
	vec4_set_result_zw(&r.z, _mm_sub_pd(vec4_load_zw(m.z), ms));
	vec4_set_result_zw(&r.w, _mm_sub_pd(vec4_load_zw(m.w), ms));
	return r;
}
#define mat4_sub_s_implemented 1
#endif

#if !mat4_sub_transpose_implemented || MATHUTIL_DETECT_CPU
sse2_func(mat4_t, mat4_sub_transpose)(mat4_t l, mat4_t r)
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

	mx_xy = vec4_load_xy(r.x);
	my_xy = vec4_load_xy(r.y);
	mz_xy = vec4_load_xy(r.z);
	mw_xy = vec4_load_xy(r.w);

	t = _mm_unpacklo_pd(mx_xy, my_xy);
	vec4_set_result_xy(&o.x, _mm_sub_pd(vec4_load_xy(l.x), t));

	t = _mm_unpacklo_pd(mz_xy, mw_xy);
	vec4_set_result_zw(&o.x, _mm_sub_pd(vec4_load_zw(l.x), t));

	t = _mm_unpackhi_pd(mx_xy, my_xy);
	vec4_set_result_xy(&o.y, _mm_sub_pd(vec4_load_xy(l.y), t));

	t = _mm_unpackhi_pd(mz_xy, mw_xy);
	vec4_set_result_zw(&o.y, _mm_sub_pd(vec4_load_zw(l.y), t));

	mx_zw = vec4_load_zw(r.x);
	my_zw = vec4_load_zw(r.y);
	mz_zw = vec4_load_zw(r.z);
	mw_zw = vec4_load_zw(r.w);

	t = _mm_unpacklo_pd(mx_zw, my_zw);
	vec4_set_result_xy(&o.z, _mm_sub_pd(vec4_load_xy(l.z), t));

	t = _mm_unpacklo_pd(mz_zw, mw_zw);
	vec4_set_result_zw(&o.z, _mm_sub_pd(vec4_load_zw(l.z), t));

	t = _mm_unpackhi_pd(mx_zw, my_zw);
	vec4_set_result_xy(&o.w, _mm_sub_pd(vec4_load_xy(l.w), t));

	t = _mm_unpackhi_pd(mz_zw, mw_zw);
	vec4_set_result_zw(&o.w, _mm_sub_pd(vec4_load_zw(l.w), t));
	return o;
}
#define mat4_sub_transpose_implemented 1
#endif

#if !mat4_mul_implemented || MATHUTIL_DETECT_CPU
sse2_func(mat4_t, mat4_mul)(mat4_t l, mat4_t r)
{
	__m128d mrx_xy = vec4_load_xy(r.x);
	__m128d mry_xy = vec4_load_xy(r.y);
	__m128d mrz_xy = vec4_load_xy(r.z);
	__m128d mrw_xy = vec4_load_xy(r.w);
	__m128d mrx_zw = vec4_load_xy(r.x);
	__m128d mry_zw = vec4_load_xy(r.y);
	__m128d mrz_zw = vec4_load_xy(r.z);
	__m128d mrw_zw = vec4_load_xy(r.w);
	__m128d mvxy;
	__m128d mvzw;

	mat4_t o;
	mvxy = vec4_load_xy(l.x);
	mvzw = vec4_load_zw(l.x);
	vec4_set_result_xy(&o.x, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), mrx_xy), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), mry_xy)),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), mrz_xy), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), mrw_xy))));
	vec4_set_result_zw(&o.x, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), mrx_zw), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), mry_zw)),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), mrz_zw), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), mrw_zw))));

	mvxy = vec4_load_xy(l.y);
	mvzw = vec4_load_zw(l.y);
	vec4_set_result_xy(&o.y, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), mrx_xy), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), mry_xy)),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), mrz_xy), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), mrw_xy))));
	vec4_set_result_zw(&o.y, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), mrx_zw), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), mry_zw)),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), mrz_zw), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), mrw_zw))));

	mvxy = vec4_load_xy(l.z);
	mvzw = vec4_load_zw(l.z);
	vec4_set_result_xy(&o.z, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), mrx_xy), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), mry_xy)),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), mrz_xy), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), mrw_xy))));
	vec4_set_result_zw(&o.z, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), mrx_zw), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), mry_zw)),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), mrz_zw), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), mrw_zw))));

	mvxy = vec4_load_xy(l.w);
	mvzw = vec4_load_zw(l.w);
	vec4_set_result_xy(&o.w, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), mrx_xy), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), mry_xy)),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), mrz_xy), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), mrw_xy))));
	vec4_set_result_zw(&o.w, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), mrx_zw), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), mry_zw)),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), mrz_zw), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), mrw_zw))));

	return o;
}
#define mat4_mul_implemented 1
#endif

#if !mat4_mul_transpose_implemented || MATHUTIL_DETECT_CPU
sse2_func(mat4_t, mat4_mul_transpose)(mat4_t l, mat4_t r)
{
	__m128d mrx_xy = vec4_load_xy(r.x);
	__m128d mry_xy = vec4_load_xy(r.y);
	__m128d mrz_xy = vec4_load_xy(r.z);
	__m128d mrw_xy = vec4_load_xy(r.w);
	__m128d mrx_zw = vec4_load_zw(r.x);
	__m128d mry_zw = vec4_load_zw(r.y);
	__m128d mrz_zw = vec4_load_zw(r.z);
	__m128d mrw_zw = vec4_load_zw(r.w);
	__m128d mvxy;
	__m128d mvzw;

	mat4_t o;
	mvxy = vec4_load_xy(l.x);
	mvzw = vec4_load_zw(l.x);
	vec4_set_result_xy(&o.x, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mrx_xy, mry_xy)), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mrx_xy, mry_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mrx_zw, mry_zw)), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mrx_zw, mry_zw)))));
	vec4_set_result_zw(&o.x, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mrz_xy, mrw_xy)), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mrz_xy, mrw_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mrz_zw, mrw_zw)), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mrz_zw, mrw_zw)))));

	mvxy = vec4_load_xy(l.y);
	mvzw = vec4_load_zw(l.y);
	vec4_set_result_xy(&o.y, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mrx_xy, mry_xy)), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mrx_xy, mry_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mrx_zw, mry_zw)), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mrx_zw, mry_zw)))));
	vec4_set_result_zw(&o.y, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mrz_xy, mrw_xy)), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mrz_xy, mrw_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mrz_zw, mrw_zw)), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mrz_zw, mrw_zw)))));

	mvxy = vec4_load_xy(l.z);
	mvzw = vec4_load_zw(l.z);
	vec4_set_result_xy(&o.z, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mrx_xy, mry_xy)), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mrx_xy, mry_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mrx_zw, mry_zw)), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mrx_zw, mry_zw)))));
	vec4_set_result_zw(&o.z, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mrz_xy, mrw_xy)), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mrz_xy, mrw_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mrz_zw, mrw_zw)), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mrz_zw, mrw_zw)))));

	mvxy = vec4_load_xy(l.w);
	mvzw = vec4_load_zw(l.w);
	vec4_set_result_xy(&o.w, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mrx_xy, mry_xy)), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mrx_xy, mry_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mrx_zw, mry_zw)), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mrx_zw, mry_zw)))));
	vec4_set_result_zw(&o.w, _mm_add_pd(
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvxy, mvxy), _mm_unpacklo_pd(mrz_xy, mrw_xy)), _mm_mul_pd(_mm_unpackhi_pd(mvxy, mvxy), _mm_unpackhi_pd(mrz_xy, mrw_xy))),
		_mm_add_pd(_mm_mul_pd(_mm_unpacklo_pd(mvzw, mvzw), _mm_unpacklo_pd(mrz_zw, mrw_zw)), _mm_mul_pd(_mm_unpackhi_pd(mvzw, mvzw), _mm_unpackhi_pd(mrz_zw, mrw_zw)))));

	return o;
}
#define mat4_mul_transpose_implemented 1
#endif

#endif // MATHUTIL_USE_DOUBLE

#endif