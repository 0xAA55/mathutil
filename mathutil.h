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

#ifndef _MATHUTIL_H_
#define _MATHUTIL_H_

#include"mathutil_conf.h"

#ifdef __cplusplus
  #include<cstdint>
  #include<cmath>
  #include<cfloat>
  namespace mathutil {
#else
  #include<stdint.h>
  #include<math.h>
  #include<float.h>
#endif

#if MATHUTIL_USE_DOUBLE
typedef double real_t, * real_p;
#define REAL_DEFINED 1

#define real_max DBL_MAX
#define real_dig DBL_DIG
#define real_mant_dig DBL_MANT_DIG
#define real_max_exp DBL_MAX_EXP
#define real_min_exp DBL_MIN_EXP
#define real_max_10_exp DBL_MAX_10_EXP
#define real_min_10_exp DBL_MIN_10_EXP

#else // !MATHUTIL_USE_DOUBLE
typedef float real_t, * real_p;
#define REAL_DEFINED 1

#define real_max FLT_MAX
#define real_dig FLT_DIG
#define real_mant_dig FLT_MANT_DIG
#define real_max_exp FLT_MAX_EXP
#define real_min_exp FLT_MIN_EXP
#define real_max_10_exp FLT_MAX_10_EXP
#define real_min_10_exp FLT_MIN_10_EXP

#endif // !MATHUTIL_USE_DOUBLE

#if !MATHUTIL_DETECT_CPU && !MATHUTIL_REFONLY

#if defined(HAVE_SSE)
#include<xmmintrin.h>
#include<immintrin.h>

#if !MATHUTIL_USE_DOUBLE
typedef union vec4_union
{
	struct
	{
		real_t x, y, z, w;
	};
	__m128 m_xyzw;
}vec4_t, * vec4_p;
#define VEC4_DEFINED 1
#define VEC4_WITH_M128_XYZW 1
#endif // !MATHUTIL_USE_DOUBLE
#endif // !defined(HAVE_SSE)

#if defined(HAVE_SSE2)
#include<emmintrin.h>
#if !defined(VEC4_DEFINED) && !defined(HAVE_AVX)
typedef union vec4_union
{
	struct
	{
		real_t x, y, z, w;
	};
	__m128d m_xy;
	__m128d m_zw;
}vec4_t, * vec4_p;
#define VEC4_DEFINED 1
#define VEC4_WITH_M128_XY_ZW 1
#endif // !defined(VEC4_DEFINED) && !defined(HAVE_AVX)
#endif // !defined(HAVE_SSE2)

#if !defined(VEC4_DEFINED) && defined(HAVE_AVX)
typedef union vec4_union
{
	struct
	{
		real_t x, y, z, w;
	};
	__m256 m_xyzw;
}vec4_t, * vec4_p;
#define VEC4_DEFINED 1
#define VEC4_WITH_M256_XYZW 1
#endif // !defined(VEC4_DEFINED) && defined(HAVE_AVX)
#endif // !MATHUTIL_DETECT_CPU && !MATHUTIL_REFONLY

#if !defined(VEC4_DEFINED)
typedef struct vec4_struct
{
	real_t x, y, z, w;
}vec4_t, * vec4_p;
#define VEC4_DEFINED 1
#endif // !defined(VEC4_DEFINED)

typedef vec4_t quat_t, * quat_p;

typedef struct mat4_struct
{
	vec4_t x, y, z, w;
}mat4_t, * mat4_p;

#define r_pi (real_t)3.141592653589793238462643383279502884
#define r_1 (real_t)1
#define r_epsilon (real_t)0.000001

#if __cplusplus
#define math_extern extern"C"
#else
#define math_extern extern
#endif

#if !defined(MATHUTIL_INTERNAL)

#define math_func(r, n) math_extern r n

// Functions for scalar numbers
math_func(real_t, r_rnd)(uint32_t* p_seed);
math_func(real_t, r_sin)(real_t x);
math_func(real_t, r_cos)(real_t x);
math_func(real_t, r_tan)(real_t x);
math_func(real_t, r_abs)(real_t x);
math_func(real_t, r_sgn)(real_t x);
math_func(real_t, r_sqr)(real_t x);
math_func(real_t, r_floor)(real_t x);
math_func(real_t, r_ceil)(real_t x);
math_func(real_t, r_atan)(real_t x);
math_func(real_t, r_exp)(real_t x);
math_func(real_t, r_log)(real_t x);
math_func(real_t, r_pow)(real_t x, real_t y);
math_func(real_t, r_mod)(real_t x, real_t y);
math_func(real_t, r_max)(real_t a, real_t b);
math_func(real_t, r_min)(real_t a, real_t b);
math_func(real_t, r_atan2)(real_t y, real_t x);
math_func(real_t, r_clamp)(real_t n, real_t min_, real_t max_);
math_func(real_t, r_lerp)(real_t a, real_t b, real_t s);
math_func(real_t, r_hermite)(real_t s);
math_func(real_t, r_slerp)(real_t a, real_t b, real_t s);

// Functions for vectors
math_func(vec4_t, vec4_t_ctor)(real_t x, real_t y, real_t z, real_t w);
math_func(vec4_t, vec4_flushcomp)(vec4_t v);
math_func(vec4_t, vec4_abs)(vec4_t v);
math_func(vec4_t, vec4_sgn)(vec4_t v);
math_func(vec4_t, vec4_invert)(vec4_t v);
math_func(real_t, vec4_length)(vec4_t v);
math_func(vec4_t, vec4_normalize)(vec4_t v);
math_func(vec4_t, vec4_scale)(vec4_t v, real_t s);
math_func(vec4_t, vec4_pow)(vec4_t v, real_t n);
math_func(real_t, vec4_dot)(vec4_t v1, vec4_t v2);
math_func(vec4_t, vec4_add)(vec4_t v1, vec4_t v2);
math_func(vec4_t, vec4_sub)(vec4_t v1, vec4_t v2);
math_func(vec4_t, vec4_mul)(vec4_t v1, vec4_t v2);
math_func(vec4_t, vec4_div)(vec4_t v1, vec4_t v2);
math_func(vec4_t, vec4_min)(vec4_t v1, vec4_t v2);
math_func(vec4_t, vec4_max)(vec4_t v1, vec4_t v2);
math_func(vec4_t, vec4_cross3)(vec4_t v1, vec4_t v2);
math_func(vec4_t, vec4_clamp)(vec4_t v, real_t min_, real_t max_);
math_func(vec4_t, vec4_rot_quat)(vec4_t v, quat_t q);
math_func(vec4_t, vec4_mul_mat4)(vec4_t v, mat4_t m);
math_func(vec4_t, vec4_mul_mat4_transpose)(vec4_t v, mat4_t m);
math_func(vec4_t, vec4_lerp)(vec4_t v1, vec4_t v2, real_t s);
math_func(vec4_t, vec4_slerp)(vec4_t v1, vec4_t v2, real_t s);

// Functions for quaternions, which can be used to describe rotation
math_func(quat_t, quat_t_ctor)(real_t x, real_t y, real_t z, real_t w);
math_func(quat_t, quat_flushcomp)(quat_t q);
math_func(quat_t, quat_rot_axis)(vec4_t axis, real_t angle);
math_func(quat_t, quat_mul)(quat_t q1, quat_t q2);
math_func(quat_t, quat_add_vec)(quat_t q, vec4_t v, real_t s);

// Functions for matrices
math_func(mat4_t, mat4_t_ctor)(vec4_t mx, vec4_t my, vec4_t mz, vec4_t mw);
math_func(mat4_t, mat4_flushcomp)(mat4_t m);
math_func(mat4_t, mat4_rot_x)(real_t angle);
math_func(mat4_t, mat4_rot_y)(real_t angle);
math_func(mat4_t, mat4_rot_z)(real_t angle);
math_func(mat4_t, mat4_rot_axis)(vec4_t axis, real_t angle);
math_func(mat4_t, mat4_rot_euler)(real_t yaw, real_t pitch, real_t roll);
math_func(mat4_t, mat4_from_quat)(quat_t q);
math_func(mat4_t, mat4_from_quat_transpose)(quat_t q);
math_func(mat4_t, mat4_transpose)(mat4_t m);
math_func(real_t, mat4_det)(mat4_t m);
math_func(int, mat4_inverse)(mat4_t m, real_p det_out, mat4_p mat_out);
math_func(mat4_t, mat4_add)(mat4_t l, mat4_t r);
math_func(mat4_t, mat4_add_s)(mat4_t m, real_t s);
math_func(mat4_t, mat4_add_transpose)(mat4_t l, mat4_t r);
math_func(mat4_t, mat4_sub)(mat4_t l, mat4_t r);
math_func(mat4_t, mat4_sub_s)(mat4_t m, real_t s);
math_func(mat4_t, mat4_sub_transpose)(mat4_t l, mat4_t r);
math_func(mat4_t, mat4_mul)(mat4_t l, mat4_t r);
math_func(mat4_t, mat4_mul_s)(mat4_t m, real_t s);
math_func(mat4_t, mat4_mul_transpose)(mat4_t l, mat4_t r);

#ifdef __cplusplus
	typedef real_t real;
	class vec4;
	class quat;
	class mat4;
	typedef vec4_t& vec4_r;
	typedef quat_t& quat_r;
	typedef mat4_t& mat4_r;
	class vec4
	{
	private:
		vec4_t m_cv;
		friend class quat;
	public:
		inline vec4();
		inline vec4(const vec4& v);
		inline vec4(vec4_t& cv);
		inline vec4(real r);
		inline vec4(real x, real y, real z, real w);
		inline vec4(real r[4]);
		inline vec4(const quat& q);

		inline vec4_r const cv();
		inline operator vec4_t ();

		inline void get_components(real& x, real& y, real& z, real& w);
		inline vec4& set_components(real x, real y, real z, real w);

		inline vec4 const abs();
		inline vec4 const sgn();
		inline vec4 const invert();
		inline real const length();
		inline vec4 const normalize();
		inline vec4 const scale(real s);
		inline vec4 const pow(real n);
		inline real const dot(vec4& v);
		inline vec4 const add(vec4& v);
		inline vec4 const sub(vec4& v);
		inline vec4 const mul(real s);
		inline vec4 const mul(vec4& v);
		inline vec4 const mul(quat& q);
		inline vec4 const mul(mat4& m);
		inline vec4 const mul_transpose(mat4& m);
		inline vec4 const div(vec4& v);
		inline vec4 const min(vec4& v);
		inline vec4 const max(vec4& v);
		inline vec4 const cross(vec4&v);
		inline vec4 const clamp(real min_, real max_);
		inline vec4 const lerp(vec4& to, real s);
		inline vec4 const slerp(vec4& to, real s);

		inline vec4 const operator +(vec4& v) {  return add(v); }
		inline vec4 const operator -(vec4& v) { return sub(v); }
		inline vec4 const operator -() { return invert(); }
		inline vec4 const operator *(vec4& v) { return mul(v); }
		inline vec4 const operator *(real s) { return mul(s); }
		inline vec4 const operator *(quat& q) { return mul(q); }
		inline vec4 const operator *(mat4& m) { return mul(m); }
		inline vec4 const operator /(vec4& v) { return div(v); }
		inline vec4 const operator %(vec4& v) { return cross(v); }
		inline vec4& operator += (vec4& v);
		inline vec4& operator -= (vec4& v);
		inline vec4& operator *= (real s);
		inline vec4& operator *= (vec4& v);
		inline vec4& operator *= (quat& q);
		inline vec4& operator *= (mat4& m);
		inline vec4& operator /= (vec4& v);
		inline vec4& operator ++ ();
		inline vec4& operator -- ();
		inline vec4 operator ++ (int) { vec4 v(*this); ++(*this); return v; }
		inline vec4 operator -- (int) { vec4 v(*this); --(*this); return v; }
	};

	inline vec4 abs(vec4& v) { return v.abs(); }
	inline vec4 sgn(vec4& v) { return v.sgn(); }
	inline vec4 invert(vec4& v) { return v.invert(); }
	inline real length(vec4& v) { return v.length(); }
	inline vec4 normalize(vec4& v) { return v.normalize(); }
	inline vec4 scale(vec4& v, real s) { return v.scale(s); }
	inline vec4 pow(vec4& v, real n) { return v.pow(n); }
	inline real dot(vec4& v1, vec4& v2) { return v1.dot(v2); }
	inline vec4 add(vec4& v1, vec4& v2) { return v1.add(v2); }
	inline vec4 sub(vec4& v1, vec4& v2) { return v1.sub(v2); }
	inline vec4 mul(vec4& v, real s) { return v.mul(s); }
	inline vec4 mul(vec4& v1, vec4& v2) { return v1.mul(v2); }
	inline vec4 mul(vec4& v, quat& q) { return v.mul(q); }
	inline vec4 mul(vec4& v, mat4& m) { return v.mul(m); }
	inline vec4 mul_transpose(vec4& v, mat4& m) { return v.mul_transpose(m); }
	inline vec4 div(vec4& v1, vec4& v2) { return v1.div(v2); }
	inline vec4 min(vec4& v1, vec4& v2) { return v1.min(v2); }
	inline vec4 max(vec4& v1, vec4& v2) { return v1.max(v2); }
	inline vec4 cross(vec4& v1, vec4& v2) { return v1.cross(v2); }
	inline vec4 clamp(vec4& v, real min_, real max_) { return v.clamp(min_, max_); }
	inline vec4 lerp(vec4& v1, vec4& v2, real s) { return v1.lerp(v2, s); }
	inline vec4 slerp(vec4& v1, vec4& v2, real s) { return v1.slerp(v2, s); }
	
	class quat
	{
	private:
		quat_t m_cq;
		friend class vec4;
	public:
		inline quat();
		inline quat(const quat& q);
		inline quat(quat_t &cq);
		inline quat(real x, real y, real z, real w);
		inline quat(real r[4]);
		inline quat(const vec4 &v);
		inline quat(vec4 &axis, real angle);

		inline quat_r const cq();
		inline operator quat_t ();
		inline operator mat4 ();

		inline void get_components(real& x, real& y, real& z, real& w);
		inline quat& set_components(real x, real y, real z, real w);

		inline static quat rot_axis(vec4 &axis, real angle);
		inline quat const mul(quat &q);
		inline quat const normalize();
		inline quat const add_vec(vec4& v, real s);

		inline quat const operator -() {  return quat(vec4(*this).invert()); }
		inline quat const operator +(vec4& v) { return add_vec(v, 1); }
		inline quat const operator -(vec4 &v) { return add_vec(v, -1); }
		inline quat const operator *(quat& q) { return mul(q); }
		inline quat& operator += (vec4 &v);
		inline quat& operator -= (vec4 &v);
		inline quat& operator *= (quat &q);
	};

	inline quat rot_axis_q(vec4 &axis, real angle) { return quat::rot_axis(axis, angle); }
	inline quat mul(quat& q1, quat& q2) { return q1.mul(q2); }
	inline quat normalize(quat& q) { return q.normalize(); }
	inline quat add_vec(quat& q, vec4& v, real s) { return q.add_vec(v, s); }

	class mat4
	{
	private:
		mat4_t m_cm;
	public:
		inline mat4();
		inline mat4(real r);
		inline mat4(const mat4& m);
		inline mat4(mat4_t& cm);
		inline mat4
		(
			real xx, real xy, real xz, real xw,
			real yx, real yy, real yz, real yw,
			real zx, real zy, real zz, real zw,
			real wx, real wy, real wz, real ww
		);
		inline mat4(real r[4][4]);
		inline mat4(vec4& x, vec4& y, vec4& z, vec4& w);

		inline mat4_r const cm();
		inline operator mat4_t ();

		inline void get_components (vec4 &x, vec4 &y, vec4 &z, vec4 &w);
		inline mat4& set_components(vec4 x, vec4 y, vec4 z, vec4 w);

		inline static mat4 rot_x(real angle);
		inline static mat4 rot_y(real angle);
		inline static mat4 rot_z(real angle);
		inline static mat4 rot_axis(vec4 &axis, real angle);
		inline static mat4 rot_euler(real yaw, real pitch, real roll);
		inline static mat4 from_quat(quat &q);
		inline static mat4 from_quat_transpose(quat& q);
		inline mat4 const transpose();
		inline real const det();
		inline mat4 const inverse();
		inline bool const inverse(mat4& out);
		inline bool const inverse(mat4& out, real& det_out);
		inline mat4 const add(real s);
		inline mat4 const add(mat4& m);
		inline mat4 const add_transpose(mat4& m);
		inline mat4 const sub(real s);
		inline mat4 const sub(mat4& m);
		inline mat4 const sub_transpose(mat4& m);
		inline mat4 const mul(real s);
		inline mat4 const mul(mat4& m);
		inline mat4 const mul_transpose(mat4& m);

		inline mat4 const operator -()
		{
			vec4 x, y, z, w;
			get_components(x, y, z, w);
			x.invert();
			y.invert();
			z.invert();
			w.invert();
			return mat4(x, y, z, w);
		}
		inline mat4 const operator +(mat4& m) { return add(m); }
		inline mat4 const operator -(mat4& m) { return sub(m); }
		inline mat4 const operator *(mat4& m) { return mul(m); }
		inline mat4 const operator +(real s) { return add(s); }
		inline mat4 const operator -(real s) { return sub(s); }
		inline mat4 const operator *(real s) { return mul(s); }

		inline mat4& operator +=(mat4& m);
		inline mat4& operator -=(mat4& m);
		inline mat4& operator *=(mat4& m);
		inline mat4& operator +=(real s);
		inline mat4& operator -=(real s);
		inline mat4& operator *=(real s);
		inline mat4& operator ++();
		inline mat4& operator --();
		inline mat4 operator ++(int) { mat4 m(*this); ++(*this); return m; }
		inline mat4 operator --(int) { mat4 m(*this); --(*this); return m; }
	};

	inline mat4 rot_x(real angle) { return mat4::rot_x(angle); }
	inline mat4 rot_y(real angle) { return mat4::rot_y(angle); }
	inline mat4 rot_z(real angle) { return mat4::rot_z(angle); }
	inline mat4 rot_axis_m(vec4& axis, real angle) { return mat4::rot_axis(axis, angle); }
	inline mat4 rot_euler(real yaw, real pitch, real roll) { return mat4::rot_euler(yaw, pitch, roll); }
	inline mat4 from_quat(quat& q) { return mat4::from_quat(q); }
	inline mat4 from_quat_transpose(quat& q) { return mat4::from_quat_transpose(q); }
	inline mat4 transpose(mat4& m) { return m.transpose(); }
	inline real det(mat4& m) { return m.det(); }
	inline mat4 inverse(mat4& m) { return m.inverse(); }
	inline bool inverse(mat4& m, mat4& out) { return m.inverse(out); }
	inline bool inverse(mat4& m, mat4& out, real& det_out) { return m.inverse(out, det_out); }
	inline mat4 add(mat4& m, real s) { return m.add(s); }
	inline mat4 add(mat4& l, mat4& r) { return l.add(r); }
	inline mat4 add_transpose(mat4& l, mat4& r) { return l.add_transpose(r); }
	inline mat4 sub(mat4& m, real s) { return m.sub(s); }
	inline mat4 sub(mat4& l, mat4& r) { return l.sub(r); }
	inline mat4 sub_transpose(mat4& l, mat4& r) { return l.sub_transpose(r); }
	inline mat4 mul(mat4& m, real s) { return m.mul(s); }
	inline mat4 mul(mat4& l, mat4& r) { return l.mul(r); }
	inline mat4 mul_transpose(mat4& l, mat4& r) { return l.mul_transpose(r); }

	inline vec4::vec4() {}
	inline vec4::vec4(real r) : m_cv(vec4_t_ctor(r, r, r, r)) {}
	inline vec4::vec4(const vec4& v) : m_cv(v.m_cv) {}
	inline vec4::vec4(vec4_t& cv) : m_cv(cv) {}
	inline vec4::vec4(real x, real y, real z, real w) : m_cv(vec4_t_ctor(x, y, z, w)) {}
	inline vec4::vec4(real r[4]) : m_cv(vec4_t_ctor(r[0], r[1], r[2], r[3])) {}
	inline vec4::vec4(const quat& q) : m_cv(q.m_cq) {}

	inline vec4_r const vec4::cv() { return m_cv; }
	inline vec4::operator vec4_t () { return m_cv; }

	inline void vec4::get_components(real & x, real& y, real& z, real& w)
	{
		m_cv = vec4_flushcomp(m_cv);
		x = m_cv.x;
		y = m_cv.y;
		z = m_cv.z;
		w = m_cv.w;
	}
	inline vec4& vec4::set_components(real x, real y, real z, real w)
	{
		m_cv = vec4_t_ctor(x, y, z, w);
		return *this;
	}

	inline vec4 const vec4::abs() { return vec4_abs(m_cv); }
	inline vec4 const vec4::sgn() { return vec4_sgn(m_cv); }
	inline vec4 const vec4::invert() { return vec4_invert(m_cv); }
	inline real const vec4::length() { return vec4_length(m_cv); }
	inline vec4 const vec4::normalize() { return vec4_normalize(m_cv); }
	inline vec4 const vec4::scale(real s) { return vec4_scale(m_cv, s); }
	inline vec4 const vec4::pow(real n) { return vec4_pow(m_cv, n); }
	inline real const vec4::dot(vec4& v) { return vec4_dot(m_cv, v); }
	inline vec4 const vec4::add(vec4& v) { return vec4_add(m_cv, v); }
	inline vec4 const vec4::sub(vec4& v) { return vec4_sub(m_cv, v); }
	inline vec4 const vec4::mul(real s) { return vec4_scale(m_cv, s); }
	inline vec4 const vec4::mul(vec4& v) { return vec4_mul(m_cv, v); }
	inline vec4 const vec4::mul(quat& q) { return vec4_rot_quat(m_cv, q); }
	inline vec4 const vec4::mul(mat4& m) { return vec4_mul_mat4(m_cv, m); }
	inline vec4 const vec4::mul_transpose(mat4& m) { return vec4_mul_mat4_transpose(m_cv, m); }
	inline vec4 const vec4::div(vec4& v) { return vec4_div(m_cv, v); }
	inline vec4 const vec4::min(vec4& v) { return vec4_min(m_cv, v); }
	inline vec4 const vec4::max(vec4& v) { return vec4_max(m_cv, v); }
	inline vec4 const vec4::cross(vec4& v) { return vec4_cross3(m_cv, v); }
	inline vec4 const vec4::clamp(real min_, real max_) { return vec4_clamp(m_cv, min_, max_); }
	inline vec4 const vec4::lerp(vec4& to, real s) { return vec4_lerp(m_cv, to, s); }
	inline vec4 const vec4::slerp(vec4& to, real s) { return vec4_slerp(m_cv, to, s); }

	inline vec4& vec4::operator += (vec4& v) { m_cv = vec4_add(m_cv, v); return *this; }
	inline vec4& vec4::operator -= (vec4& v) { m_cv = vec4_sub(m_cv, v); return *this; }
	inline vec4& vec4::operator *= (real s) { m_cv = vec4_scale(m_cv, s); return *this; }
	inline vec4& vec4::operator *= (vec4& v) { m_cv = vec4_mul(m_cv, v); return *this; }
	inline vec4& vec4::operator *= (quat& q) { m_cv = vec4_rot_quat(m_cv, q); return *this; }
	inline vec4& vec4::operator *= (mat4& m) { m_cv = vec4_mul_mat4(m_cv, m); return *this; }
	inline vec4& vec4::operator /= (vec4& v) { m_cv = vec4_div(m_cv, v); return *this; }
	inline vec4& vec4::operator ++ () { m_cv = vec4_add(m_cv, vec4_t_ctor(1, 1, 1, 1)); return *this; }
	inline vec4& vec4::operator -- () { m_cv = vec4_sub(m_cv, vec4_t_ctor(1, 1, 1, 1)); return *this; }

	inline quat::quat() {}
	inline quat::quat(const quat& q) : m_cq(q.m_cq) {}
	inline quat::quat(quat_t& cq) : m_cq(cq) {}
	inline quat::quat(real x, real y, real z, real w) : m_cq(quat_t_ctor(x, y, z, w)) {}
	inline quat::quat(real r[4]) : m_cq(quat_t_ctor(r[0], r[1], r[2], r[3])) {}
	inline quat::quat(const vec4& v) : m_cq(v.m_cv) {}
	inline quat::quat(vec4& axis, real angle) : m_cq(quat_rot_axis(axis, angle)) {}

	inline quat_r const quat::cq() { return m_cq; }
	inline quat::operator quat_t () { return m_cq; }
	inline quat::operator mat4 () { return mat4_from_quat(m_cq); }

	inline void quat::get_components(real& x, real& y, real& z, real& w)
	{
		m_cq = quat_flushcomp(m_cq);
		x = m_cq.x;
		y = m_cq.y;
		z = m_cq.z;
		w = m_cq.w;
	}
	inline quat& quat::set_components(real x, real y, real z, real w)
	{
		m_cq = quat_t_ctor(x, y, z, w);
		return *this;
	}

	inline quat quat::rot_axis(vec4& axis, real angle) { return quat_rot_axis(axis, angle); }
	inline quat const quat::mul(quat& q) { return quat_mul(m_cq, q); }
	inline quat const quat::normalize() { return vec4_normalize(m_cq); }
	inline quat const quat::add_vec(vec4& v, real s) { return quat_add_vec(m_cq, v, s); }

	inline quat& quat::operator += (vec4& v) { m_cq = quat_add_vec(m_cq, v, 1); return *this; }
	inline quat& quat::operator -= (vec4& v) { m_cq = quat_add_vec(m_cq, v, -1); return *this; }
	inline quat& quat::operator *= (quat& q) { m_cq = quat_mul(m_cq, q); return *this; }

	inline mat4::mat4() {}
	inline mat4::mat4(real r) : m_cm(mat4_t_ctor(
		vec4_t_ctor(r, r, r, r),
		vec4_t_ctor(r, r, r, r),
		vec4_t_ctor(r, r, r, r),
		vec4_t_ctor(r, r, r, r)
	)) {}
	inline mat4::mat4(const mat4& m) : m_cm(m.m_cm) {}
	inline mat4::mat4(mat4_t& cm) : m_cm(cm) {}
	inline mat4::mat4
	(
		real xx, real xy, real xz, real xw,
		real yx, real yy, real yz, real yw,
		real zx, real zy, real zz, real zw,
		real wx, real wy, real wz, real ww
	) : m_cm(mat4_t_ctor(
			vec4_t_ctor(xx, xy, xz, xw),
			vec4_t_ctor(yx, yy, yz, yw),
			vec4_t_ctor(zx, zy, zz, zw),
			vec4_t_ctor(wx, wy, wz, ww)
	)) {}
	inline mat4::mat4(real r[4][4]) : m_cm(mat4_t_ctor(
		vec4_t_ctor(r[0][0], r[0][1], r[0][2], r[0][3]),
		vec4_t_ctor(r[1][0], r[1][1], r[1][2], r[1][3]),
		vec4_t_ctor(r[2][0], r[2][1], r[2][2], r[2][3]),
		vec4_t_ctor(r[3][0], r[3][1], r[3][2], r[3][3])
	)) {}
	inline mat4::mat4(vec4& x, vec4& y, vec4& z, vec4& w) : m_cm(mat4_t_ctor(x, y, z, w)){}

	inline mat4_r const mat4::cm() { return m_cm; }
	inline mat4::operator mat4_t () { return m_cm; }

	inline void mat4::get_components(vec4& x, vec4& y, vec4& z, vec4& w)
	{
		m_cm = mat4_flushcomp(m_cm);
		x = m_cm.x;
		y = m_cm.y;
		z = m_cm.z;
		w = m_cm.w;
	}
	inline mat4& mat4::set_components(vec4 x, vec4 y, vec4 z, vec4 w)
	{
		m_cm = mat4_t_ctor(x, y, z, w);
		return *this;
	}

	inline mat4 mat4::rot_x(real angle) { return mat4_rot_x(angle); }
	inline mat4 mat4::rot_y(real angle) { return mat4_rot_y(angle); }
	inline mat4 mat4::rot_z(real angle) { return mat4_rot_z(angle); }
	inline mat4 mat4::rot_axis(vec4& axis, real angle) { return mat4_rot_axis(axis, angle); }
	inline mat4 mat4::rot_euler(real yaw, real pitch, real roll) { return mat4_rot_euler(yaw, pitch, roll); }
	inline mat4 mat4::from_quat(quat& q) { return mat4_from_quat(q); }
	inline mat4 mat4::from_quat_transpose(quat& q) { return mat4_from_quat_transpose(q); }
	inline mat4 const mat4::transpose() { return mat4_transpose(m_cm); }
	inline real const mat4::det() { return mat4_det(m_cm); }
	inline mat4 const mat4::inverse() { mat4_t r; if (mat4_inverse(m_cm, NULL, &r)) return r; else return mat4((real)0); };
	inline bool const mat4::inverse(mat4& out) { return static_cast<bool>(mat4_inverse(m_cm, NULL, &out.m_cm)); }
	inline bool const mat4::inverse(mat4& out, real& det_out) { return static_cast<bool>(mat4_inverse(m_cm, &det_out, &out.m_cm)); }
	inline mat4 const mat4::add(real s) { return mat4_add_s(m_cm, s); }
	inline mat4 const mat4::add(mat4& m) { return mat4_add(m_cm, m); }
	inline mat4 const mat4::add_transpose(mat4& m) { return mat4_add_transpose(m_cm, m); }
	inline mat4 const mat4::sub(real s) { return mat4_sub_s(m_cm, s); }
	inline mat4 const mat4::sub(mat4& m) { return mat4_sub(m_cm, m); }
	inline mat4 const mat4::sub_transpose(mat4& m) { return mat4_sub_transpose(m_cm, m); }
	inline mat4 const mat4::mul(real s) { return mat4_mul_s(m_cm, s); }
	inline mat4 const mat4::mul(mat4& m) { return mat4_mul(m_cm, m); }
	inline mat4 const mat4::mul_transpose(mat4& m) { return mat4_mul_transpose(m_cm, m); }

	inline mat4& mat4::operator +=(mat4& m) { m_cm = mat4_add(m_cm, m); return *this; }
	inline mat4& mat4::operator -=(mat4& m) { m_cm = mat4_sub(m_cm, m); return *this; }
	inline mat4& mat4::operator *=(mat4& m) { m_cm = mat4_mul(m_cm, m); return *this; }
	inline mat4& mat4::operator +=(real s) { m_cm = mat4_add_s(m_cm, s); return *this; }
	inline mat4& mat4::operator -=(real s) { m_cm = mat4_sub_s(m_cm, s); return *this; }
	inline mat4& mat4::operator *=(real s) { m_cm = mat4_mul_s(m_cm, s); return *this; }
	inline mat4& mat4::operator ++() { m_cm = mat4_add_s(m_cm, 1); return *this; }
	inline mat4& mat4::operator --() { m_cm = mat4_sub_s(m_cm, 1); return *this; }
}
#else
#  define vec4 vec4_t_ctor
#  define quat quat_t_ctor
#  define mat4 mat4_t_ctor
#endif

#endif

#endif // !_MATHUTIL_H_
