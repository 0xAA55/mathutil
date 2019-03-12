#include"mathutil_ref.h"
#include<stdlib.h>

#if MATHUTIL_DETECT_SIMD
#define base_func(r,n) r n ## _ref
#else
#define base_func(r,n) r n
#endif

base_func(real_t,r_rnd)(uint32_t *p_seed)
{
	return (real_t)(((*p_seed = *p_seed * 214013L + 2531011L) >> 16) & 0x7fff) / (real_t)RAND_MAX;
}

base_func(real_t,r_sin)(real_t x)
{
	return (real_t)sin(x);
}

base_func(real_t,r_cos)(real_t x)
{
	return (real_t)cos(x);
}

base_func(real_t,r_abs)(real_t x)
{
	return (real_t)fabs(x);
}

base_func(real_t,r_sgn)(real_t x)
{
	return x >= 0 ? r_1 : -r_1;
}

base_func(real_t,r_sqr)(real_t x)
{
	return (real_t)sqrt(x);
}

base_func(real_t,r_floor)(real_t x)
{
	return (real_t)floor(x);
}

base_func(real_t,r_ceil)(real_t x)
{
	return (real_t)ceil(x);
}

base_func(real_t,r_atan)(real_t x)
{
	return (real_t)atan(x);
}

base_func(real_t,r_exp)(real_t x)
{
	return (real_t)exp(x);
}

base_func(real_t,r_log)(real_t x)
{
	return (real_t)log(x);
}

base_func(real_t,r_pow)(real_t x, real_t y)
{
	return (real_t)pow(x, y);
}

base_func(real_t,r_max)(real_t a, real_t b)
{
	return a > b ? a : b;
}

base_func(real_t,r_min)(real_t a, real_t b)
{
	return a < b ? a : b;
}

base_func(real_t,r_mod)(real_t x, real_t y)
{
	return (real_t)fmod(x, y);
}

base_func(real_t,r_atan2)(real_t x, real_t y)
{
	return (real_t)atan2(x, y);
}

base_func(real_t,r_clamp)(real_t n, real_t min, real_t max)
{
	if(n < min) n = min;
	if(n > max) n = max;
	return n;
}

base_func(real_t,r_lerp)(real_t a, real_t b, real_t s)
{
	return a + (b - a) * s;
}

base_func(real_t,r_hermite)(real_t s)
{
	return s * s * (3 - 2 * s);
}

base_func(real_t,r_slerp)(real_t a, real_t b, real_t s)
{
	return r_lerp(a, b, r_hermite(r_clamp(s, 0, 1)));
}

base_func(vec4_t,vec4)(real_t x, real_t y, real_t z, real_t w)
{
	vec4_t v = {x, y, z, w};
	return v;
}

base_func(vec4_t,vec4_abs)(vec4_t v)
{
	return vec4
	(
		r_abs(v.x),
		r_abs(v.y),
		r_abs(v.z),
		r_abs(v.w)
	);
}

base_func(vec4_t,vec4_sgn)(vec4_t v)
{
	return vec4
	(
		r_sgn(v.x),
		r_sgn(v.y),
		r_sgn(v.z),
		r_sgn(v.w)
	);
}

base_func(vec4_t,vec4_invert)(vec4_t v)
{
	return vec4
	(
		-(v.x),
		-(v.y),
		-(v.z),
		-(v.w)
	);
}

base_func(real_t,vec4_length)(vec4_t v)
{
	return r_sqr(vec4_dot(v, v));
}

base_func(vec4_t,vec4_normalize)(vec4_t v)
{
	real_t l = vec4_length(v);
	if(l < r_epsilon) return vec4(0,0,0,0);
	else return vec4_scale(v, r_1 / l);
}

base_func(vec4_t,vec4_scale)(vec4_t v, real_t s)
{
	return vec4
	(
		s * (v.x),
		s * (v.y),
		s * (v.z),
		s * (v.w)
	);
}

base_func(vec4_t,vec4_pow)(vec4_t v, real_t n)
{
	return vec4
	(
		r_pow(v.x, n),
		r_pow(v.y, n),
		r_pow(v.z, n),
		r_pow(v.w, n)
	);
}

base_func(real_t,vec4_dot)(vec4_t v1, vec4_t v2)
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z + v1.w * v2.w;
}

base_func(vec4_t,vec4_add)(vec4_t v1, vec4_t v2)
{
	return vec4
	(
		v1.x + v2.x,
		v1.y + v2.y,
		v1.z + v2.z,
		v1.w + v2.w
	);
}

base_func(vec4_t,vec4_sub)(vec4_t v1, vec4_t v2)
{
	return vec4
	(
		v1.x - v2.x,
		v1.y - v2.y,
		v1.z - v2.z,
		v1.w - v2.w
	);
}

base_func(vec4_t,vec4_mul)(vec4_t v1, vec4_t v2)
{
	return vec4
	(
		v1.x * v2.x,
		v1.y * v2.y,
		v1.z * v2.z,
		v1.w * v2.w
	);
}

base_func(vec4_t,vec4_div)(vec4_t v1, vec4_t v2)
{
	return vec4
	(
		v1.x / v2.x,
		v1.y / v2.y,
		v1.z / v2.z,
		v1.w / v2.w
	);
}

base_func(vec4_t,vec4_min)(vec4_t v1, vec4_t v2)
{
	return vec4
	(
		r_min(v1.x, v2.x),
		r_min(v1.y, v2.y),
		r_min(v1.z, v2.z),
		r_min(v1.w, v2.w)
	);
}

base_func(vec4_t,vec4_max)(vec4_t v1, vec4_t v2)
{
	return vec4
	(
		r_max(v1.x, v2.x),
		r_max(v1.y, v2.y),
		r_max(v1.z, v2.z),
		r_max(v1.w, v2.w)
	);
}

base_func(vec4_t,vec4_cross)(vec4_t v1, vec4_t v2)
{
	return vec4
	(
		v1.y * v2.z - v1.z * v2.y,
		v1.z * v2.x - v1.x * v2.z,
		v1.x * v2.y - v1.y * v2.x,
		0
	);
}

base_func(vec4_t,vec4_mul_mat4)(vec4_t v, mat4_t m)
{
	return vec4
	(
		v.x * m.x.x + v.y * m.y.x + v.z * m.z.x + v.w * m.w.x,
		v.x * m.x.y + v.y * m.y.y + v.z * m.z.y + v.w * m.w.y,
		v.x * m.x.z + v.y * m.y.z + v.z * m.z.z + v.w * m.w.z,
		v.x * m.x.w + v.y * m.y.w + v.z * m.z.w + v.w * m.w.w
	);
}

base_func(vec4_t,vec4_mul_mat4_transpose)(vec4_t v, mat4_t m)
{
	return vec4
	(
		v.x * m.x.x + v.y * m.x.y + v.z * m.x.z + v.w * m.x.w,
		v.x * m.y.x + v.y * m.y.y + v.z * m.y.z + v.w * m.y.w,
		v.x * m.z.x + v.y * m.z.y + v.z * m.z.z + v.w * m.z.w,
		v.x * m.w.x + v.y * m.w.y + v.z * m.w.z + v.w * m.w.w
	);
}

base_func(vec4_t,vec4_lerp)(vec4_t v1, vec4_t v2, real_t s)
{
	return vec4
	(
		r_lerp(v1.x, v2.x, s),
		r_lerp(v1.y, v2.y, s),
		r_lerp(v1.z, v2.z, s),
		r_lerp(v1.w, v2.w, s)
	);
}

base_func(vec4_t,vec4_clamp)(vec4_t v, real_t min_, real_t max_)
{
	return vec4
	(
		r_clamp(v.x, min_, max_),
		r_clamp(v.y, min_, max_),
		r_clamp(v.z, min_, max_),
		r_clamp(v.w, min_, max_)
	);
}

base_func(vec4_t,vec4_slerp)(vec4_t v1, vec4_t v2, real_t s)
{
	s = r_hermite(r_clamp(s, 0, 1));
	return vec4
	(
		r_lerp(v1.x, v2.x, s),
		r_lerp(v1.y, v2.y, s),
		r_lerp(v1.z, v2.z, s),
		r_lerp(v1.w, v2.w, s)
	);
}

base_func(vec4_t,vec4_rot_quat)(vec4_t v, quat_t q)
{
	return quat_mul(quat_mul(q, vec4(v.x, v.y, v.z, 0)), vec4(-q.x, -q.y, -q.z, q.w));
}

base_func(quat_t,quat)(real_t x, real_t y, real_t z, real_t w)
{
	quat_t r = {x, y, z, w};
	return r;
}

base_func(quat_t,quat_rot_axis)(vec4_t axis, real_t angle)
{
	real_t ha = angle / 2;
	real_t sin_ha = r_sin(ha);
	return quat
	(
		axis.x * sin_ha,
		axis.y * sin_ha,
		axis.z * sin_ha,
		r_cos(ha)
	);
}

base_func(quat_t,quat_mul)(quat_t q1, quat_t q2)
{
	return quat
	(
		q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y,
		q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x,
		q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w,
		q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z
	);
}

base_func(quat_t,quat_add_vec)(quat_t q, vec4_t v, real_t s)
{
	return vec4_add(vec4_scale(quat_mul(vec4(v.x * s, v.y * s, v.z * s, 0), q), (real_t)0.5), q);
}

base_func(mat4_t,mat4)(vec4_t x, vec4_t y, vec4_t z, vec4_t w)
{
	mat4_t r;
	r.x = x;
	r.y = y;
	r.z = z;
	r.w = w;
	return r;
}

base_func(mat4_t,mat4_rot_x)(real_t angle)
{
	real_t sa = r_sin(angle);
	real_t ca = r_cos(angle);
	return mat4
	(
		vec4(1,  0,  0, 0),
		vec4(0, ca, sa, 0),
		vec4(0,-sa, ca, 0),
		vec4(0,  0,  0, 1)
	);
}

base_func(mat4_t,mat4_rot_y)(real_t angle)
{
	real_t sa = r_sin(angle);
	real_t ca = r_cos(angle);
	return mat4
	(
		vec4(ca, 0,-sa, 0),
		vec4(0,  1,  0, 0),
		vec4(sa, 0, ca, 0),
		vec4(0,  0,  0, 1)
	);
}

base_func(mat4_t,mat4_rot_z)(real_t angle)
{
	real_t sa = r_sin(angle);
	real_t ca = r_cos(angle);
	return mat4
	(
		vec4( ca, sa, 0, 0),
		vec4(-sa, ca, 0, 0),
		vec4(  0,  0, 1, 0),
		vec4(  0,  0, 0, 1)
	);
}

base_func(mat4_t,mat4_rot_axis)(vec4_t axis, real_t angle)
{
	real_t sa = r_sin(angle);
	real_t ca = r_cos(angle);
	vec4_t v = vec4_normalize(axis);
	return mat4
	(
		vec4((1 - ca) * v.x * v.x + ca, (1 - ca) * v.y * v.x + sa * v.z, (1 - ca) * v.z * v.x - sa * v.y, 0),
		vec4((1 - ca) * v.x * v.y - sa * v.z, (1 - ca) * v.y * v.y + ca, (1 - ca) * v.z * v.y + sa * v.x, 0),
		vec4((1 - ca) * v.x * v.z + sa * v.y, (1 - ca) * v.y * v.z - sa * v.x, (1 - ca) * v.z * v.z + ca, 0),
		vec4(0, 0, 0, 1)
	);
}

base_func(mat4_t,mat4_rot_euler)(real_t yaw, real_t pitch, real_t roll)
{
	real_t
		sr = r_sin(yaw),
		cr = r_cos(yaw),
		sp = r_sin(pitch),
		cp = r_cos(pitch),
		sy = r_sin(roll),
		cy = r_cos(roll);

	real_t
		srcp = sr * cp,
		srsp = sr * sp,
		crcp = cr * cp,
		crsp = cr * sp;

	return mat4
	(
		vec4( cr * cy + srsp * sy, srcp,-sy * cr + srsp * cy, 0),
		vec4(-sr * cy + crsp * sy, crcp, sy * sr + crsp * cy, 0),
		vec4(sy * cp, -sp, cp * cy, 0),
		vec4(0, 0, 0, 1)
	);
}

base_func(mat4_t,mat4_from_quat)(quat_t q)
{
	return mat4
	(
		vec4(1 - 2 * q.y * q.y - 2 * q.z * q.z, 0 + 2 * q.x * q.y + 2 * q.z * q.w, 0 + 2 * q.x * q.z - 2 * q.y * q.w, 0),
		vec4(0 + 2 * q.x * q.y - 2 * q.z * q.w, 1 - 2 * q.x * q.x - 2 * q.z * q.z, 0 + 2 * q.y * q.z + 2 * q.x * q.w, 0),
		vec4(0 + 2 * q.x * q.z + 2 * q.y * q.w, 0 + 2 * q.y * q.z - 2 * q.x * q.w, 1 - 2 * q.x * q.x - 2 * q.y * q.y, 0),
		vec4(0, 0, 0, 1)
	);
}

base_func(mat4_t,mat4_from_quat_transpose)(quat_t q)
{
	return mat4
	(
		vec4(1 - 2 * q.y * q.y - 2 * q.z * q.z, 0 + 2 * q.x * q.y - 2 * q.z * q.w, 0 + 2 * q.x * q.z + 2 * q.y * q.w, 0),
		vec4(0 + 2 * q.x * q.y + 2 * q.z * q.w, 1 - 2 * q.x * q.x - 2 * q.z * q.z, 0 + 2 * q.y * q.z - 2 * q.x * q.w, 0),
		vec4(0 + 2 * q.x * q.z - 2 * q.y * q.w, 0 + 2 * q.y * q.z + 2 * q.x * q.w, 1 - 2 * q.x * q.x - 2 * q.y * q.y, 0),
		vec4(0, 0, 0, 1)
	);
}

base_func(mat4_t,mat4_transpose)(mat4_t m)
{
	return mat4
	(
		vec4(m.x.x, m.y.x, m.z.x, m.w.x),
		vec4(m.x.y, m.y.y, m.z.y, m.w.y),
		vec4(m.x.z, m.y.z, m.z.z, m.w.z),
		vec4(m.x.w, m.y.w, m.z.w, m.w.w)
	);
}

base_func(mat4_t,mat4_add)(mat4_t l, mat4_t r)
{
	return mat4
	(
		vec4_add(l.x, r.x),
		vec4_add(l.y, r.y),
		vec4_add(l.z, r.z),
		vec4_add(l.w, r.w)
	);
}

base_func(mat4_t,mat4_add_s)(mat4_t m, real_t s)
{
	return mat4
	(
		vec4_add(m.x, vec4(s,s,s,s)),
		vec4_add(m.y, vec4(s,s,s,s)),
		vec4_add(m.z, vec4(s,s,s,s)),
		vec4_add(m.w, vec4(s,s,s,s))
	);
}

base_func(mat4_t,mat4_add_transpose)(mat4_t l, mat4_t r)
{
	return mat4
	(
		vec4(l.x.x + r.x.x, l.x.y + r.y.x, l.x.z + r.z.x, l.x.w + r.w.x),
		vec4(l.y.x + r.x.y, l.y.y + r.y.y, l.y.z + r.z.y, l.y.w + r.w.y),
		vec4(l.z.x + r.x.z, l.z.y + r.y.z, l.z.z + r.z.z, l.z.w + r.w.z),
		vec4(l.w.x + r.x.w, l.w.y + r.y.w, l.w.z + r.z.w, l.w.w + r.w.w)
	);
}

base_func(mat4_t,mat4_sub)(mat4_t l, mat4_t r)
{
	return mat4
	(
		vec4_sub(l.x, r.x),
		vec4_sub(l.y, r.y),
		vec4_sub(l.z, r.z),
		vec4_sub(l.w, r.w)
	);
}

base_func(mat4_t,mat4_sub_s)(mat4_t m, real_t s)
{
	return mat4
	(
		vec4_sub(m.x, vec4(s,s,s,s)),
		vec4_sub(m.y, vec4(s,s,s,s)),
		vec4_sub(m.z, vec4(s,s,s,s)),
		vec4_sub(m.w, vec4(s,s,s,s))
	);
}

base_func(mat4_t,mat4_sub_transpose)(mat4_t l, mat4_t r)
{
	return mat4
	(
		vec4(l.x.x - r.x.x, l.x.y - r.y.x, l.x.z - r.z.x, l.x.w - r.w.x),
		vec4(l.y.x - r.x.y, l.y.y - r.y.y, l.y.z - r.z.y, l.y.w - r.w.y),
		vec4(l.z.x - r.x.z, l.z.y - r.y.z, l.z.z - r.z.z, l.z.w - r.w.z),
		vec4(l.w.x - r.x.w, l.w.y - r.y.w, l.w.z - r.z.w, l.w.w - r.w.w)
	);
}

base_func(mat4_t,mat4_mul)(mat4_t l, mat4_t r)
{
	return mat4
	(
		vec4_mul_mat4(l.x, r),
		vec4_mul_mat4(l.y, r),
		vec4_mul_mat4(l.z, r),
		vec4_mul_mat4(l.w, r)
	);
}

base_func(mat4_t,mat4_mul_s)(mat4_t m, real_t s)
{
	return mat4
	(
		vec4_scale(m.x, s),
		vec4_scale(m.y, s),
		vec4_scale(m.z, s),
		vec4_scale(m.w, s)
	);
}

base_func(mat4_t,mat4_mul_transpose)(mat4_t l, mat4_t r)
{
	return mat4
	(
		vec4_mul_mat4_transpose(l.x, r),
		vec4_mul_mat4_transpose(l.y, r),
		vec4_mul_mat4_transpose(l.z, r),
		vec4_mul_mat4_transpose(l.w, r)
	);
}


