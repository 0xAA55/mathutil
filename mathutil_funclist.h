
#if defined(_MATHUTIL_H_)

#ifndef math_func
#  define math_func(r,n,arg) r n arg
#endif

math_func(real_t, r_rnd, (uint32_t *p_seed));
math_func(real_t, r_sin, (real_t x));
math_func(real_t, r_cos, (real_t x));
math_func(real_t, r_tan, (real_t x));
math_func(real_t, r_abs, (real_t x));
math_func(real_t, r_sgn, (real_t x));
math_func(real_t, r_sqr, (real_t x));
math_func(real_t, r_floor, (real_t x));
math_func(real_t, r_ceil, (real_t x));
math_func(real_t, r_atan, (real_t x));
math_func(real_t, r_exp, (real_t x));
math_func(real_t, r_log, (real_t x));
math_func(real_t, r_pow, (real_t x, real_t y));
math_func(real_t, r_mod, (real_t x, real_t y));
math_func(real_t, r_max, (real_t a, real_t b));
math_func(real_t, r_min, (real_t a, real_t b));
math_func(real_t, r_atan2, (real_t y, real_t x));
math_func(real_t, r_clamp, (real_t n, real_t min_, real_t max_));
math_func(real_t, r_lerp, (real_t a, real_t b, real_t s));
math_func(real_t, r_hermite, (real_t s));
math_func(real_t, r_slerp, (real_t a, real_t b, real_t s));

math_func(vec4_t, vec4, (real_t x, real_t y, real_t z, real_t w));
math_func(vec4_t, vec4_flushcomp, (vec4_t v));
math_func(vec4_t, vec4_abs, (vec4_t v));
math_func(vec4_t, vec4_sgn, (vec4_t v));
math_func(vec4_t, vec4_invert, (vec4_t v));
math_func(real_t, vec4_length, (vec4_t v));
math_func(vec4_t, vec4_normalize, (vec4_t v));
math_func(vec4_t, vec4_scale, (vec4_t v, real_t s));
math_func(vec4_t, vec4_pow, (vec4_t v, real_t n));
math_func(real_t, vec4_dot, (vec4_t v1, vec4_t v2));
math_func(vec4_t, vec4_add, (vec4_t v1, vec4_t v2));
math_func(vec4_t, vec4_sub, (vec4_t v1, vec4_t v2));
math_func(vec4_t, vec4_mul, (vec4_t v1, vec4_t v2));
math_func(vec4_t, vec4_div, (vec4_t v1, vec4_t v2));
math_func(vec4_t, vec4_min, (vec4_t v1, vec4_t v2));
math_func(vec4_t, vec4_max, (vec4_t v1, vec4_t v2));
math_func(vec4_t, vec4_cross3, (vec4_t v1, vec4_t v2));
math_func(vec4_t, vec4_clamp, (vec4_t v, real_t min_, real_t max_));
math_func(vec4_t, vec4_rot_quat, (vec4_t v, quat_t q));
math_func(vec4_t, vec4_mul_mat4, (vec4_t v, mat4_t m));
math_func(vec4_t, vec4_mul_mat4_transpose, (vec4_t v, mat4_t m));
math_func(vec4_t, vec4_lerp, (vec4_t v1, vec4_t v2, real_t s));
math_func(vec4_t, vec4_slerp, (vec4_t v1, vec4_t v2, real_t s));

math_func(quat_t, quat, (real_t x, real_t y, real_t z, real_t w));
math_func(quat_t, quat_flushcomp, (quat_t q));
math_func(quat_t, quat_rot_axis, (vec4_t axis, real_t angle));
math_func(quat_t, quat_mul, (quat_t q1, quat_t q2));
math_func(quat_t, quat_add_vec, (quat_t q, vec4_t v, real_t s));

math_func(mat4_t, mat4, (vec4_t mx, vec4_t my, vec4_t mz, vec4_t mw));
math_func(mat4_t, mat4_flushcomp, (mat4_t m));
math_func(mat4_t, mat4_rot_x, (real_t angle));
math_func(mat4_t, mat4_rot_y, (real_t angle));
math_func(mat4_t, mat4_rot_z, (real_t angle));
math_func(mat4_t, mat4_rot_axis, (vec4_t axis, real_t angle));
math_func(mat4_t, mat4_rot_euler, (real_t yaw, real_t pitch, real_t roll));
math_func(mat4_t, mat4_from_quat, (quat_t q));
math_func(mat4_t, mat4_from_quat_transpose, (quat_t q));
math_func(mat4_t, mat4_transpose, (mat4_t m));
math_func(mat4_t, mat4_add, (mat4_t l, mat4_t r));
math_func(mat4_t, mat4_add_s, (mat4_t m, real_t s));
math_func(mat4_t, mat4_add_transpose, (mat4_t l, mat4_t r));
math_func(mat4_t, mat4_sub, (mat4_t l, mat4_t r));
math_func(mat4_t, mat4_sub_s, (mat4_t m, real_t s));
math_func(mat4_t, mat4_sub_transpose, (mat4_t l, mat4_t r));
math_func(mat4_t, mat4_mul, (mat4_t l, mat4_t r));
math_func(mat4_t, mat4_mul_s, (mat4_t m, real_t s));
math_func(mat4_t, mat4_mul_transpose, (mat4_t l, mat4_t r));

#endif // !_MATHUTIL_CONF_H_
