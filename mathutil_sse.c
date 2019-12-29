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

#define MATHUTIL_INTERNAL 1

#include"mathutil.h"
#include"cpudetect.h"

#if MATHUTIL_DETECT_CPU
#define sse_func(r,n) r n ## _sse
#include"mathutil_sse.h"
#endif

#if !MATHUTIL_USE_DOUBLE

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

#if MATHUTIL_DETECT_CPU
int mathutil_sse_implements()
{
	if(!CPUID_SSE()) return 0;

	// Yeah

	return 1;
}
#endif // MATHUTIL_DETECT_CPU

#endif // MATHUTIL_USE_DOUBLE
