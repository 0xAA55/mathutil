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
#define sse2_func(r,n) r n ## _sse2
#include"mathutil_sse2.h"
#endif

#if !MATHUTIL_USE_DOUBLE

#if MATHUTIL_DETECT_CPU
int mathutil_sse2_implements()
{
	if(!CPUID_SSE2()) return 0;
	
	vec4_abs = vec4_abs_sse2;
	vec4_sgn = vec4_sgn_sse2;
	vec4_invert = vec4_invert_sse2;

	// Boo

	return 1;
}
#endif // MATHUTIL_DETECT_CPU

#else // MATHUTIL_USE_DOUBLE

#if MATHUTIL_DETECT_CPU
int mathutil_sse2_implements()
{
	if(!CPUID_SSE2()) return 0;

	vec4_abs = vec4_abs_sse2;
	vec4_sgn = vec4_sgn_sse2;
	vec4_invert = vec4_invert_sse2;
	vec4_length = vec4_length_sse2;
	vec4_scale = vec4_scale_sse2;
	vec4_clamp = vec4_clamp_sse2;
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
#endif // MATHUTIL_DETECT_CPU

#endif // !MATHUTIL_USE_DOUBLE

