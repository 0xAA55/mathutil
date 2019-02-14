#define _CRT_SECURE_NO_WARNINGS

#include"mathutil_bench.h"
#include"mathutil.h"
#include"mathutil_ref.h"
#include"mathutil_sse.h"
#include"mathutil_sse2.h"
#include"mathutil_sse3.h"
#include"mathutil_ssse3.h"
#include"mathutil_sse41.h"
#include"rttimer.h"
#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<time.h>

#if MATHUTIL_USE_DOUBLE
#define type_of_real "double"
#else // !MATHUTIL_USE_DOUBLE
#define type_of_real "float"
#endif

#if _WIN32
#  include<Windows.h>
#  if defined(_WIN64)
#  define bench_suffix "_x64_" type_of_real
#  else // !defined(_WIN64)
#  define bench_suffix "_x86_" type_of_real
#  endif // !defined(_WIN64)
#elif defined(__linux__)
#  define _GNU_SOURCE
#  include<sched.h>
#  include<pthread.h>
#  define bench_suffix "_" type_of_real
#endif

#define BENCH_DURATION 0.05
#define BENCH_VAR ALIGNED_(16)

static FILE *bench_fp = NULL;

static void bench_printf_init()
{
#if _WIN32
	AllocConsole();

	freopen("CONIN$", "r", stdin);
	freopen("CONOUT$", "w", stdout);
	freopen("CONOUT$", "w", stderr);

#endif // _WIN32
	bench_fp = fopen("mathutil_bench" bench_suffix ".txt","w");
}

static void bench_printf_cleanup()
{
#if _WIN32
	FreeConsole();
#endif // _WIN32
	if(bench_fp) fclose(bench_fp);
	bench_fp = NULL;
}

void bench_printf(const char *fmt, ...)
{
	va_list ap;
	va_start(ap, fmt);

	vprintf(fmt, ap);
	if(bench_fp) vfprintf(bench_fp, fmt, ap);

	va_end(ap);
}

static real_t bench_rnd(uint32_t *p_seed)
{
	return r_rnd(p_seed);
}

static uint32_t bench_init_rs = 0;

static void bench_bind_cpu()
{
#if _WIN32
	if(SetThreadAffinityMask(GetCurrentThread(), 1))
	{
		bench_printf("SetThreadAffinityMask(GetCurrentThread(), 1) success.\n");
	}
	else
	{
		bench_printf("[WARN] SetThreadAffinityMask(GetCurrentThread(), 1) failed: GetLastError() = %u\n",
			GetLastError());
	}
#elif defined(__linux__)
	int s, j;
	unsigned num_cpu;
	cpu_set_t cpuset;
	pthread_t thread;

	thread = pthread_self();
	CPU_ZERO(&cpuset);
	CPU_SET(0, &cpuset);

	s = pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
	if(s) bench_printf("[WARN] pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset) failed: returned = %d\n", s);

	s = pthread_getaffinity_np(thread, sizeof(cpu_set_t), &cpuset);
	if(s) bench_printf("[WARN] pthread_getaffinity_np(thread, sizeof(cpu_set_t), &cpuset) failed: returned = %d\n", s);

	num_cpu = 0;
	for (j = 0; j < CPU_SETSIZE; j++)
	{
        if(CPU_ISSET(j, &cpuset)) num_cpu++;
	}
	if(num_cpu != 1) bench_printf("[WARN] There's %u CPUs in the cpuset.\n", num_cpu);
#endif
}

static void bench_init()
{
	bench_printf_init();
	bench_bind_cpu();
	bench_init_rs = (uint32_t)time(NULL) % 10000;
}

static void bench_finish()
{
	bench_printf_cleanup();
}

#if MATHUTIL_DETECT_SIMD
#  define math_func(r_hint,r,f,arg,carg) static r(*f ## _fp)arg = f ## _ref; static r r_hint ## _ ## f; static double f ## _score = 0;  static double f ## _score_ref = 0; static const char *prev_suffix_ ## f = "ref";
#  include"mathutil_funclist.h"
#  undef math_func
#else // !MATHUTIL_DETECT_SIMD
#  define math_func(r_hint,r,f,arg,carg) static double f ## _score = 0;
#  include"mathutil_funclist.h"
#  undef math_func
#endif // !MATHUTIL_DETECT_SIMD

#define test_batch 0x100000

#if MATHUTIL_DETECT_SIMD
#define sr_cmp(sr,f) \
	{ \
		if(r_abs(sr - sr_ ## f) >= r_epsilon) \
		{ \
			bench_printf("[WARN] Scalar result is different\n" \
				"\t\tnew: %lf\n" \
				"\t\told: %lf\n", \
				(double)sr, \
				(double)sr_ ## f); \
		} \
	}

static int vec4_is_same(vec4_t v1, vec4_t v2)
{
	return
		r_abs(v1.x - v2.x) < r_epsilon &&
		r_abs(v1.y - v2.y) < r_epsilon &&
		r_abs(v1.z - v2.z) < r_epsilon &&
		r_abs(v1.w - v2.w) < r_epsilon;
}

static int quat_is_same(quat_t q1, quat_t q2)
{
	return
		r_abs(q1.x - q2.x) < r_epsilon &&
		r_abs(q1.y - q2.y) < r_epsilon &&
		r_abs(q1.z - q2.z) < r_epsilon &&
		r_abs(q1.w - q2.w) < r_epsilon;
}

static int mat4_is_same(mat4_t m1, mat4_t m2)
{
	return
		vec4_is_same(m1.x, m2.x) &&
		vec4_is_same(m1.y, m2.y) &&
		vec4_is_same(m1.z, m2.z) &&
		vec4_is_same(m1.w, m2.w);
}

#define vr_cmp(vr,f) \
	{ \
		if(!vec4_is_same(vr, vr_ ## f)) \
		{ \
			bench_printf("[WARN] Vector result is different\n" \
				"\t\tnew: (%lf, %lf, %lf, %lf)\n" \
				"\t\told: (%lf, %lf, %lf, %lf)\n", \
				(double)vr.x, (double)vr.y, (double)vr.z, (double)vr.w, \
				(double)vr_ ## f.x, (double)vr_ ## f.y, (double)vr_ ## f.z, (double)vr_ ## f.w); \
		} \
	}

#define qr_cmp(qr,f) \
	{ \
		if(!quat_is_same(qr, qr_ ## f)) \
		{ \
			bench_printf("[WARN] Quaternion result is different\n" \
				"\t\tnew: (%lf, %lf, %lf, %lf)\n" \
				"\t\told: (%lf, %lf, %lf, %lf)\n", \
				(double)qr.x, (double)qr.y, (double)qr.z, (double)qr.w, \
				(double)qr_ ## f.x, (double)qr_ ## f.y, (double)qr_ ## f.z, (double)qr_ ## f.w); \
		} \
	}

#define mr_cmp(mr,f) \
	{ \
		if(!mat4_is_same(mr, mr_ ## f)) \
		{ \
			bench_printf("[WARN] Vector result is different\n" \
				"\t\tnew:\n" \
				"\t\t    %lf, %lf, %lf, %lf\n" \
				"\t\t    %lf, %lf, %lf, %lf\n" \
				"\t\t    %lf, %lf, %lf, %lf\n" \
				"\t\t    %lf, %lf, %lf, %lf\n" \
				"\t\told:\n" \
				"\t\t    %lf, %lf, %lf, %lf\n" \
				"\t\t    %lf, %lf, %lf, %lf\n" \
				"\t\t    %lf, %lf, %lf, %lf\n" \
				"\t\t    %lf, %lf, %lf, %lf\n", \
				(double)mr.x.x, (double)mr.x.y, (double)mr.x.z, (double)mr.x.w, \
				(double)mr.y.x, (double)mr.y.y, (double)mr.y.z, (double)mr.y.w, \
				(double)mr.z.x, (double)mr.z.y, (double)mr.z.z, (double)mr.z.w, \
				(double)mr.w.x, (double)mr.w.y, (double)mr.w.z, (double)mr.w.w, \
				(double)mr_ ## f.x.x, (double)mr_ ## f.x.y, (double)mr_ ## f.x.z, (double)mr_ ## f.x.w, \
				(double)mr_ ## f.y.x, (double)mr_ ## f.y.y, (double)mr_ ## f.y.z, (double)mr_ ## f.y.w, \
				(double)mr_ ## f.z.x, (double)mr_ ## f.z.y, (double)mr_ ## f.z.z, (double)mr_ ## f.z.w, \
				(double)mr_ ## f.w.x, (double)mr_ ## f.w.y, (double)mr_ ## f.w.z, (double)mr_ ## f.w.w); \
		} \
	}

#define test_function(r,f,carg) if(NULL != (void*)f) \
	{ \
		uint64_t times = 0; \
		unsigned batch; \
		rttimer_t tmr; \
		double tv; \
		bench_printf("Benchmarking " # f # carg ";\n"); \
		rttimer_init(&tmr, 1); \
		rttimer_start(&tmr); \
		do \
		{ \
			for(batch = 0; batch < test_batch; batch++) r = f carg;  \
			times += test_batch; \
			tv = rttimer_gettime(&tmr); \
		} while(tv < BENCH_DURATION); \
		score = times / tv; \
		f ## _score = score; \
		print_result(times, tv); \
		f ## _fp = f; \
	} \
	else \
	{ \
		bench_printf("Function pointer of " # f "() is a NULL pointer!\n"); \
	}

#define test_func_ref(r,f,carg) { \
		double score = 0; \
		test_function(r,f,carg); \
		r ## _ ## f = r; \
		f ## _score_ref = score; \
	}

#define test_func_if_overrided(r,f,carg) if(f ## _fp != f) \
	{ \
		double score = 0; \
		double prev = f ## _score; \
		test_function(r,f,carg); \
		bench_printf("\tCompare to ref: %lf %%\n", 100 * f ## _score / f ## _score_ref); \
		if(f ## _score < f ## _score_ref)  bench_printf("[SLOWER THAN REF] (cur: " # f "%s).\n",suffix); \
		if(strcmp(prev_suffix_ ## f, "ref")) \
		{ \
			bench_printf("\tCompare to prev: %lf %%\n", 100 * f ## _score / prev); \
			if(f ## _score < prev)  bench_printf("[SLOWER THAN PREV] (prev: " # f "%s).\n", prev_suffix_ ## f); \
		} \
		r ## _cmp(r,f); \
		prev_suffix_ ## f = suffix; \
	}
	
#else // !MATHUTIL_DETECT_SIMD
	#define test_function(r,f,carg) if(NULL != (void*)f) \
	{ \
		uint64_t times = 0; \
		unsigned batch; \
		rttimer_t tmr; \
		double tv; \
		bench_printf("Benchmarking " # f # carg ";\n"); \
		rttimer_init(&tmr, 1); \
		rttimer_start(&tmr); \
		do \
		{ \
			for(batch = 0; batch < test_batch; batch++) r = f carg;  \
			times += test_batch; \
			tv = rttimer_gettime(&tmr); \
		} while(tv < BENCH_DURATION); \
		f ## _score = times / tv; \
		print_result(times, tv); \
	} \
	else \
	{ \
		bench_printf("Function pointer of " # f "() is a NULL pointer!\n"); \
	}
	
#endif // !MATHUTIL_DETECT_SIMD

static void print_result(uint64_t count, double cost)
{
	bench_printf("\tTime cost: %lf\n", cost);
	if(count > test_batch)
		bench_printf("\tRun count: %llu (%llu rounds)\n", count, count / test_batch);
	else
		bench_printf("\tRun count: %llu\n", count);
	bench_printf("\tAverage freq: %lf times per us\n", ((double)count / cost) * 0.000001);
}

static void test_all_funcs()
{
	uint32_t tmp_rs = bench_init_rs;
	BENCH_VAR uint32_t rs = 1;
	uint32_t *p_seed = &rs;
	BENCH_VAR real_t x = bench_rnd(&tmp_rs) * 100;
	BENCH_VAR real_t y = bench_rnd(&tmp_rs) * 100;
	BENCH_VAR real_t z = bench_rnd(&tmp_rs) * 100;
	BENCH_VAR real_t w = bench_rnd(&tmp_rs) * 100;
	BENCH_VAR real_t a = bench_rnd(&tmp_rs) * 100;
	BENCH_VAR real_t b = bench_rnd(&tmp_rs) * 100;
	BENCH_VAR real_t n = bench_rnd(&tmp_rs) * 100;
	BENCH_VAR real_t s = bench_rnd(&tmp_rs);
	BENCH_VAR real_t min_ = -bench_rnd(&tmp_rs);
	BENCH_VAR real_t max_ = bench_rnd(&tmp_rs);
	BENCH_VAR vec4_t v = {bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10};
	BENCH_VAR vec4_t mx = {bench_rnd(&tmp_rs),bench_rnd(&tmp_rs),bench_rnd(&tmp_rs),bench_rnd(&tmp_rs)};
	BENCH_VAR vec4_t my = {bench_rnd(&tmp_rs),bench_rnd(&tmp_rs),bench_rnd(&tmp_rs),bench_rnd(&tmp_rs)};
	BENCH_VAR vec4_t mz = {bench_rnd(&tmp_rs),bench_rnd(&tmp_rs),bench_rnd(&tmp_rs),bench_rnd(&tmp_rs)};
	BENCH_VAR vec4_t mw = {bench_rnd(&tmp_rs),bench_rnd(&tmp_rs),bench_rnd(&tmp_rs),bench_rnd(&tmp_rs)};
	BENCH_VAR vec4_t v1 = {bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10};
	BENCH_VAR vec4_t v2 = {bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10};
	BENCH_VAR quat_t q = {0,0,0,1};
	BENCH_VAR quat_t q1 = {0,0,0,1};
	BENCH_VAR quat_t q2 = {0,0,0,1};
	BENCH_VAR mat4_t m = {0};
	BENCH_VAR mat4_t l = {0};
	BENCH_VAR mat4_t r = {0};
	BENCH_VAR real_t angle = r_pi;
	BENCH_VAR real_t yaw = r_pi * (real_t)0.123;
	BENCH_VAR real_t pitch = r_pi * (real_t)0.456;
	BENCH_VAR real_t roll = r_pi * (real_t)0.789;
	BENCH_VAR vec4_t axis = {0,1,0,0};
	BENCH_VAR real_t sr = 0;
	BENCH_VAR vec4_t vr = {0,0,0,0};
	BENCH_VAR quat_t qr = {0,0,0,1};
	BENCH_VAR mat4_t mr = {0};
	
	(void)sr;
	(void)vr;
	(void)qr;
	(void)mr;
	
	q = quat_rot_axis(vec4(0,0,1,0), r_pi * bench_rnd(&tmp_rs));
	q1 = quat_rot_axis(vec4(0,1,0,0), r_pi * bench_rnd(&tmp_rs));
	q2 = quat_rot_axis(vec4(1,0,0,0), r_pi * bench_rnd(&tmp_rs));
	
#if MATHUTIL_DETECT_SIMD
#	define math_func(r_hint,r,n,arg,carg) test_func_ref(r_hint,n,carg)
#	include"mathutil_funclist.h"
#	undef math_func
#else // !MATHUTIL_DETECT_SIMD
#	define math_func(r_hint,r,n,arg,carg) test_function(r_hint,n,carg)
#	include"mathutil_funclist.h"
#	undef math_func
#endif // !MATHUTIL_DETECT_SIMD
}

#if MATHUTIL_DETECT_SIMD
static void test_overrided_funcs(const char*suffix)
{
	uint32_t tmp_rs = bench_init_rs;
	BENCH_VAR uint32_t rs = 1;
	uint32_t *p_seed = &rs;
	BENCH_VAR real_t x = bench_rnd(&tmp_rs) * 100;
	BENCH_VAR real_t y = bench_rnd(&tmp_rs) * 100;
	BENCH_VAR real_t z = bench_rnd(&tmp_rs) * 100;
	BENCH_VAR real_t w = bench_rnd(&tmp_rs) * 100;
	BENCH_VAR real_t a = bench_rnd(&tmp_rs) * 100;
	BENCH_VAR real_t b = bench_rnd(&tmp_rs) * 100;
	BENCH_VAR real_t n = bench_rnd(&tmp_rs) * 100;
	BENCH_VAR real_t s = bench_rnd(&tmp_rs);
	BENCH_VAR real_t min_ = -bench_rnd(&tmp_rs);
	BENCH_VAR real_t max_ = bench_rnd(&tmp_rs);
	BENCH_VAR vec4_t v = {bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10};
	BENCH_VAR vec4_t mx = {bench_rnd(&tmp_rs),bench_rnd(&tmp_rs),bench_rnd(&tmp_rs),bench_rnd(&tmp_rs)};
	BENCH_VAR vec4_t my = {bench_rnd(&tmp_rs),bench_rnd(&tmp_rs),bench_rnd(&tmp_rs),bench_rnd(&tmp_rs)};
	BENCH_VAR vec4_t mz = {bench_rnd(&tmp_rs),bench_rnd(&tmp_rs),bench_rnd(&tmp_rs),bench_rnd(&tmp_rs)};
	BENCH_VAR vec4_t mw = {bench_rnd(&tmp_rs),bench_rnd(&tmp_rs),bench_rnd(&tmp_rs),bench_rnd(&tmp_rs)};
	BENCH_VAR vec4_t v1 = {bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10};
	BENCH_VAR vec4_t v2 = {bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10, bench_rnd(&tmp_rs) * 10};
	BENCH_VAR quat_t q = {0,0,0,1};
	BENCH_VAR quat_t q1 = {0,0,0,1};
	BENCH_VAR quat_t q2 = {0,0,0,1};
	BENCH_VAR mat4_t m = {0};
	BENCH_VAR mat4_t l = {0};
	BENCH_VAR mat4_t r = {0};
	BENCH_VAR real_t angle = r_pi;
	BENCH_VAR real_t yaw = r_pi * (real_t)0.123;
	BENCH_VAR real_t pitch = r_pi * (real_t)0.456;
	BENCH_VAR real_t roll = r_pi * (real_t)0.789;
	BENCH_VAR vec4_t axis = {0,1,0,0};
	BENCH_VAR real_t sr = 0;
	BENCH_VAR vec4_t vr = {0,0,0,0};
	BENCH_VAR quat_t qr = {0,0,0,1};
	BENCH_VAR mat4_t mr = {0};
	
	q = quat_rot_axis(vec4(0,0,1,0), r_pi * bench_rnd(&tmp_rs));
	q1 = quat_rot_axis(vec4(0,1,0,0), r_pi * bench_rnd(&tmp_rs));
	q2 = quat_rot_axis(vec4(1,0,0,0), r_pi * bench_rnd(&tmp_rs));

#	define math_func(r_hint,r,n,arg,carg) test_func_if_overrided(r_hint,n,carg)
#	include"mathutil_funclist.h"
#	undef math_func
}
#endif // MATHUTIL_DETECT_SIMD

#if MATHUTIL_DETECT_SIMD
void mathutil_ref_implements();

void mathutil_bench()
{
	bench_init();

	mathutil_ref_implements();
	bench_printf(
		"Benchmarking for base implements\n"
		"================================\n");
	test_all_funcs();
	
	mathutil_sse_implements();
	bench_printf(
		"Benchmarking for SSE implements\n"
		"===============================\n");
	test_overrided_funcs("_sse");
	
	mathutil_sse2_implements();
	bench_printf(
		"Benchmarking for SSE2 implements\n"
		"===============================\n");
	test_overrided_funcs("_sse2");
	
	mathutil_sse3_implements();
	bench_printf(
		"Benchmarking for SSE3 implements\n"
		"===============================\n");
	test_overrided_funcs("_sse3");
	
	mathutil_ssse3_implements();
	bench_printf(
		"Benchmarking for SSSE3 implements\n"
		"===============================\n");
	test_overrided_funcs("_ssse3");
	
	mathutil_sse41_implements();
	bench_printf(
		"Benchmarking for SSE4.1 implements\n"
		"===============================\n");
	test_overrided_funcs("_sse41");
	
	bench_finish();
}


#else // !MATHUTIL_DETECT_SIMD

void mathutil_bench()
{
	bench_init();
	test_all_funcs();
	bench_finish();
}

#endif // !MATHUTIL_DETECT_SIMD

int main()
{
	mathutil_bench();
	
	return EXIT_SUCCESS;
}

