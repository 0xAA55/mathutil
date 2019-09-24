#include"mathutil.h"
#include"mathutil_ref.h"
#include"mathutil_sse.h"
#include"mathutil_sse2.h"
#include"mathutil_sse3.h"
#include"mathutil_ssse3.h"
#include"mathutil_sse41.h"

#if MATHUTIL_DETECT_CPU
void mathutil_ref_implements()
{
#	define math_func(r,n,arg) n = n ## _ref
#	include"mathutil_funclist.h"
#	undef math_func
}

static int g_mathutil_initialized = 0;

static void mathutil_init()
{
	if(g_mathutil_initialized) return;

	mathutil_ref_implements();
	mathutil_sse_implements();
	mathutil_sse2_implements();
	mathutil_sse3_implements();
	mathutil_ssse3_implements();
	mathutil_sse41_implements();

	g_mathutil_initialized = 1;
}

#	define math_func(r,n,arg) r n ## _first arg{mathutil_init(); return n carg;} r(*n)arg = n ## _first
#	include"mathutil_funclist.h"
#	undef math_func

#endif
