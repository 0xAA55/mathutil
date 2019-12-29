#ifndef _MATHUTIL_SSE_COMMON_H_
#define _MATHUTIL_SSE_COMMON_H_ 1

#include"mathutil_conf.h"
#include"mathutil.h"

#if COMPILER_FLAVOR == 1
#  if !__INTELLISENSE__
#    define _compiler_barrier _ReadWriteBarrier()
#  else
#    define _compiler_barrier
#  endif
#elif COMPILER_FLAVOR == 2
#  define _compiler_barrier asm volatile("" ::: "memory")
#else // COMPILER_FLAVOR == else
#  define _compiler_barrier (void)1
#endif // !COMPILER_FLAVOR

#include<xmmintrin.h>
#include<immintrin.h>

#if COMPILER_FLAVOR != 1
#  define _mm_castps_si128(a) (__m128i)(a)
#  define _mm_castsi128_ps(a) (__m128)(a)
#  define _mm_castpd_si128(a) (__m128i)(a)
#  define _mm_castsi128_pd(a) (__m128d)(a)
#endif

#if MATHUTIL_VAR_NOT_ALIGNED && !MATHUTIL_VAR_ASSUME_ALIGNED
#define _mm_load_ps_a _mm_loadu_ps
#define _mm_load_pd_a _mm_loadu_pd
#define _mm_store_ps_a _mm_storeu_ps
#define _mm_store_pd_a _mm_storeu_pd
#define _mm_load_si128_a _mm_loadu_si128
#define _mm_store_si128_a _mm_storeu_si128
#else // !MATHUTIL_VAR_NOT_ALIGNED || MATHUTIL_VAR_ASSUME_ALIGNED
#define _mm_load_ps_a _mm_load_ps
#define _mm_load_pd_a _mm_load_pd
#define _mm_store_ps_a _mm_store_ps
#define _mm_store_pd_a _mm_store_pd
#define _mm_load_si128_a _mm_load_si128
#define _mm_store_si128_a _mm_store_si128
#endif

#if !MATHUTIL_USE_DOUBLE

#if VEC4_WITH_M128_XYZW
#define vec4_load(v) ((v).m_xyzw)
#define vec4_loadi(v) _mm_castps_si128((v).m_xyzw)
#define quat_load(v) ((v).m_xyzw)
#define vec4_set_result(pv, m) ((pv)->m_xyzw = m)
#define vec4_set_iresult(pv, m) ((pv)->m_xyzw = _mm_castsi128_ps(m))
#define quat_set_result(pq, m) ((pq)->m_xyzw = m)
#else
#define vec4_load(v) _mm_load_ps_a(&((v).x))
#define vec4_loadi(v) _mm_castps_si128(_mm_load_ps_a(&((v).x)))
#define quat_load(v) _mm_load_ps_a(&((v).x))
#define vec4_set_result(pv, m) _mm_store_ps_a(&((pv)->x), m)
#define vec4_set_iresult(pv, m) _mm_store_ps_a(&((pv)->x), _mm_castsi128_ps(m))
#define quat_set_result(pq, m) _mm_store_ps_a(&((pq)->x), m)
#endif

#else // MATHUTIL_USE_DOUBLE

#if VEC4_WITH_M128_XY_ZW
#define vec4_load_xy(v) ((v).m_xy)
#define vec4_load_zw(v) ((v).m_zw)
#define vec4_loadi_xy(v) _mm_castpd_si128((v).m_xy)
#define vec4_loadi_zw(v) _mm_castpd_si128((v).m_zw)
#define vec4_set_result_xy(pv, m) ((pv)->m_xy = m)
#define vec4_set_result_zw(pv, m) ((pv)->m_zw = m)
#define vec4_set_iresult_xy(pv, m) ((pv)->m_xy = _mm_castsi128_pd(m))
#define vec4_set_iresult_zw(pv, m) ((pv)->m_zw = _mm_castsi128_pd(m))
#define quat_load_xy(q) ((q).m_xy)
#define quat_load_zw(q) ((q).m_zw)
#define quat_set_result_xy(pq, m) ((pq)->m_xy = m)
#define quat_set_result_zw(pq, m) ((pq)->m_zw = m)
#else // !VEC4_WITH_M128_XY_ZW
#define vec4_load_xy(v) _mm_load_pd_a(&((v).x))
#define vec4_load_zw(v) _mm_load_pd_a(&((v).z))
#define vec4_loadi_xy(v) _mm_castpd_si128(_mm_load_pd_a(&((v).x)))
#define vec4_loadi_zw(v) _mm_castpd_si128(_mm_load_pd_a(&((v).z)))
#define vec4_set_result_xy(pv, m) _mm_store_pd_a(&((pv)->x), m)
#define vec4_set_result_zw(pv, m) _mm_store_pd_a(&((pv)->z), m)
#define vec4_set_iresult_xy(pv, m) _mm_store_pd_a(&((pv)->x), _mm_castsi128_pd(m))
#define vec4_set_iresult_zw(pv, m) _mm_store_pd_a(&((pv)->z), _mm_castsi128_pd(m))
#define quat_load_xy(q) _mm_load_pd_a(&((q).x))
#define quat_load_zw(q) _mm_load_pd_a(&((q).z))
#define quat_set_result_xy(pq, m) _mm_store_pd_a(&((pq)->x), m)
#define quat_set_result_zw(pq, m) _mm_store_pd_a(&((pq)->z), m)
#endif // !VEC4_WITH_M128_XY_ZW

#endif // MATHUTIL_USE_DOUBLE

#endif // _MATHUTIL_SSE_COMMON_H_
