#include"cpudetect.h"

#if MATHUTIL_DETECT_CPU

//  Misc.
static int HW_MMX = 0;
static int HW_x64 = 0;
static int HW_ABM = 0;      // Advanced Bit Manipulation
static int HW_RDRAND = 0;
static int HW_BMI1 = 0;
static int HW_BMI2 = 0;
static int HW_ADX = 0;
static int HW_PREFETCHWT1 = 0;

//  SIMD: 128-bit
static int HW_SSE = 0;
static int HW_SSE2 = 0;
static int HW_SSE3 = 0;
static int HW_SSSE3 = 0;
static int HW_SSE41 = 0;
static int HW_SSE42 = 0;
static int HW_SSE4a = 0;
static int HW_AES = 0;
static int HW_SHA = 0;

//  SIMD: 256-bit
static int HW_AVX = 0;
static int HW_XOP = 0;
static int HW_FMA3 = 0;
static int HW_FMA4 = 0;
static int HW_AVX2 = 0;

//  SIMD: 512-bit
static int HW_AVX512F = 0;    //  AVX512 Foundation
static int HW_AVX512CD = 0;   //  AVX512 Conflict Detection
static int HW_AVX512PF = 0;   //  AVX512 Prefetch
static int HW_AVX512ER = 0;   //  AVX512 Exponential + Reciprocal
static int HW_AVX512VL = 0;   //  AVX512 Vector Length Extensions
static int HW_AVX512BW = 0;   //  AVX512 Byte + Word
static int HW_AVX512DQ = 0;   //  AVX512 Doubleword + Quadword
static int HW_AVX512IFMA = 0; //  AVX512 Integer 52-bit Fused Multiply-Add
static int HW_AVX512VBMI = 0; //  AVX512 Vector Byte Manipulation Instructions

static int _OS_x64 = 0;
static int _OS_AVX = 0;
static int _OS_AVX512 = 0;

static void cpudetect_init();

#define _queryfunc_hw(f) int CPUID_ ## f(){cpudetect_init(); return HW_ ## f;}
#define _queryfunc_os(f) int OS_ ## f(){cpudetect_init(); return _OS_ ## f;}

_queryfunc_hw(MMX);
_queryfunc_hw(x64);
_queryfunc_hw(ABM);
_queryfunc_hw(RDRAND);
_queryfunc_hw(BMI1);
_queryfunc_hw(BMI2);
_queryfunc_hw(ADX);
_queryfunc_hw(PREFETCHWT1);

_queryfunc_hw(SSE);
_queryfunc_hw(SSE2);
_queryfunc_hw(SSE3);
_queryfunc_hw(SSSE3);
_queryfunc_hw(SSE41);
_queryfunc_hw(SSE42);
_queryfunc_hw(SSE4a);
_queryfunc_hw(AES);
_queryfunc_hw(SHA);

_queryfunc_hw(AVX);
_queryfunc_hw(XOP);
_queryfunc_hw(FMA3);
_queryfunc_hw(FMA4);
_queryfunc_hw(AVX2);

_queryfunc_hw(AVX512F);
_queryfunc_hw(AVX512CD);
_queryfunc_hw(AVX512PF);
_queryfunc_hw(AVX512ER);
_queryfunc_hw(AVX512VL);
_queryfunc_hw(AVX512BW);
_queryfunc_hw(AVX512DQ);
_queryfunc_hw(AVX512IFMA);
_queryfunc_hw(AVX512VBMI);

_queryfunc_os(x64);
_queryfunc_os(AVX);
_queryfunc_os(AVX512);

#include<stdint.h>
#include<signal.h>

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
#if _WIN32
#include<Windows.h>
#include<intrin.h>
#define cpuid(info, x) __cpuidex(info, x, 0)

#if _WIN32 && !defined(_M_X64)
static int IsWow64()
{
	BOOL (WINAPI *LPFN_ISWOW64PROCESS)(HANDLE, PBOOL);
	BOOL bIsWow64 = FALSE;

	LPFN_ISWOW64PROCESS = (BOOL(WINAPI*)(HANDLE, PBOOL))GetProcAddress(GetModuleHandle(TEXT("kernel32.dll")), "IsWow64Process");
	if(!LPFN_ISWOW64PROCESS) return 0;
	if(!LPFN_ISWOW64PROCESS(GetCurrentProcess(), &bIsWow64))
	{
		bIsWow64 = FALSE;
	}
	return bIsWow64 != FALSE;
}
#endif

#elif defined(__GNUC__) || defined(_clang_)
#include<cpuid.h>
static void cpuid(int info[4], int InfoType)
{
    __cpuid_count(InfoType, 0, info[0], info[1], info[2], info[3]);
}
#endif // !(defined(__GNUC__) || defined(_clang_))

#ifndef _XCR_XFEATURE_ENABLED_MASK
#define _XCR_XFEATURE_ENABLED_MASK 0
#endif

static size_t g_num_ids = 0;
static size_t g_num_exids = 0;

static int g_cpudetect_initialized = 0;

static void cpudetect_OS_AVX()
{
    //  Copied from: http://stackoverflow.com/a/22521619/922184
	int info[4];
	int osUsesXSAVE_XRSTORE;
	int cpuAVXSuport;

	_OS_AVX = 0;

    cpuid(info, 1);

    osUsesXSAVE_XRSTORE = (info[2] & (1 << 27)) != 0;
    cpuAVXSuport = (info[2] & (1 << 28)) != 0;

    if (osUsesXSAVE_XRSTORE && cpuAVXSuport)
    {
        uint64_t xcrFeatureMask = _xgetbv(_XCR_XFEATURE_ENABLED_MASK);
        _OS_AVX = (xcrFeatureMask & 0x6) == 0x6;
    }
}

static void cpudetect_OS_AVX512()
{
	uint64_t xcrFeatureMask;

	_OS_AVX512 = 0;
    if(!_OS_AVX) return;

    xcrFeatureMask = _xgetbv(_XCR_XFEATURE_ENABLED_MASK);
    _OS_AVX512 = (xcrFeatureMask & 0xe6) == 0xe6;
}

static void cpudetect_init_x86()
{
	int info[4];

	cpuid(info, 0);
	g_num_ids = info[0];

	cpuid(info, 0x80000000);
	g_num_exids = info[0];

	//  Detect Features
	if (g_num_ids >= 0x00000001)
	{
		cpuid(info,0x00000001);
		HW_MMX    = (info[3] & ((int)1 << 23)) != 0;
		HW_SSE    = (info[3] & ((int)1 << 25)) != 0;
		HW_SSE2   = (info[3] & ((int)1 << 26)) != 0;
		HW_SSE3   = (info[2] & ((int)1 <<  0)) != 0;

		HW_SSSE3  = (info[2] & ((int)1 <<  9)) != 0;
		HW_SSE41  = (info[2] & ((int)1 << 19)) != 0;
		HW_SSE42  = (info[2] & ((int)1 << 20)) != 0;
		HW_AES    = (info[2] & ((int)1 << 25)) != 0;

		HW_AVX    = (info[2] & ((int)1 << 28)) != 0;
		HW_FMA3   = (info[2] & ((int)1 << 12)) != 0;

		HW_RDRAND = (info[2] & ((int)1 << 30)) != 0;
	}
	if (g_num_ids >= 0x00000007)
	{
		cpuid(info,0x00000007);
		HW_AVX2   = (info[1] & ((int)1 <<  5)) != 0;

		HW_BMI1        = (info[1] & ((int)1 <<  3)) != 0;
		HW_BMI2        = (info[1] & ((int)1 <<  8)) != 0;
		HW_ADX         = (info[1] & ((int)1 << 19)) != 0;
		HW_SHA         = (info[1] & ((int)1 << 29)) != 0;
		HW_PREFETCHWT1 = (info[2] & ((int)1 <<  0)) != 0;

		HW_AVX512F     = (info[1] & ((int)1 << 16)) != 0;
		HW_AVX512CD    = (info[1] & ((int)1 << 28)) != 0;
		HW_AVX512PF    = (info[1] & ((int)1 << 26)) != 0;
		HW_AVX512ER    = (info[1] & ((int)1 << 27)) != 0;
		HW_AVX512VL    = (info[1] & ((int)1 << 31)) != 0;
		HW_AVX512BW    = (info[1] & ((int)1 << 30)) != 0;
		HW_AVX512DQ    = (info[1] & ((int)1 << 17)) != 0;
		HW_AVX512IFMA  = (info[1] & ((int)1 << 21)) != 0;
		HW_AVX512VBMI  = (info[2] & ((int)1 <<  1)) != 0;
	}
	if (g_num_exids >= 0x80000001){
		cpuid(info,0x80000001);
		HW_x64   = (info[3] & ((int)1 << 29)) != 0;
		HW_ABM   = (info[2] & ((int)1 <<  5)) != 0;
		HW_SSE4a = (info[2] & ((int)1 <<  6)) != 0;
		HW_FMA4  = (info[2] & ((int)1 << 16)) != 0;
		HW_XOP   = (info[2] & ((int)1 << 11)) != 0;
	}


#	if _WIN32
#	ifdef _M_X64
	_OS_x64 = 1;
#	else
	_OS_x64 = (IsWow64() != 0);
#	endif
#	else
	_OS_x64 = 1;
#	endif

	cpudetect_OS_AVX();
	cpudetect_OS_AVX512();
}

#endif // !(defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86))

static void cpudetect_init()
{
	if(g_cpudetect_initialized) return;
	
#	if defined(__x86_64__) || defined(_M_X64) || defined(__i386) || defined(_M_IX86)
	cpudetect_init_x86();
#	endif
	
	g_cpudetect_initialized = 1;
}

#endif // !MATHUTIL_DETECT_CPU
