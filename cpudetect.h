#ifndef _CPUDETECT_H_
#define _CPUDETECT_H_ 1

#include"mathutil_conf.h"

#if MATHUTIL_DETECT_CPU

int CPUID_MMX();
int CPUID_x64();
int CPUID_ABM();
int CPUID_RDRAND();
int CPUID_BMI1();
int CPUID_BMI2();
int CPUID_ADX();
int CPUID_PREFETCHWT1();

int CPUID_SSE();
int CPUID_SSE2();
int CPUID_SSE3();
int CPUID_SSSE3();
int CPUID_SSE41();
int CPUID_SSE42();
int CPUID_SSE4a();
int CPUID_AES();
int CPUID_SHA();

int CPUID_AVX();
int CPUID_XOP();
int CPUID_FMA3();
int CPUID_FMA4();
int CPUID_AVX2();

int CPUID_AVX512F();
int CPUID_AVX512CD();
int CPUID_AVX512PF();
int CPUID_AVX512ER();
int CPUID_AVX512VL();
int CPUID_AVX512BW();
int CPUID_AVX512DQ();
int CPUID_AVX512IFMA();
int CPUID_AVX512VBMI();

int OS_x64();
int OS_AVX();
int OS_AVX512();

#endif // !MATHUTIL_DETECT_CPU

#endif // !_CPUDETECT_H_
