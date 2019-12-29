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
#define ssse3_func(r,n) r n ## _ssse3
#include"mathutil_ssse3.h"
#endif

#if !MATHUTIL_USE_DOUBLE

#if MATHUTIL_DETECT_CPU
int mathutil_ssse3_implements()
{
	if(!CPUID_SSSE3()) return 0;

	// Boo

	return 1;
}
#endif // MATHUTIL_DETECT_CPU

#else // MATHUTIL_USE_DOUBLE

#if MATHUTIL_DETECT_CPU
int mathutil_ssse3_implements()
{
	if(!CPUID_SSSE3()) return 0;
	
	// Boo

	return 1;
}
#endif // MATHUTIL_DETECT_CPU

#endif // !MATHUTIL_USE_DOUBLE
