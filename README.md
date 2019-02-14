# mathutil
A utility for 3D graphic programming. You can use 4-d vectors, quaternions, 4x4 matrices, planes, etc.

mathutil_ref.c is the source file that defines all of the math functions as the reference.

The point of this project is: it can detect the instruction set of the CPU, and use SIMD instructions to accelerate some of the functions.

The CPU detection and function pointer replacement is default turned on for MSVC compilers and off for GCC, because I cannot use the AVX intrinsics without using -mavx option in GCC, and if I use it, GCC will optimize all of the ref implements to use AVX instructions, which destroys the dynamic code selection concept of this library that wants to be compatible with variant CPUs, but it optimizes better than MSVC. But this also allow GCC to optimize the code more complete, as we can produce different binaries for each architecture of the CPUs, and use another way to choose them.

The acceleration design is based on x86 system, if you want to use it on ARM, the SIMD instruction acceleration is unavailable.

Tested on Windows 7, compiled with MinGW-w64, no error or warning reported.

