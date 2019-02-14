# mathutil
A utility for 3D graphic programming. You can use 4-d vectors, quaternions, 4x4 matrices, planes, etc.

mathutil_ref.c is the source file that defines all of the math functions as the reference.

The point of this project is: it can detect the instruction set of the CPU, and use SIMD instructions to accelerate some of the functions.

The acceleration design is based on x86 system, if you want to use it on ARM, the SIMD instruction acceleration is unavailable.

Tested on Windows 7, compiled with MinGW-w64.

