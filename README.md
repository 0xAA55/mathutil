# mathutil
## Introduction
### What is mathutil?
*mathutil* is a cross-platform library for **3D game developing** written in ***C***. It providing functions for **3D、4D vector, plane, quarternion, and 4x4 matrix** related calculations, and some basic math functions like **sin**, **cos**, etc.

### Why use mathutil?
*mathutil* is a replacement for something like *d3dx9.h*, but it also provides **source code**. If you want to make your project **portable**, then don't use d3dx9. Use *mathutil* instead.

### Optimization
*mathutil* uses SIMD instructions for optimization if available. If *mathutil* is compiled by MSVC, it detects **CPUID** and replaces functions with one that uses SIMD instructions and compatible with the current CPU.

Currently, for MSVC compilers, *mathutil* only provides **SSE optimized functions** on x86-64 platform.

But if *mathutil* is compiled by GCC, it will never try to detect CPUID because GCC is better on optimization than MSVC, and GCC doesn't allow using intrinsics which doesn't match the current architecture.

## Usage

### Compile *mathutil*
The makefile will help you to create a static library with GCC optimization flags set. To build a library for MSVC, open the project file located on 'msvc' directory to run build.
***NOTE*** that the optimization flags for GCC included -mavx. You may need to change this for other architecture CPUs.

### Function references
See [**manual**](doc/mathutil_reference_manual.md) more informations.
