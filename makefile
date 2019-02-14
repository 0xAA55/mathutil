CC=gcc
LD=gcc
AR = $(GCC_PREFIX)gcc-ar
RANLIB = $(GCC_PREFIX)gcc-ranlib
OPTIMIZATIONS=-O3 -fdata-sections -ffunction-sections -fmerge-all-constants -flto -fuse-linker-plugin -ffat-lto-objects
CFLAGS=-Wall $(OPTIMIZATIONS) -mavx

libmathutil.a: cpudetect.o mathutil.o mathutil_ref.o mathutil_sse.o mathutil_sse2.o mathutil_sse3.o mathutil_ssse3.o mathutil_sse41.o
	$(AR) rcu $@ $+
	$(RANLIB) $@

clean:
	del *.o *.a
