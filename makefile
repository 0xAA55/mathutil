CC=gcc
LD=gcc
AR=gcc-ar
RANLIB=gcc-ranlib
CFLAGS=-Wall -O3 -flto -ffat-lto-objects -mavx -I..

OBJS=cpudetect.o
OBJS+=mathutil.o
OBJS+=mathutil_ref.o
OBJS+=mathutil_sse.o
OBJS+=mathutil_sse2.o
OBJS+=mathutil_sse3.o
OBJS+=mathutil_ssse3.o
OBJS+=mathutil_sse41.o

all: libmathutil.a
	
libmathutil.a: $(OBJS)
	$(AR) rcu $@ $+
	$(RANLIB) $@

clean:
	rm *.o *.a
