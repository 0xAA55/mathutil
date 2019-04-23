CC=gcc
LD=gcc
AR = $(GCC_PREFIX)gcc-ar
RANLIB = $(GCC_PREFIX)gcc-ranlib
OPTIMIZATIONS=-g -O3 -fdata-sections -ffunction-sections -fmerge-all-constants -flto -fuse-linker-plugin -ffat-lto-objects
CFLAGS=-Wall $(OPTIMIZATIONS) -mavx -I..

OBJS=cpudetect.o
OBJS+=mathutil.o
OBJS+=mathutil_ref.o
OBJS+=mathutil_sse.o
OBJS+=mathutil_sse2.o
OBJS+=mathutil_sse3.o
OBJS+=mathutil_ssse3.o
OBJS+=mathutil_sse41.o

all: libmathutil.a

-include $(OBJS:.o=.d)

%.o: %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o
	$(CC) -MM $(CFLAGS) $*.c > $*.d
	
libmathutil.a: $(OBJS)
	$(AR) rcu $@ $+
	$(RANLIB) $@

clean:
	del *.o *.a *.d
