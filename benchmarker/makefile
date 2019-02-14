CC=gcc
LD=gcc
OPTIMIZATIONS=-O3 -fdata-sections -ffunction-sections -fmerge-all-constants -flto -fuse-linker-plugin -Wl,--gc-sections
CFLAGS=-Wall $(OPTIMIZATIONS) -I../
LDLIBS=-lmathutil
LDFLAGS=$(OPTIMIZATIONS)

mathutil_bench.exe: rttimer.o mathutil_bench.o
	$(LD) $(LDFLAGS) $^ $(LDLIBS) -o $@

clean:
	del *.o *.exe
