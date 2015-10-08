CFLAGS=-O3 -march=core-avx2 -mtune=haswell -std=c89 -Wall -Wextra -pedantic -Wno-long-long

all: check bench

bench: mro.c
	$(CC) $(CFLAGS) -o bench ../../utils/bench.c mro.c
	./bench
	rm bench

check: mro.c
	$(CC) $(CFLAGS) -I../ -o check ../../utils/check.c mro.c
	./check
	rm check
