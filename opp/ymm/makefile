CFLAGS=-O3 -march=core-avx2 -mtune=haswell -std=c89 -Wall -Wextra -pedantic -Wno-long-long

all: check bench

bench: opp.c
	$(CC) $(CFLAGS) -o bench ../../utils/bench.c opp.c ../api.c
	./bench
	rm bench

check: opp.c
	$(CC) $(CFLAGS) -I../ -o check ../../utils/check.c opp.c ../api.c
	./check
	rm check
