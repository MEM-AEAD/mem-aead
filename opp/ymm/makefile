CFLAGS=-O3 -march=core-avx2 -mtune=haswell -std=c89 -Wall -Wextra -pedantic -Wno-long-long

all: check bench

debug:
	@$(CC) $(CFLAGS) -DOPP_DEBUG -o debug ../../utils/debug.c ../../utils/io.c opp.c
	@./debug
	@rm debug

bench: opp.c
	@$(CC) $(CFLAGS) -o bench ../../utils/bench.c opp.c
	@./bench
	@rm bench

check: opp.c
	@$(CC) $(CFLAGS) -I../ -o check ../../utils/check.c opp.c
	@./check
	@rm check
