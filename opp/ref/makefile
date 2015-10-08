CFLAGS=-O3 -march=native -std=c89 -Wall -Wextra -pedantic -Wno-long-long

all: check bench

debug:
	@$(CC) $(CFLAGS) -DOPP_DEBUG -o debug ../../utils/debug.c ../../utils/io.c opp.c
	@./debug
	@rm debug

bench:
	@$(CC) $(CFLAGS) -o bench ../../utils/bench.c opp.c
	@./bench
	@rm bench

check:
	@$(CC) $(CFLAGS) -I../ -o check ../../utils/check.c opp.c
	@./check
	@rm check

genkat:
	@$(CC) $(CFLAGS) -o genkat ../../utils/genkat.c opp.c
	@./genkat > ../kat.h
	@rm genkat
