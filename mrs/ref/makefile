CFLAGS=-O3 -march=native -std=c89 -Wall -Wextra -pedantic -Wno-long-long

all: check bench

debug:
	@$(CC) $(CFLAGS) -DMRS_DEBUG -o debug ../../utils/debug.c mrs.c
	@./debug
	@rm debug

bench:
	@$(CC) $(CFLAGS) -o bench ../../utils/bench.c mrs.c
	@./bench
	@rm bench

check:
	@$(CC) $(CFLAGS) -I../ -o check ../../utils/check.c mrs.c
	@./check
	@rm check

genkat:
	@$(CC) $(CFLAGS) -o genkat ../../utils/genkat.c mrs.c
	@./genkat > ../kat.h
	@rm genkat
