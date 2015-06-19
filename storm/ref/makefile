CFLAGS=-O3 -march=native -std=c89 -Wall -Wextra -pedantic -Wno-long-long

all: check bench

bench:
	@$(CC) $(CFLAGS) -o bench ../../utils/bench.c storm.c
	@./bench
	@rm bench

check:
	@$(CC) $(CFLAGS) -I. -o check ../../utils/check.c storm.c
	@./check
	@rm check

genkat:
	@$(CC) $(CFLAGS) -o genkat ../../utils/genkat.c storm.c
	@./genkat > kat.h
	@rm genkat
