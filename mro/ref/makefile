CFLAGS=-O3 -std=c89 -Wall -Wextra -pedantic -Wno-long-long

all: check bench

debug:
	@$(CC) $(CFLAGS) -DMRO_DEBUG -o debug ../../utils/debug.c ../../utils/io.c mro.c
	@./debug
	@rm debug

bench:
	@$(CC) $(CFLAGS) -o bench ../../utils/bench.c mro.c
	@./bench
	@rm bench

check:
	@$(CC) $(CFLAGS) -I../ -o check ../../utils/check.c mro.c
	@./check
	@rm check

genkat:
	@$(CC) $(CFLAGS) -o genkat ../../utils/genkat.c mro.c
	@./genkat > ../kat.h
	@rm genkat
