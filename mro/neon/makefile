CFLAGS=-O3 -mcpu=cortex-a8 -mfpu=neon -std=c99 -Wall -Wextra -pedantic -Wno-long-long

all: check bench

debug:
	@$(CC) $(CFLAGS) -o debug ../../utils/debug.c mro.c
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

