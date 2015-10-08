CFLAGS=-O3 -mcpu=cortex-a8 -mfpu=neon -std=c99 -Wall -Wextra -pedantic -Wno-long-long

all: check bench

bench:
	@$(CC) $(CFLAGS) -o bench ../../utils/bench.c mrs.c
	@./bench
	@rm bench

check:
	@$(CC) $(CFLAGS) -I../ -o check ../../utils/check.c mrs.c
	@./check
	@rm check
