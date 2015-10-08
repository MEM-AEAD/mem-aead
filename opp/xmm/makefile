CFLAGS=-O3 -march=native -std=c89 -Wall -Wextra -pedantic -Wno-long-long

all: check bench

debug:
	@$(CC) $(CFLAGS) -o debug ../../utils/debug.c opp.c api.c
	@./debug
	@rm debug

bench:
	@$(CC) $(CFLAGS) -o bench ../../utils/bench.c opp.c api.c
	@./bench
	@rm bench

check:
	@$(CC) $(CFLAGS) -I../ -o check ../../utils/check.c opp.c api.c
	@./check
	@rm check
