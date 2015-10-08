#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

void crypto_aead_encrypt(
        unsigned char *c, size_t *clen,
        const unsigned char *h, size_t hlen,
        const unsigned char *p, size_t plen,
        const unsigned char *nonce,
        const unsigned char *key);

int crypto_aead_decrypt(
        unsigned char *p, size_t *plen,
        const unsigned char *h, size_t hlen,
        const unsigned char *c, size_t clen,
        const unsigned char *nonce,
        const unsigned char *key);


static int bench_cmp( const void *x, const void *y )
{
  const int64_t *ix = ( const int64_t * )x;
  const int64_t *iy = ( const int64_t * )y;
  return *ix - *iy;
}

#if   defined(__i386__)
static unsigned long long cpucycles( void )
{
  unsigned long long result;
  __asm__ __volatile__
  (
    ".byte 15;.byte 49"
    : "=A" ( result )
  );
  return result;
}
#elif defined(__x86_64__)
static unsigned long long cpucycles( void )
{
  unsigned long long result;
  __asm__ __volatile__
  (
    ".byte 15;.byte 49\n"
    "shlq $32,%%rdx\n"
    "orq %%rdx,%%rax"
    : "=a" ( result ) ::  "%rdx"
  );
  return result;
}
#elif defined(__arm__)
#include <unistd.h>
#include <sys/types.h>
#include <sys/syscall.h>
#include <linux/perf_event.h>
#include <errno.h>

static int fddev = -1;

__attribute__((constructor)) static void
init(void)
{
  static struct perf_event_attr attr;
  attr.type = PERF_TYPE_HARDWARE;
  attr.config = PERF_COUNT_HW_CPU_CYCLES;
  fddev = syscall(__NR_perf_event_open, &attr, 0, -1, -1, 0);
}

__attribute__((destructor)) static void
fini(void)
{
  close(fddev);
}

static unsigned long long cpucycles(void)
{
  unsigned long long result = 0;
  if (read(fddev, &result, sizeof(result)) < sizeof(result))
  {
    printf("error\n");
    return 0;
  }
  return result;
}
#elif TARGET_OS_IPHONE
#include <mach/mach_time.h>
#if 1
#define GHZ 1.4 /* iPad Air */
#endif
#if 0
#define GHZ 1.3 /* iPhone 5s & iPad Mini Retina */
#endif
static mach_timebase_info_data_t sTimebaseInfo;
static unsigned long long cpucycles( void )
{
    if ( sTimebaseInfo.denom == 0 ) {
        (void) mach_timebase_info(&sTimebaseInfo);
    }
    return ( ( ( double )sTimebaseInfo.numer / ( double )sTimebaseInfo.denom ) * GHZ ) * mach_absolute_time();
}
#else
#error "Don't know how to count cycles!"
#endif

void frequency()
{
	uint64_t t;
	printf("Estimating cycle counter frequency...");
	t = cpucycles();
	sleep(1);
	t = cpucycles() - t;
	printf("%f GHz\n", t/1e9);
}

void bench()
{
#define BENCH_SIZE     4096
#define BENCH_TRIALS     32
#define BENCH_MAXLEN   1536
  static unsigned char  in[BENCH_SIZE];
  static unsigned char out[BENCH_SIZE+32];
  static unsigned char  ad[BENCH_SIZE];
  static unsigned char   n[32];
  static unsigned char   k[32];
  static size_t outlen;
  static size_t adlen = 0;
  static unsigned long long median[BENCH_SIZE + 1];
  int i, j;

  printf( "#bytes  median  per byte\n" );

  /* 1 ... BENCH_MAXLEN */
  for( j = 0; j <= BENCH_SIZE; ++j )
  {
    uint64_t cycles[BENCH_TRIALS + 1];

    for( i = 0; i <= BENCH_TRIALS; ++i )
    {
      cycles[i] = cpucycles();
      crypto_aead_encrypt(out, &outlen, ad, adlen, in, j, n, k);
    }

    for( i = 0; i < BENCH_TRIALS; ++i )
      cycles[i] = cycles[i + 1] - cycles[i];

    qsort( cycles, BENCH_TRIALS, sizeof( uint64_t ), bench_cmp );
    median[j] = cycles[BENCH_TRIALS / 2];
  }

  for( j = 0; j <= BENCH_SIZE; j += 8 )
    printf( "%5d, %7.2f\n", j, ( double )median[j] / j );

  printf( "#%u   %6llu   %7.2f\n", BENCH_SIZE / 2, median[BENCH_SIZE/2], ( double )median[BENCH_SIZE/2] / (BENCH_SIZE / 2) );
  printf( "#%u   %6llu   %7.2f\n", BENCH_SIZE, median[BENCH_SIZE], ( double )median[BENCH_SIZE] / BENCH_SIZE );
  printf( "#long     long   %7.2f\n", ( double )( median[BENCH_SIZE] - median[BENCH_SIZE/2] ) / (BENCH_SIZE / 2.0) );
}

int main(int argc, char **argv)
{
  if( argc > 1 && 0 == strcmp(argv[1], "-f") )
  	frequency();
  bench();
  return 0;
}

