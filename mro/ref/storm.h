/*
   STORM reference source code package - reference C implementations

   Written in 2015 by Philipp Jovanovic <philipp@jovanovic.io>

   To the extent possible under law, the author(s) have dedicated all copyright
   and related and neighboring rights to this software to the public domain
   worldwide. This software is distributed without any warranty.

   You should have received a copy of the CC0 Public Domain Dedication along with
   this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/
#ifndef STORM_H
#define STORM_H

#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#define STORM_W 64             /* word size */
#define STORM_L 4              /* round number */
#define STORM_T (STORM_W *  4) /* tag size */
#define STORM_N (STORM_W *  2) /* nonce size */
#define STORM_K (STORM_W *  4) /* key size */
#define STORM_B (STORM_W * 16) /* permutation width */

typedef uint64_t storm_word_t;

typedef struct state__
{
    storm_word_t S[16];
} storm_state_t[1];

typedef enum tag__
{
    ABS_AD  = 0x00,
    ABS_MSG = 0x01
} tag_t;

/* workaround for C89 compilers */
#if !defined(__cplusplus) && (!defined(__STDC_VERSION__) || __STDC_VERSION__ < 199901L)
  #if   defined(_MSC_VER)
    #define STORM_INLINE __inline
  #elif defined(__GNUC__)
    #define STORM_INLINE __inline__
  #else
    #define STORM_INLINE
  #endif
#else
  #define STORM_INLINE inline
#endif

#define BITS(x) (sizeof(x) * CHAR_BIT)
#define BYTES(x) (((x) + 7) / 8)
#define WORDS(x) (((x) + (STORM_W-1)) / STORM_W)

static STORM_INLINE uint64_t load64(const void * in)
{
#if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
    uint64_t v;
    memcpy(&v, in, sizeof v);
    return v;
#else
    const uint8_t * p = (const uint8_t *)in;
    return ((uint64_t)p[0] <<  0) |
           ((uint64_t)p[1] <<  8) |
           ((uint64_t)p[2] << 16) |
           ((uint64_t)p[3] << 24) |
           ((uint64_t)p[4] << 32) |
           ((uint64_t)p[5] << 40) |
           ((uint64_t)p[6] << 48) |
           ((uint64_t)p[7] << 56);
#endif
}

static STORM_INLINE void store64(void * out, const uint64_t v)
{
#if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
    memcpy(out, &v, sizeof v);
#else
    uint8_t * p = (uint8_t *)out;
    p[0] = (uint8_t)(v >>  0);
    p[1] = (uint8_t)(v >>  8);
    p[2] = (uint8_t)(v >> 16);
    p[3] = (uint8_t)(v >> 24);
    p[4] = (uint8_t)(v >> 32);
    p[5] = (uint8_t)(v >> 40);
    p[6] = (uint8_t)(v >> 48);
    p[7] = (uint8_t)(v >> 56);
#endif
}

#define LOAD load64
#define STORE store64

static void* (* const volatile burn)(void*, int, size_t) = memset;

#if defined(STORM_DEBUG)

#include <inttypes.h>
#include <stdio.h>

#define FMT "016" PRIX64

static void print_state(storm_state_t state)
{
    static const char fmt[] = "%" FMT " "
                              "%" FMT " "
                              "%" FMT " "
                              "%" FMT "\n";
    const storm_word_t * S = state->S;
    printf(fmt, S[ 0],S[ 1],S[ 2],S[ 3]);
    printf(fmt, S[ 4],S[ 5],S[ 6],S[ 7]);
    printf(fmt, S[ 8],S[ 9],S[10],S[11]);
    printf(fmt, S[12],S[13],S[14],S[15]);
    printf("\n");
}

static void print_bytes(const uint8_t * in, size_t inlen)
{
    size_t i;
    for (i = 0; i < inlen; ++i)
    {
        printf("%02X ", in[i]);
        if (i%16 == 15)
        {
            printf("\n");
        }
    }
    printf("\n");
}
#endif

/* high-level operations */
void storm_aead_encrypt(
        unsigned char *c, size_t *clen,
        const unsigned char *h, size_t hlen,
        const unsigned char *m, size_t mlen,
        const unsigned char *nonce,
        const unsigned char *key);

int storm_aead_decrypt(
        unsigned char *m, size_t *mlen,
        const unsigned char *h, size_t hlen,
        const unsigned char *c, size_t clen,
        const unsigned char *nonce,
        const unsigned char *key);

/*foo*/

#endif
