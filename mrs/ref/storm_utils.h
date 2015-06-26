/*
   STORM reference source code package - reference C implementations

   Written in 2015 by Philipp Jovanovic <philipp@jovanovic.io>

   To the extent possible under law, the author(s) have dedicated all copyright
   and related and neighboring rights to this software to the public domain
   worldwide. This software is distributed without any warranty.

   You should have received a copy of the CC0 Public Domain Dedication along with
   this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/
#ifndef UTILS_H
#define UTILS_H

#include <limits.h>
#include <stddef.h>
#include <string.h>
#include <stdint.h>

/* Workaround for C89 compilers */
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

static void* (* const volatile burn)(void*, int, size_t) = memset;

#endif
