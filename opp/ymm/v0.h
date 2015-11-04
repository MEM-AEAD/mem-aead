/*
    OPP - MEM AEAD source code package

    :copyright: (c) 2015 by Philipp Jovanovic and Samuel Neves
    :license: Creative Commons CC0 1.0
*/
#ifndef OPP_YMM_V0_H
#define OPP_YMM_V0_H

#include <stdint.h>
#include <string.h>
#include <immintrin.h>

#define BYTES(X) (((X) + 7) / 8)
#define WORDS(X) (((X) + (OPP_W - 1)) / OPP_W)

#define ADD256(A, B) _mm256_add_epi64((A), (B))
#define SUB256(A, B) _mm256_sub_epi64((A), (B))
#define XOR256(A, B) _mm256_xor_si256((A), (B))
#define  OR256(A, B)  _mm256_or_si256((A), (B))
#define AND256(A, B) _mm256_and_si256((A), (B))
#define SHL256(A, B) _mm256_slli_epi64((A), (B))
#define SHR256(A, B) _mm256_srli_epi64((A), (B))

#define ADD128(A, B) _mm_add_epi64((A), (B))
#define SUB128(A, B) _mm_sub_epi64((A), (B))
#define XOR128(A, B) _mm_xor_si128((A), (B))
#define  OR128(A, B)  _mm_or_si128((A), (B))
#define AND128(A, B) _mm_and_si128((A), (B))
#define SHL128(A, B) _mm_slli_epi64((A), (B))
#define SHR128(A, B) _mm_srli_epi64((A), (B))

#define ROT64(X, C) ( ((X) >> (C)) | ((X) << (64 - (C))) )

#define R16_128 _mm_setr_epi8( 2, 3, 4, 5, 6, 7, 0, 1, 10, 11, 12, 13, 14, 15, 8, 9 )
#define R24_128 _mm_setr_epi8( 3, 4, 5, 6, 7, 0, 1, 2, 11, 12, 13, 14, 15, 8, 9, 10 )

#define R16_256 _mm256_broadcastsi128_si256(R16_128)
#define R24_256 _mm256_broadcastsi128_si256(R24_128)

#define ROT128(X, C) (                                            \
      ((C) == 32) ? _mm_shuffle_epi32(X, _MM_SHUFFLE(2,3,0,1))    \
    : ((C) == 24) ? _mm_shuffle_epi8(X, R24_128)                  \
    : ((C) == 16) ? _mm_shuffle_epi8(X, R16_128)                  \
    : ((C) == 63) ? OR128(SHR128(X, C), ADD128(X, X))             \
    :               OR128(SHR128(X, C), SHL128(X, 64-(C)))        \
)

#define ROT256(X, C) (                                               \
      ((C) == 32) ? _mm256_shuffle_epi32(X, _MM_SHUFFLE(2,3,0,1))    \
    : ((C) == 24) ? _mm256_shuffle_epi8(X, R24_256)                  \
    : ((C) == 16) ? _mm256_shuffle_epi8(X, R16_256)                  \
    : ((C) == 63) ? OR256(SHR256(X, C), ADD256(X, X))                \
    :               OR256(SHR256(X, C), SHL256(X, 64-(C)))           \
)

#define LOADU128(X) _mm_loadu_si128((const __m128i *)((X)))
#define STOREU128(X, V) _mm_storeu_si128((__m128i *)((X)), (V))

#define LOADU256(X) _mm256_loadu_si256((const __m256i *)((X)))
#define STOREU256(X, V) _mm256_storeu_si256((__m256i *)((X)), (V))

#endif
