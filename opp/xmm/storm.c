/*
   STORM reference source code package - reference C implementations

   Written in 2015 by Philipp Jovanovic <philipp@jovanovic.io>

   To the extent possible under law, the author(s) have dedicated all copyright
   and related and neighboring rights to this software to the public domain
   worldwide. This software is distributed without any warranty.

   You should have received a copy of the CC0 Public Domain Dedication along with
   this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/
#include <string.h>
#include "storm.h"

#if defined(_MSC_VER)
    #include <intrin.h>
#else
    #include <x86intrin.h>
#endif

#define STORM_N (STORM_W *  2)   /* nonce size */
#define STORM_K (STORM_W *  4)   /* key size */
#define STORM_B (STORM_W * 16)   /* permutation width */

#define BYTES(x) (((x) + 7) / 8)
#define WORDS(x) (((x) + (STORM_W-1)) / STORM_W)

#if defined(_MSC_VER)
    #define ALIGN(x) __declspec(align(x))
#else
    #define ALIGN(x) __attribute__((aligned(x)))
#endif

#define XOR(A, B) _mm_xor_si128((A), (B))
#define AND(A, B) _mm_and_si128((A), (B))
#define ADD(A, B) _mm_add_epi64((A), (B))
#define SUB(A, B) _mm_sub_epi64((A), (B))

/* constants for ROTR */
#define R0 32
#define R1 24
#define R2 16
#define R3 63

/* constants for ROTL */
#define L0 (STORM_W - R0)
#define L1 (STORM_W - R1)
#define L2 (STORM_W - R2)
#define L3 (STORM_W - R3)

#define LOAD(in) _mm_load_si128((__m128i*)(in))
#define STORE(out, x) _mm_store_si128((__m128i*)(out), (x))
#define LOADU(in) _mm_loadu_si128((__m128i*)(in))
#define STOREU(out, x) _mm_storeu_si128((__m128i*)(out), (x))
#define LOADL(in) _mm_loadl_epi64((__m128i*)(in))

#  if defined(__XOP__)
#define ROT(X, C) _mm_roti_epi64((X), -(C))
#elif defined(__SSSE3__)
#define ROT(X, C)                                                                                           \
(                                                                                                           \
        (C) == 32 ? _mm_shuffle_epi8((X), _mm_set_epi8(11,10, 9, 8, 15,14,13,12,  3, 2, 1, 0,  7, 6, 5, 4)) \
    :   (C) == 24 ? _mm_shuffle_epi8((X), _mm_set_epi8(10, 9, 8,15, 14,13,12,11,  2, 1, 0, 7,  6, 5, 4, 3)) \
    :   (C) == 16 ? _mm_shuffle_epi8((X), _mm_set_epi8( 9, 8,15,14, 13,12,11,10,  1, 0, 7, 6,  5, 4, 3, 2)) \
    :   (C) == 63 ? _mm_or_si128(_mm_add_epi64((X), (X)), _mm_srli_epi64((X), 63))                          \
    :   /* else */  _mm_or_si128(_mm_srli_epi64((X), (C)), _mm_slli_epi64((X), 64 - (C)))                   \
)
#else
#define ROT(X, C)                                                                         \
(                                                                                         \
        (C) == 63 ? _mm_or_si128(_mm_add_epi64((X), (X)), _mm_srli_epi64((X), 63))        \
    :   /* else */  _mm_or_si128(_mm_srli_epi64((X), (C)), _mm_slli_epi64((X), 64 - (C))) \
)
#endif

/* quarter round */
#define G(S)                                           \
do                                                     \
{                                                      \
    S[0] = ADD(S[0], S[2]);    S[1] = ADD(S[1], S[3]); \
    S[6] = XOR(S[6], S[0]);    S[7] = XOR(S[7], S[1]); \
    S[6] = ROT(S[6],   R0);    S[7] = ROT(S[7],   R0); \
                                                       \
    S[4] = ADD(S[4], S[6]);    S[5] = ADD(S[5], S[7]); \
    S[2] = XOR(S[2], S[4]);    S[3] = XOR(S[3], S[5]); \
    S[2] = ROT(S[2],   R1);    S[3] = ROT(S[3],   R1); \
                                                       \
    S[0] = ADD(S[0], S[2]);    S[1] = ADD(S[1], S[3]); \
    S[6] = XOR(S[6], S[0]);    S[7] = XOR(S[7], S[1]); \
    S[6] = ROT(S[6],   R2);    S[7] = ROT(S[7],   R2); \
                                                       \
    S[4] = ADD(S[4], S[6]);    S[5] = ADD(S[5], S[7]); \
    S[2] = XOR(S[2], S[4]);    S[3] = XOR(S[3], S[5]); \
    S[2] = ROT(S[2],   R3);    S[3] = ROT(S[3],   R3); \
} while(0)

/* inverse quarter round */
#define GI(S)                                          \
do                                                     \
{                                                      \
    S[2] = ROT(S[2],  L3);     S[3] = ROT(S[3],   L3); \
    S[2] = XOR(S[2], S[4]);    S[3] = XOR(S[3], S[5]); \
    S[4] = SUB(S[4], S[6]);    S[5] = SUB(S[5], S[7]); \
                                                       \
    S[6] = ROT(S[6],   L2);    S[7] = ROT(S[7],   L2); \
    S[6] = XOR(S[6], S[0]);    S[7] = XOR(S[7], S[1]); \
    S[0] = SUB(S[0], S[2]);    S[1] = SUB(S[1], S[3]); \
                                                       \
    S[2] = ROT(S[2],   L1);    S[3] = ROT(S[3],   L1); \
    S[2] = XOR(S[2], S[4]);    S[3] = XOR(S[3], S[5]); \
    S[4] = SUB(S[4], S[6]);    S[5] = SUB(S[5], S[7]); \
                                                       \
    S[6] = ROT(S[6],   L0);    S[7] = ROT(S[7],   L0); \
    S[6] = XOR(S[6], S[0]);    S[7] = XOR(S[7], S[1]); \
    S[0] = SUB(S[0], S[2]);    S[1] = SUB(S[1], S[3]); \
} while(0)

#if defined(__SSSE3__)
#define DIAGONALIZE(S)                     \
do                                         \
{                                          \
    __m128i T[2];                          \
                                           \
    T[0] = _mm_alignr_epi8(S[3], S[2], 8); \
    T[1] = _mm_alignr_epi8(S[2], S[3], 8); \
    S[2] = T[0];                           \
    S[3] = T[1];                           \
                                           \
    T[0] = S[4];                           \
    S[4] = S[5];                           \
    S[5] = T[0];                           \
                                           \
    T[0] = _mm_alignr_epi8(S[7], S[6], 8); \
    T[1] = _mm_alignr_epi8(S[6], S[7], 8); \
    S[6] = T[1];                           \
    S[7] = T[0];                           \
} while(0)

#define UNDIAGONALIZE(S)                   \
do                                         \
{                                          \
    __m128i T[2];                          \
                                           \
    T[0] = _mm_alignr_epi8(S[2], S[3], 8); \
    T[1] = _mm_alignr_epi8(S[3], S[2], 8); \
    S[2] = T[0];                           \
    S[3] = T[1];                           \
                                           \
    T[0] = S[4];                           \
    S[4] = S[5];                           \
    S[5] = T[0];                           \
                                           \
    T[0] = _mm_alignr_epi8(S[6], S[7], 8); \
    T[1] = _mm_alignr_epi8(S[7], S[6], 8); \
    S[6] = T[1];                           \
    S[7] = T[0];                           \
} while(0)

#else

#define DIAGONALIZE(S)                                               \
do                                                                   \
{                                                                    \
    __m128i T[2];                                                    \
                                                                     \
    T[0] = S[6]; T[1] = S[2];                                        \
    S[6] = S[4]; S[4] = S[5]; S[5] = S[6];                           \
    S[6] = _mm_unpackhi_epi64(S[7], _mm_unpacklo_epi64(T[0], T[0])); \
    S[7] = _mm_unpackhi_epi64(T[0], _mm_unpacklo_epi64(S[7], S[7])); \
    S[2] = _mm_unpackhi_epi64(S[2], _mm_unpacklo_epi64(S[3], S[3])); \
    S[3] = _mm_unpackhi_epi64(S[3], _mm_unpacklo_epi64(T[1], T[1])); \
} while(0)

#define UNDIAGONALIZE(S)                                             \
do                                                                   \
{                                                                    \
    __m128i T[2];                                                    \
                                                                     \
    T[0] = S[4]; S[4] = S[5]; S[5] = T[0];                           \
    T[0] = S[2]; T[1] = S[6];                                        \
    S[2] = _mm_unpackhi_epi64(S[3], _mm_unpacklo_epi64(S[2], S[2])); \
    S[3] = _mm_unpackhi_epi64(T[0], _mm_unpacklo_epi64(S[3], S[3])); \
    S[6] = _mm_unpackhi_epi64(S[6], _mm_unpacklo_epi64(S[7], S[7])); \
    S[7] = _mm_unpackhi_epi64(S[7], _mm_unpacklo_epi64(T[1], T[1])); \
} while(0)

#endif

#define F(S)          \
do                    \
{                     \
    G(S);             \
    DIAGONALIZE(S);   \
    G(S);             \
    UNDIAGONALIZE(S); \
} while(0)

#define FI(S)         \
do                    \
{                     \
    DIAGONALIZE(S);   \
    GI(S);            \
    UNDIAGONALIZE(S); \
    GI(S);            \
} while(0)

#define PERMUTE(S)               \
do                               \
{                                \
    size_t i;                    \
    for(i = 0; i < STORM_R; ++i) \
    {                            \
        F(S);                    \
    }                            \
} while(0)

#define PERMUTE_INVERSE(S)       \
do                               \
{                                \
    size_t i;                    \
    for(i = 0; i < STORM_R; ++i) \
    {                            \
        FI(S);                   \
    }                            \
} while(0)

#define PAD(BLOCK, BLOCKLEN, IN, INLEN) \
do                                      \
{                                       \
    memset(BLOCK, 0, BLOCKLEN);         \
    memcpy(BLOCK, IN, INLEN);           \
    BLOCK[INLEN] = 0x01;                \
    BLOCK[BLOCKLEN - 1] |= 0x80;        \
} while(0)

#define LOAD_KEY(K, KEY, IV, IVLEN)                \
do                                                 \
{                                                  \
    size_t i;                                      \
    K[0] = _mm_set_epi64x(0,0);                    \
    K[1] = _mm_set_epi64x(0,0);                    \
    for (i = 0; i < IVLEN; ++i)                    \
    {                                              \
        K[i] = LOADU(IV + i * 2 * BYTES(STORM_W)); \
    }                                              \
    K[2] = LOADU(KEY +  0);                        \
    K[3] = LOADU(KEY + 16);                        \
    K[4] = _mm_set_epi64x(0,0);                    \
    K[5] = _mm_set_epi64x(0,0);                    \
    K[6] = _mm_set_epi64x(STORM_R, STORM_W);       \
    K[7] = _mm_set_epi64x(0, STORM_T);             \
    PERMUTE(K);                                    \
} while(0)

#if defined(M4)
#define UPDATE_KEY(K)                                                        \
do                                                                           \
{                                                                            \
    __m128i T00, T01, T10, T11, T20, T21, T30, T31;                          \
    T00 = _mm_xor_si128(_mm_slli_epi64(K[0], 53), _mm_srli_epi64(K[0], 11)); \
    T01 = _mm_xor_si128(_mm_slli_epi64(K[1], 53), _mm_srli_epi64(K[1], 11)); \
    T10 = _mm_alignr_epi8(K[3], K[2], 8);                                    \
    T11 = _mm_alignr_epi8(K[4], K[3], 8);                                    \
    T20 = _mm_slli_epi64(T10, 13);                                           \
    T21 = _mm_slli_epi64(T11, 13);                                           \
    T30 = _mm_xor_si128(T20, T00);                                           \
    T31 = _mm_xor_si128(T21, T01);                                           \
    K[0] = K[2]; K[1] = K[3];                                                \
    K[2] = K[4]; K[3] = K[5];                                                \
    K[4] = K[6]; K[5] = K[7];                                                \
    K[6] = T30;  K[7] = T31;                                                 \
} while(0)
#else
#define UPDATE_KEY(K)                                                                                   \
do                                                                                                      \
{                                                                                                       \
    __m128i T = XOR(ROT(_mm_set_epi64x(0, K[0][0]), 9), _mm_srli_epi64(_mm_set_epi64x(0, K[4][1]), 7)); \
    K[0] = _mm_set_epi64x(K[1][0], K[0][1]);                                                            \
    K[1] = _mm_set_epi64x(K[2][0], K[1][1]);                                                            \
    K[2] = _mm_set_epi64x(K[3][0], K[2][1]);                                                            \
    K[3] = _mm_set_epi64x(K[4][0], K[3][1]);                                                            \
    K[4] = _mm_set_epi64x(K[5][0], K[4][1]);                                                            \
    K[5] = _mm_set_epi64x(K[6][0], K[5][1]);                                                            \
    K[6] = _mm_set_epi64x(K[7][0], K[6][1]);                                                            \
    K[7] = _mm_set_epi64x(T[0],    K[7][1]);                                                            \
} while(0)
#endif

#define ABSORB_BLOCK(S, X, IN)         \
do                                     \
{                                      \
    __m128i B[8];                      \
    UPDATE_KEY(X);                     \
    B[0] = XOR(X[0], LOADU(IN +   0)); \
    B[1] = XOR(X[1], LOADU(IN +  16)); \
    B[2] = XOR(X[2], LOADU(IN +  32)); \
    B[3] = XOR(X[3], LOADU(IN +  48)); \
    B[4] = XOR(X[4], LOADU(IN +  64)); \
    B[5] = XOR(X[5], LOADU(IN +  80)); \
    B[6] = XOR(X[6], LOADU(IN +  96)); \
    B[7] = XOR(X[7], LOADU(IN + 112)); \
    PERMUTE(B);                        \
    S[0] = XOR(S[0], B[0]);            \
    S[1] = XOR(S[1], B[1]);            \
    S[2] = XOR(S[2], B[2]);            \
    S[3] = XOR(S[3], B[3]);            \
    S[4] = XOR(S[4], B[4]);            \
    S[5] = XOR(S[5], B[5]);            \
    S[6] = XOR(S[6], B[6]);            \
    S[7] = XOR(S[7], B[7]);            \
} while(0)

#define ABSORB_LASTBLOCK(S, X, IN, INLEN)          \
do                                                 \
{                                                  \
    ALIGN(64) unsigned char BLOCK[BYTES(STORM_B)]; \
    PAD(BLOCK, sizeof BLOCK, IN, INLEN);           \
    ABSORB_BLOCK(S, X, BLOCK);                     \
} while(0)

#define ENCRYPT_BLOCK(S, X, OUT, IN)   \
do                                     \
{                                      \
    __m128i B[8];                      \
    UPDATE_KEY(X);                     \
    B[0] = XOR(X[0], LOADU(IN +   0)); \
    B[1] = XOR(X[1], LOADU(IN +  16)); \
    B[2] = XOR(X[2], LOADU(IN +  32)); \
    B[3] = XOR(X[3], LOADU(IN +  48)); \
    B[4] = XOR(X[4], LOADU(IN +  64)); \
    B[5] = XOR(X[5], LOADU(IN +  80)); \
    B[6] = XOR(X[6], LOADU(IN +  96)); \
    B[7] = XOR(X[7], LOADU(IN + 112)); \
    PERMUTE(B);                        \
    STORE(OUT +   0, XOR(B[0], X[0])); \
    STORE(OUT +  16, XOR(B[1], X[1])); \
    STORE(OUT +  32, XOR(B[2], X[2])); \
    STORE(OUT +  48, XOR(B[3], X[3])); \
    STORE(OUT +  64, XOR(B[4], X[4])); \
    STORE(OUT +  80, XOR(B[5], X[5])); \
    STORE(OUT +  96, XOR(B[6], X[6])); \
    STORE(OUT + 112, XOR(B[7], X[7])); \
    S[0] = XOR(S[0], LOADU(IN +   0)); \
    S[1] = XOR(S[1], LOADU(IN +  16)); \
    S[2] = XOR(S[2], LOADU(IN +  32)); \
    S[3] = XOR(S[3], LOADU(IN +  48)); \
    S[4] = XOR(S[4], LOADU(IN +  64)); \
    S[5] = XOR(S[5], LOADU(IN +  80)); \
    S[6] = XOR(S[6], LOADU(IN +  96)); \
    S[7] = XOR(S[7], LOADU(IN + 112)); \
} while(0)

#define ENCRYPT_LASTBLOCK(S, X, OUT, IN, INLEN)                       \
do                                                                    \
{                                                                     \
    if (INLEN > 0)                                                    \
    {                                                                 \
        __m128i B[8];                                                 \
        ALIGN(64) unsigned char BLOCK[BYTES(STORM_B)];                \
        UPDATE_KEY(X);                                                \
        B[0] = X[0];                                                  \
        B[1] = X[1];                                                  \
        B[2] = X[2];                                                  \
        B[3] = X[3];                                                  \
        B[4] = X[4];                                                  \
        B[5] = X[5];                                                  \
        B[6] = X[6];                                                  \
        B[7] = XOR(X[7], _mm_set_epi64x(INLEN, 0));                   \
        PERMUTE(B);                                                   \
        PAD(BLOCK, sizeof BLOCK, IN, INLEN);                          \
        STORE(BLOCK +   0, XOR(B[0], XOR(X[0], LOADU(BLOCK +   0)))); \
        STORE(BLOCK +  16, XOR(B[1], XOR(X[1], LOADU(BLOCK +  16)))); \
        STORE(BLOCK +  32, XOR(B[2], XOR(X[2], LOADU(BLOCK +  32)))); \
        STORE(BLOCK +  48, XOR(B[3], XOR(X[3], LOADU(BLOCK +  48)))); \
        STORE(BLOCK +  64, XOR(B[4], XOR(X[4], LOADU(BLOCK +  64)))); \
        STORE(BLOCK +  80, XOR(B[5], XOR(X[5], LOADU(BLOCK +  80)))); \
        STORE(BLOCK +  96, XOR(B[6], XOR(X[6], LOADU(BLOCK +  96)))); \
        STORE(BLOCK + 112, XOR(B[7], XOR(X[7], LOADU(BLOCK + 112)))); \
        memcpy(OUT, BLOCK, INLEN);                                    \
        PAD(BLOCK, sizeof BLOCK, IN, INLEN);                          \
        S[0] = XOR(S[0], LOADU(BLOCK +   0));                         \
        S[1] = XOR(S[1], LOADU(BLOCK +  16));                         \
        S[2] = XOR(S[2], LOADU(BLOCK +  32));                         \
        S[3] = XOR(S[3], LOADU(BLOCK +  48));                         \
        S[4] = XOR(S[4], LOADU(BLOCK +  64));                         \
        S[5] = XOR(S[5], LOADU(BLOCK +  80));                         \
        S[6] = XOR(S[6], LOADU(BLOCK +  96));                         \
        S[7] = XOR(S[7], LOADU(BLOCK + 112));                         \
    }                                                                 \
} while(0)

#define DECRYPT_BLOCK(S, X, OUT, IN)    \
do                                      \
{                                       \
    __m128i B[8];                       \
    UPDATE_KEY(X);                      \
    B[0] = XOR(X[0], LOADU(IN +   0));  \
    B[1] = XOR(X[1], LOADU(IN +  16));  \
    B[2] = XOR(X[2], LOADU(IN +  32));  \
    B[3] = XOR(X[3], LOADU(IN +  48));  \
    B[4] = XOR(X[4], LOADU(IN +  64));  \
    B[5] = XOR(X[5], LOADU(IN +  80));  \
    B[6] = XOR(X[6], LOADU(IN +  96));  \
    B[7] = XOR(X[7], LOADU(IN + 112));  \
    PERMUTE_INVERSE(B);                 \
    STORE(OUT +   0, XOR(B[0], X[0]));  \
    STORE(OUT +  16, XOR(B[1], X[1]));  \
    STORE(OUT +  32, XOR(B[2], X[2]));  \
    STORE(OUT +  48, XOR(B[3], X[3]));  \
    STORE(OUT +  64, XOR(B[4], X[4]));  \
    STORE(OUT +  80, XOR(B[5], X[5]));  \
    STORE(OUT +  96, XOR(B[6], X[6]));  \
    STORE(OUT + 112, XOR(B[7], X[7]));  \
    S[0] = XOR(S[0], LOADU(OUT +   0)); \
    S[1] = XOR(S[1], LOADU(OUT +  16)); \
    S[2] = XOR(S[2], LOADU(OUT +  32)); \
    S[3] = XOR(S[3], LOADU(OUT +  48)); \
    S[4] = XOR(S[4], LOADU(OUT +  64)); \
    S[5] = XOR(S[5], LOADU(OUT +  80)); \
    S[6] = XOR(S[6], LOADU(OUT +  96)); \
    S[7] = XOR(S[7], LOADU(OUT + 112)); \
} while(0)

#define DECRYPT_LASTBLOCK(S, X, OUT, IN, INLEN)                       \
do                                                                    \
{                                                                     \
    if (INLEN > 0)                                                    \
    {                                                                 \
        __m128i B[8];                                                 \
        ALIGN(64) unsigned char BLOCK[BYTES(STORM_B)];                \
        UPDATE_KEY(X);                                                \
        B[0] = X[0];                                                  \
        B[1] = X[1];                                                  \
        B[2] = X[2];                                                  \
        B[3] = X[3];                                                  \
        B[4] = X[4];                                                  \
        B[5] = X[5];                                                  \
        B[6] = X[6];                                                  \
        B[7] = XOR(X[7], _mm_set_epi64x(INLEN, 0));                   \
        PERMUTE(B);                                                   \
        PAD(BLOCK, sizeof BLOCK, IN, INLEN);                          \
        STORE(BLOCK +   0, XOR(B[0], XOR(X[0], LOADU(BLOCK +   0)))); \
        STORE(BLOCK +  16, XOR(B[1], XOR(X[1], LOADU(BLOCK +  16)))); \
        STORE(BLOCK +  32, XOR(B[2], XOR(X[2], LOADU(BLOCK +  32)))); \
        STORE(BLOCK +  48, XOR(B[3], XOR(X[3], LOADU(BLOCK +  48)))); \
        STORE(BLOCK +  64, XOR(B[4], XOR(X[4], LOADU(BLOCK +  64)))); \
        STORE(BLOCK +  80, XOR(B[5], XOR(X[5], LOADU(BLOCK +  80)))); \
        STORE(BLOCK +  96, XOR(B[6], XOR(X[6], LOADU(BLOCK +  96)))); \
        STORE(BLOCK + 112, XOR(B[7], XOR(X[7], LOADU(BLOCK + 112)))); \
        memcpy(OUT, BLOCK, INLEN);                                    \
        PAD(BLOCK, sizeof BLOCK, OUT, INLEN);                         \
        S[0] = XOR(S[0], LOADU(BLOCK +   0));                         \
        S[1] = XOR(S[1], LOADU(BLOCK +  16));                         \
        S[2] = XOR(S[2], LOADU(BLOCK +  32));                         \
        S[3] = XOR(S[3], LOADU(BLOCK +  48));                         \
        S[4] = XOR(S[4], LOADU(BLOCK +  64));                         \
        S[5] = XOR(S[5], LOADU(BLOCK +  80));                         \
        S[6] = XOR(S[6], LOADU(BLOCK +  96));                         \
        S[7] = XOR(S[7], LOADU(BLOCK + 112));                         \
    }                                                                 \
} while(0)

#define ABSORB_DATA(S, X, IN, INLEN)                        \
do                                                          \
{                                                           \
    if(INLEN > 0)                                           \
    {                                                       \
        size_t i = 0;                                       \
        size_t l = INLEN;                                   \
        while(l >= BYTES(STORM_B))                          \
        {                                                   \
            ABSORB_BLOCK(S, X, IN + i * BYTES(STORM_B));    \
            i += 1; l -= BYTES(STORM_B);                    \
        }                                                   \
        ABSORB_LASTBLOCK(S, X, IN + i * BYTES(STORM_B), l); \
    }                                                       \
} while(0)

#define ENCRYPT_DATA(S, X, OUT, IN, INLEN)                                             \
do                                                                                     \
{                                                                                      \
    if(INLEN > 0)                                                                      \
    {                                                                                  \
        size_t i = 0;                                                                  \
        size_t l = INLEN;                                                              \
        while(l >= BYTES(STORM_B))                                                     \
        {                                                                              \
            ENCRYPT_BLOCK(S, X, OUT + i * BYTES(STORM_B), IN + i * BYTES(STORM_B));    \
            i += 1; l -= BYTES(STORM_B);                                               \
        }                                                                              \
        ENCRYPT_LASTBLOCK(S, X, OUT + i * BYTES(STORM_B), IN + i * BYTES(STORM_B), l); \
    }                                                                                  \
} while(0)

#define DECRYPT_DATA(S, X, OUT, IN, INLEN)                                             \
do                                                                                     \
{                                                                                      \
    if(INLEN > 0)                                                                      \
    {                                                                                  \
        size_t i = 0;                                                                  \
        size_t l = INLEN;                                                              \
        while(l >= BYTES(STORM_B))                                                     \
        {                                                                              \
            DECRYPT_BLOCK(S, X, OUT + i * BYTES(STORM_B), IN + i * BYTES(STORM_B));    \
            i += 1; l -= BYTES(STORM_B);                                               \
        }                                                                              \
        DECRYPT_LASTBLOCK(S, X, OUT + i * BYTES(STORM_B), IN + i * BYTES(STORM_B), l); \
    }                                                                                  \
} while(0)

#define SQUEEZE_TAG(S, TAG)               \
do                                        \
{                                         \
    PERMUTE(S);                           \
    STOREU(TAG, S[0]);                    \
    STOREU(TAG + BYTES(STORM_T)/2, S[1]); \
} while(0)

static void* (* const volatile burn)(void*, int, size_t) = memset;

void storm_aead_encrypt(
    unsigned char *c, size_t *clen,
    const unsigned char *h, size_t hlen,
    const unsigned char *m, size_t mlen,
    const unsigned char *nonce,
    const unsigned char *key
    )
{
    __m128i S[8], X[8];

    /* compute X and init state */
    LOAD_KEY(X, key, nonce, WORDS(STORM_N)/2);
    memcpy(S, X, 8 * sizeof(__m128i));

    /* absorb header */
    ABSORB_DATA(S, X, h, hlen);

    /* encrypt message */
    ENCRYPT_DATA(S, X, c, m, mlen);
    *clen = mlen + BYTES(STORM_T);

    /* extract tag */
    SQUEEZE_TAG(S, c + mlen);
}

int storm_aead_decrypt(
    unsigned char *m, size_t *mlen,
    const unsigned char *h, size_t hlen,
    const unsigned char *c, size_t clen,
    const unsigned char *nonce,
    const unsigned char *key
    )
{
    int result = -1;
    ALIGN(64) unsigned char tag[BYTES(STORM_T)];
    __m128i S[8], X[8];

    if (clen < BYTES(STORM_T))
        return result;

    /* compute X and init state */
    LOAD_KEY(X, key, nonce, WORDS(STORM_N)/2);
    memcpy(S, X, 8 * sizeof(__m128i));

    /* absorb header */
    ABSORB_DATA(S, X, h, hlen);

    /* decrypt message */
    DECRYPT_DATA(S, X, m, c, clen - BYTES(STORM_T));
    *mlen = clen - BYTES(STORM_T);

    /* extract tag */
    SQUEEZE_TAG(S, tag);

    /* verify tag */
    S[0] = _mm_cmpeq_epi8(LOADU(tag +                0), LOADU(c + clen - BYTES(STORM_T)  ));
    S[1] = _mm_cmpeq_epi8(LOADU(tag + BYTES(STORM_T)/2), LOADU(c + clen - BYTES(STORM_T)/2));
    result = (((unsigned long)_mm_movemask_epi8(AND(S[0], S[1])) + 1) >> 16) - 1;

    /* burn decrypted plaintext on authentication failure */
    if(result != 0)
    {
        burn(m, 0, *mlen);
    }

    return result;
}
