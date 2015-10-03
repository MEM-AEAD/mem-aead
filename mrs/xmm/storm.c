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
#include <stdio.h>
#include "storm.h"

#if defined(_MSC_VER)
    #include <intrin.h>
#else
    #include <x86intrin.h>
#endif

#define STORM_W 64                  /* word size */
#define STORM_L 4                   /* round number */
#define STORM_T (STORM_W *  4)      /* tag size */
#define STORM_N (STORM_W *  2)      /* nonce size */
#define STORM_K (STORM_W *  4)      /* key size */
#define STORM_B (STORM_W * 16)      /* permutation width */
#define STORM_C (STORM_W *  4)      /* capacity */
#define STORM_R (STORM_B - STORM_C) /* rate */

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

#define R0 32
#define R1 24
#define R2 16
#define R3 63

#define LOAD(in) _mm_load_si128((__m128i*)(in))
#define STORE(out, x) _mm_store_si128((__m128i*)(out), (x))
#define LOADU(in) _mm_loadu_si128((__m128i*)(in))
#define STOREU(out, x) _mm_storeu_si128((__m128i*)(out), (x))

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

#define G(S)                                            \
do                                                      \
{                                                       \
    S[0] = ADD(S[0], S[2]);    S[1] = ADD(S[1], S[3]);  \
    S[6] = XOR(S[6], S[0]);    S[7] = XOR(S[7], S[1]);  \
    S[6] = ROT(S[6],   R0);    S[7] = ROT(S[7],   R0);  \
                                                        \
    S[4] = ADD(S[4], S[6]);    S[5] = ADD(S[5], S[7]);  \
    S[2] = XOR(S[2], S[4]);    S[3] = XOR(S[3], S[5]);  \
    S[2] = ROT(S[2],   R1);    S[3] = ROT(S[3],   R1);  \
                                                        \
    S[0] = ADD(S[0], S[2]);    S[1] = ADD(S[1], S[3]);  \
    S[6] = XOR(S[6], S[0]);    S[7] = XOR(S[7], S[1]);  \
    S[6] = ROT(S[6],   R2);    S[7] = ROT(S[7],   R2);  \
                                                        \
    S[4] = ADD(S[4], S[6]);    S[5] = ADD(S[5], S[7]);  \
    S[2] = XOR(S[2], S[4]);    S[3] = XOR(S[3], S[5]);  \
    S[2] = ROT(S[2],   R3);    S[3] = ROT(S[3],   R3);  \
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

#define PERMUTE(S)               \
do                               \
{                                \
    int r;                       \
    for(r = 0; r < STORM_L; ++r) \
    {                            \
        F(S);                    \
    }                            \
} while(0)

#define PAD(OUT, OUTLEN, IN, INLEN) \
do                                  \
{                                   \
    memset(OUT, 0, OUTLEN);         \
    memcpy(OUT, IN, INLEN);         \
} while(0)

#define ABSORB_BLOCK(S, IN)                                   \
do                                                            \
{                                                             \
    size_t j;                                                 \
    PERMUTE(S);                                               \
    for (j = 0; j < WORDS(STORM_B)/2; ++j)                    \
    {                                                         \
        S[j] = XOR(S[j], LOADU(IN + j * 2 * BYTES(STORM_W))); \
    }                                                         \
} while(0)


#define ABSORB_LASTBLOCK(S, IN, INLEN)                 \
do                                                     \
{                                                      \
    ALIGN(64) unsigned char lastblock[BYTES(STORM_B)]; \
    PAD(lastblock, sizeof lastblock, IN, INLEN);       \
    ABSORB_BLOCK(S, lastblock);                        \
} while(0)

#define ENCRYPT_BLOCK(S, OUT, IN)                             \
do                                                            \
{                                                             \
    size_t j;                                                 \
    PERMUTE(S);                                               \
    for (j = 0; j < WORDS(STORM_R)/2; ++j)                    \
    {                                                         \
        S[j] = XOR(S[j], LOADU(IN + j * 2 * BYTES(STORM_W))); \
        STOREU(OUT + j * 2 * BYTES(STORM_W), S[j]);           \
    }                                                         \
} while(0)

#define ENCRYPT_LASTBLOCK(S, OUT, IN, INLEN)           \
do                                                     \
{                                                      \
    ALIGN(64) unsigned char lastblock[BYTES(STORM_R)]; \
    PAD(lastblock, sizeof lastblock, IN, INLEN);       \
    ENCRYPT_BLOCK(S, lastblock, lastblock);            \
    memcpy(OUT, lastblock, INLEN);                     \
} while(0)

#define DECRYPT_BLOCK(S, OUT, IN)                           \
do                                                          \
{                                                           \
    size_t j;                                               \
    PERMUTE(S);                                             \
    for (j = 0; j < WORDS(STORM_R)/2; ++j)                  \
    {                                                       \
        __m128i T = LOADU(IN + j * 2 * BYTES(STORM_W));     \
        STOREU(OUT + j * 2 * BYTES(STORM_W), XOR(S[j], T)); \
        S[j] = T;                                           \
    }                                                       \
} while(0)

#define DECRYPT_LASTBLOCK(S, OUT, IN, INLEN)                      \
do                                                                \
{                                                                 \
    size_t j;                                                     \
    ALIGN(64) unsigned char lastblock[BYTES(STORM_R)];            \
    PERMUTE(S);                                                   \
    for (j = 0; j < WORDS(STORM_R)/2; ++j)                        \
    {                                                             \
        STOREU(lastblock + j * 2 * BYTES(STORM_W), S[j]);         \
    }                                                             \
    memcpy(lastblock, IN, INLEN);                                 \
    for (j = 0; j < WORDS(STORM_R)/2; ++j)                        \
    {                                                             \
        __m128i T = LOADU(lastblock + j * 2 * BYTES(STORM_W));    \
        STOREU(lastblock + j * 2 * BYTES(STORM_W), XOR(S[j], T)); \
        S[j] = T;                                                 \
    }                                                             \
    memcpy(OUT, lastblock, INLEN);                                \
} while(0)

#define INIT(S, KEY, IV, IVLEN, TAG)               \
do                                                 \
{                                                  \
    size_t j;                                      \
    memset(S, 0,  8 * sizeof(__m128i));            \
    for (j = 0; j < IVLEN; ++j)                    \
    {                                              \
        S[j] = LOADU(IV + j * 2 * BYTES(STORM_W)); \
    }                                              \
    S[4] = _mm_set_epi64x(STORM_L, 0);             \
    S[5] = _mm_set_epi64x(TAG, STORM_T);           \
    S[6] = LOADU(KEY + 0 * 2 * BYTES(STORM_W));    \
    S[7] = LOADU(KEY + 1 * 2 * BYTES(STORM_W));    \
} while(0)

#define ABSORB_DATA(S, IN, INLEN)                 \
do                                                \
{                                                 \
    size_t i = 0;                                 \
    size_t l = INLEN;                             \
    while (l >= BYTES(STORM_B))                   \
    {                                             \
        ABSORB_BLOCK(S, IN + i);                  \
        i += BYTES(STORM_B); l -= BYTES(STORM_B); \
    }                                             \
    if (l > 0)                                    \
    {                                             \
        ABSORB_LASTBLOCK(S, IN + i, l);           \
    }                                             \
} while(0)

#define ENCRYPT_DATA(S, OUT, IN, INLEN)           \
do                                                \
{                                                 \
    size_t i = 0;                                 \
    size_t l = INLEN;                             \
    while (l >= BYTES(STORM_R))                   \
    {                                             \
        ENCRYPT_BLOCK(S, OUT + i, IN + i);        \
        i += BYTES(STORM_R); l -= BYTES(STORM_R); \
    }                                             \
    if (l > 0)                                    \
    {                                             \
        ENCRYPT_LASTBLOCK(S, OUT + i, IN + i, l); \
    }                                             \
} while(0)

#define DECRYPT_DATA(S, OUT, IN, INLEN)           \
do                                                \
{                                                 \
    size_t i = 0;                                 \
    size_t l = INLEN;                             \
    while (l >= BYTES(STORM_R))                   \
    {                                             \
        DECRYPT_BLOCK(S, OUT + i, IN + i);        \
        i += BYTES(STORM_R), l -= BYTES(STORM_R); \
    }                                             \
    if (l > 0)                                    \
    {                                             \
        DECRYPT_LASTBLOCK(S, OUT + i, IN + i, l); \
    }                                             \
} while(0)

#define FINALISE(S, HLEN, MLEN, TAG)              \
do                                                \
{                                                 \
    PERMUTE(S);                                   \
    S[0] = XOR(S[0], _mm_set_epi64x(MLEN, HLEN)); \
    PERMUTE(S);                                   \
    STOREU(TAG,                    S[0]);         \
    STOREU(TAG + BYTES(STORM_T)/2, S[1]);         \
} while(0)

typedef enum tag__
{
    ABS_TAG = 0x00,
    ENC_TAG = 0x01
} tag_t;

static void* (* const volatile burn)(void*, int, size_t) = memset;

void storm_aead_encrypt(
    unsigned char *c, size_t *clen,
    const unsigned char *h, size_t hlen,
    const unsigned char *m, size_t mlen,
    const unsigned char *nonce,
    const unsigned char *key
    )
{
    __m128i S[8];

    /* absorption phase */
    INIT(S, key, nonce, WORDS(STORM_N)/2, ABS_TAG);
    ABSORB_DATA(S, h, hlen);
    ABSORB_DATA(S, m, mlen);
    FINALISE(S, hlen, mlen, c + mlen);
    *clen = mlen + BYTES(STORM_T);

    /* encryption phase */
    INIT(S, key, c + mlen, WORDS(STORM_T)/2, ENC_TAG);
    ENCRYPT_DATA(S, c, m, mlen);
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
    __m128i S[8];

    if (clen < BYTES(STORM_T)) { return -1; }

    /* decryption phase */
    INIT(S, key, c + clen - BYTES(STORM_T), WORDS(STORM_T)/2, ENC_TAG);
    DECRYPT_DATA(S, m, c, clen - BYTES(STORM_T));
    *mlen = clen - BYTES(STORM_T);

    /* absorption phase */
    INIT(S, key, nonce, WORDS(STORM_N)/2, ABS_TAG);
    ABSORB_DATA(S, h, hlen);
    ABSORB_DATA(S, m, *mlen);
    FINALISE(S, hlen, *mlen, tag);

    /* verify tag */
    S[0] = _mm_cmpeq_epi8(LOADU(tag +                0), LOADU(c + clen - BYTES(STORM_T)  ));
    S[1] = _mm_cmpeq_epi8(LOADU(tag + BYTES(STORM_T)/2), LOADU(c + clen - BYTES(STORM_T)/2));
    result = (((unsigned long)_mm_movemask_epi8(AND(S[0], S[1])) + 1) >> 16) - 1;

    /* burn decrypted plaintext on authentication failure */
    if(result != 0) { burn(m, 0, *mlen); }

    return result;
}
