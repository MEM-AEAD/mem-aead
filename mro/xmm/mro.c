/*
    MRO - MEM AEAD source code package

    :copyright: (c) 2015 by Philipp Jovanovic and Samuel Neves
    :license: Creative Commons CC0 1.0
*/
#include <string.h>
#include "mro.h"

#if defined(_MSC_VER)
    #include <intrin.h>
#else
    #include <x86intrin.h>
#endif

#define MRO_W 64             /* word size */
#define MRO_L 4              /* round number */
#define MRO_T (MRO_W *  4) /* tag size */
#define MRO_N (MRO_W *  2) /* nonce size */
#define MRO_K (MRO_W *  4) /* key size */
#define MRO_B (MRO_W * 16) /* permutation width */

#define BYTES(x) (((x) + 7) / 8)
#define WORDS(x) (((x) + (MRO_W-1)) / MRO_W)

#if defined(_MSC_VER)
    #define ALIGN(x) __declspec(align(x))
#else
    #define ALIGN(x) __attribute__((aligned(x)))
#endif

#define XOR(A, B) _mm_xor_si128((A), (B))
#define AND(A, B) _mm_and_si128((A), (B))
#define ADD(A, B) _mm_add_epi64((A), (B))
#define ZERO _mm_setzero_si128()

#define R0 32
#define R1 24
#define R2 16
#define R3 63

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

#define PERMUTE(S)             \
do                             \
{                              \
    size_t i;                  \
    for(i = 0; i < MRO_L; ++i) \
    {                          \
        F(S);                  \
    }                          \
} while(0)

#define PAD(BLOCK, BLOCKLEN, IN, INLEN) \
do                                      \
{                                       \
    memset(BLOCK, 0, BLOCKLEN);         \
    memcpy(BLOCK, IN, INLEN);           \
    BLOCK[INLEN] = 0x01;                \
} while(0)

#define INIT_MASK(L, KEY, NONCE)         \
do                                       \
{                                        \
    L[0] = LOADU(NONCE + 0);             \
    L[1] = ZERO;                         \
    L[2] = ZERO;                         \
    L[3] = ZERO;                         \
    L[4] = ZERO;                         \
    L[5] = _mm_set_epi64x(MRO_T, MRO_L); \
    L[6] = LOADU(KEY +  0);              \
    L[7] = LOADU(KEY + 16);              \
    PERMUTE(L);                          \
} while(0)

#define ALPHA(L)                                                                                          \
do                                                                                                        \
{                                                                                                         \
    __m128i T = XOR(ROT(_mm_set_epi64x(0, L[0][0]), 11), _mm_slli_epi64(_mm_set_epi64x(0, L[2][1]), 13)); \
    L[0] = _mm_set_epi64x(L[1][0], L[0][1]);                                                              \
    L[1] = _mm_set_epi64x(L[2][0], L[1][1]);                                                              \
    L[2] = _mm_set_epi64x(L[3][0], L[2][1]);                                                              \
    L[3] = _mm_set_epi64x(L[4][0], L[3][1]);                                                              \
    L[4] = _mm_set_epi64x(L[5][0], L[4][1]);                                                              \
    L[5] = _mm_set_epi64x(L[6][0], L[5][1]);                                                              \
    L[6] = _mm_set_epi64x(L[7][0], L[6][1]);                                                              \
    L[7] = _mm_set_epi64x(T[0],    L[7][1]);                                                              \
} while(0)

#define BETA(L)                                                                                           \
do                                                                                                        \
{                                                                                                         \
    __m128i T = XOR(ROT(_mm_set_epi64x(0, L[0][0]), 11), _mm_slli_epi64(_mm_set_epi64x(0, L[2][1]), 13)); \
    L[0] = XOR(L[0], _mm_set_epi64x(L[1][0], L[0][1]));                                                   \
    L[1] = XOR(L[1], _mm_set_epi64x(L[2][0], L[1][1]));                                                   \
    L[2] = XOR(L[2], _mm_set_epi64x(L[3][0], L[2][1]));                                                   \
    L[3] = XOR(L[3], _mm_set_epi64x(L[4][0], L[3][1]));                                                   \
    L[4] = XOR(L[4], _mm_set_epi64x(L[5][0], L[4][1]));                                                   \
    L[5] = XOR(L[5], _mm_set_epi64x(L[6][0], L[5][1]));                                                   \
    L[6] = XOR(L[6], _mm_set_epi64x(L[7][0], L[6][1]));                                                   \
    L[7] = XOR(L[7], _mm_set_epi64x(T[0],    L[7][1]));                                                   \
} while(0)

#define GAMMA(L)                                                                                                      \
do                                                                                                                    \
{                                                                                                                     \
    __m128i T = XOR(ROT(_mm_set_epi64x(L[0][1], L[0][0]), 11), _mm_slli_epi64(_mm_set_epi64x(L[3][0], L[2][1]), 13)); \
    L[0] = XOR(L[0], XOR(L[1], _mm_set_epi64x(L[1][0], L[0][1])));                                                    \
    L[1] = XOR(L[1], XOR(L[2], _mm_set_epi64x(L[2][0], L[1][1])));                                                    \
    L[2] = XOR(L[2], XOR(L[3], _mm_set_epi64x(L[3][0], L[2][1])));                                                    \
    L[3] = XOR(L[3], XOR(L[4], _mm_set_epi64x(L[4][0], L[3][1])));                                                    \
    L[4] = XOR(L[4], XOR(L[5], _mm_set_epi64x(L[5][0], L[4][1])));                                                    \
    L[5] = XOR(L[5], XOR(L[6], _mm_set_epi64x(L[6][0], L[5][1])));                                                    \
    L[6] = XOR(L[6], XOR(L[7], _mm_set_epi64x(L[7][0], L[6][1])));                                                    \
    L[7] = XOR(L[7], XOR(T,    _mm_set_epi64x(T[0],    L[7][1])));                                                    \
} while(0)

#define ABSORB_BLOCK(S, L, IN)         \
do                                     \
{                                      \
    __m128i B[8];                      \
    B[0] = XOR(L[0], LOADU(IN +   0)); \
    B[1] = XOR(L[1], LOADU(IN +  16)); \
    B[2] = XOR(L[2], LOADU(IN +  32)); \
    B[3] = XOR(L[3], LOADU(IN +  48)); \
    B[4] = XOR(L[4], LOADU(IN +  64)); \
    B[5] = XOR(L[5], LOADU(IN +  80)); \
    B[6] = XOR(L[6], LOADU(IN +  96)); \
    B[7] = XOR(L[7], LOADU(IN + 112)); \
    PERMUTE(B);                        \
    S[0] = XOR(S[0], XOR(L[0], B[0])); \
    S[1] = XOR(S[1], XOR(L[1], B[1])); \
    S[2] = XOR(S[2], XOR(L[2], B[2])); \
    S[3] = XOR(S[3], XOR(L[3], B[3])); \
    S[4] = XOR(S[4], XOR(L[4], B[4])); \
    S[5] = XOR(S[5], XOR(L[5], B[5])); \
    S[6] = XOR(S[6], XOR(L[6], B[6])); \
    S[7] = XOR(S[7], XOR(L[7], B[7])); \
} while(0)

#define ABSORB_LASTBLOCK(S, L, IN, INLEN)        \
do                                               \
{                                                \
    ALIGN(64) unsigned char BLOCK[BYTES(MRO_B)]; \
    PAD(BLOCK, sizeof BLOCK, IN, INLEN);         \
    ABSORB_BLOCK(S, L, BLOCK);                   \
} while(0)

#define ENCRYPT_BLOCK(L, T, BLOCK_NR, OUT, IN)               \
do                                                           \
{                                                            \
    __m128i B[8];                                            \
    B[0] = XOR(L[0], T[0]);                                  \
    B[1] = XOR(L[1], T[1]);                                  \
    B[2] = L[2];                                             \
    B[3] = L[3];                                             \
    B[4] = L[4];                                             \
    B[5] = L[5];                                             \
    B[6] = L[6];                                             \
    B[7] = XOR(L[7], _mm_set_epi64x(BLOCK_NR, 0));           \
    PERMUTE(B);                                              \
    STORE(OUT +   0, XOR(B[0], XOR(L[0], LOADU(IN +   0)))); \
    STORE(OUT +  16, XOR(B[1], XOR(L[1], LOADU(IN +  16)))); \
    STORE(OUT +  32, XOR(B[2], XOR(L[2], LOADU(IN +  32)))); \
    STORE(OUT +  48, XOR(B[3], XOR(L[3], LOADU(IN +  48)))); \
    STORE(OUT +  64, XOR(B[4], XOR(L[4], LOADU(IN +  64)))); \
    STORE(OUT +  80, XOR(B[5], XOR(L[5], LOADU(IN +  80)))); \
    STORE(OUT +  96, XOR(B[6], XOR(L[6], LOADU(IN +  96)))); \
    STORE(OUT + 112, XOR(B[7], XOR(L[7], LOADU(IN + 112)))); \
} while(0)

#define ENCRYPT_LASTBLOCK(L, T, BLOCK_NR, OUT, IN, INLEN) \
do                                                        \
{                                                         \
    ALIGN(64) unsigned char BLOCK[BYTES(MRO_B)];          \
    memset(BLOCK, 0, BYTES(MRO_B));                       \
    memcpy(BLOCK, IN, INLEN);                             \
    ENCRYPT_BLOCK(L, T, BLOCK_NR, BLOCK, BLOCK);          \
    memcpy(OUT, BLOCK, INLEN);                            \
} while(0)

#define ABSORB_DATA(S, L, IN, INLEN, FLAG)                \
do                                                        \
{                                                         \
    size_t i = 0;                                         \
    size_t l = INLEN;                                     \
    if(FLAG)                                              \
    {                                                     \
        BETA(L);                                          \
    }                                                     \
    while(l >= BYTES(MRO_B))                              \
    {                                                     \
        ABSORB_BLOCK(S, L, IN + i * BYTES(MRO_B));        \
        i += 1; l -= BYTES(MRO_B);                        \
        ALPHA(L);                                         \
    }                                                     \
    if(l > 0)                                             \
    {                                                     \
        ABSORB_LASTBLOCK(S, L, IN + i * BYTES(MRO_B), l); \
    }                                                     \
} while(0)

#define ENCRYPT_DATA(L, T, OUT, IN, INLEN)                                            \
do                                                                                    \
{                                                                                     \
    size_t i = 0;                                                                     \
    size_t l = INLEN;                                                                 \
    GAMMA(L);                                                                         \
    while(l >= BYTES(MRO_B))                                                          \
    {                                                                                 \
        ENCRYPT_BLOCK(L, T, i, OUT + i * BYTES(MRO_B), IN + i * BYTES(MRO_B));        \
        i += 1; l -= BYTES(MRO_B);                                                    \
    }                                                                                 \
    if(l > 0)                                                                         \
    {                                                                                 \
        ENCRYPT_LASTBLOCK(L, T, i, OUT + i * BYTES(MRO_B), IN + i * BYTES(MRO_B), l); \
    }                                                                                 \
} while(0)

#define DECRYPT_DATA(L, T, OUT, IN, INLEN) \
do                                         \
{                                          \
    ENCRYPT_DATA(L, T, OUT, IN, INLEN);    \
} while(0)

#define FINALISE(S, L, HLEN, MLEN)                \
do                                                \
{                                                 \
    BETA(L);                                      \
    BETA(L);                                      \
    S[0] = XOR(S[0], L[0]);                       \
    S[1] = XOR(S[1], L[1]);                       \
    S[2] = XOR(S[2], L[2]);                       \
    S[3] = XOR(S[3], L[3]);                       \
    S[4] = XOR(S[4], L[4]);                       \
    S[5] = XOR(S[5], L[5]);                       \
    S[6] = XOR(S[6], L[6]);                       \
    S[7] = XOR(S[7], L[7]);                       \
    S[7] = XOR(S[7], _mm_set_epi64x(MLEN, HLEN)); \
    PERMUTE(S);                                   \
    S[0] = XOR(S[0], L[0]);                       \
    S[1] = XOR(S[1], L[1]);                       \
    S[2] = XOR(S[2], L[2]);                       \
    S[3] = XOR(S[3], L[3]);                       \
    S[4] = XOR(S[4], L[4]);                       \
    S[5] = XOR(S[5], L[5]);                       \
    S[6] = XOR(S[6], L[6]);                       \
    S[7] = XOR(S[7], L[7]);                       \
} while(0)

typedef enum tag__
{
    ABS_AD     = 0x00,
    ABS_MSG    = 0x01
} tag_t;

static void* (* const volatile burn)(void*, int, size_t) = memset;

void crypto_aead_encrypt(
    unsigned char *c, size_t *clen,
    const unsigned char *h, size_t hlen,
    const unsigned char *m, size_t mlen,
    const unsigned char *nonce,
    const unsigned char *key
    )
{
    __m128i S[8] = {0};
    __m128i LA[8] = {0};
    __m128i LE[8] = {0};

    INIT_MASK(LE, key, nonce);

    /* absorb header and message */
    memcpy(LA, LE, 8 * sizeof(__m128i));
    ABSORB_DATA(S, LA, h, hlen, ABS_AD);

    memcpy(LA, LE, 8 * sizeof(__m128i));
    ABSORB_DATA(S, LA, m, mlen, ABS_MSG);

    memcpy(LA, LE, 8 * sizeof(__m128i));
    FINALISE(S, LA, hlen, mlen);

    /* extract tag */
    STOREU(c + mlen, S[0]);
    STOREU(c + mlen + BYTES(MRO_T)/2, S[1]);
    *clen = mlen + BYTES(MRO_T);

    /* encrypt message */
    ENCRYPT_DATA(LE, S, c, m, mlen);
}


int crypto_aead_decrypt(
    unsigned char *m, size_t *mlen,
    const unsigned char *h, size_t hlen,
    const unsigned char *c, size_t clen,
    const unsigned char *nonce,
    const unsigned char *key
    )
{
    int result = -1;
    __m128i S[8] = {0};
    __m128i LA[8] = {0};
    __m128i LE[8] = {0};

    if (clen < BYTES(MRO_T)) { return result; }

    INIT_MASK(LE, key, nonce);
    memcpy(LA, LE, 8 * sizeof(__m128i));

    *mlen = clen - BYTES(MRO_T);

    /* store received tag temporarily in the first 2 state words */
    S[0] = LOADU(c + *mlen);
    S[1] = LOADU(c + *mlen + BYTES(MRO_T)/2);

    /* decrypt message */
    DECRYPT_DATA(LE, S, m, c, clen - BYTES(MRO_T));

    /* reset state */
    memset(S, 0, 8 * sizeof(__m128i));

    /* absorb header and message */
    memcpy(LE, LA, 8 * sizeof(__m128i));
    ABSORB_DATA(S, LA, h, hlen, ABS_AD);

    memcpy(LA, LE, 8 * sizeof(__m128i));
    ABSORB_DATA(S, LA, m, *mlen, ABS_MSG);

    memcpy(LA, LE, 8 * sizeof(__m128i));
    FINALISE(S, LA, hlen, *mlen);

    /* verify tag */
    S[0] = _mm_cmpeq_epi8(S[0], LOADU(c + clen - BYTES(MRO_T)  ));
    S[1] = _mm_cmpeq_epi8(S[1], LOADU(c + clen - BYTES(MRO_T)/2));
    result = (((_mm_movemask_epi8(AND(S[0], S[1])) & 0xFFFFUL) + 1) >> 16) - 1;

    /* burn decrypted plaintext on authentication failure */
    if (result != 0) { burn(m, 0, *mlen); }

    return result;
}
