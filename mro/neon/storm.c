/*
   STORM reference source code package - reference C implementations

   Written in 2015 by Philipp Jovanovic <philipp@jovanovic.io>

   To the extent possible under law, the author(s) have dedicated all copyright
   and related and neighboring rights to this software to the public domain
   worldwide. This software is distributed without any warranty.

   You should have received a copy of the CC0 Public Domain Dedication along with
   this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/
#include "storm.h"
#include "storm_config.h"
#include <string.h>
#include <arm_neon.h>

#define STORM_N (STORM_W *  2) /* nonce size */
#define STORM_K (STORM_W *  4) /* key size */
#define STORM_B (STORM_W * 16) /* permutation width */

/* constants for ROTR */
#define R0 32
#define R1 24
#define R2 16
#define R3 63

#define BYTES(X) (((X) + 7) / 8)
#define WORDS(X) (((X) + (STORM_W - 1)) / STORM_W)

#define ALIGN(X) __attribute((aligned(X)))

#define U32TOU64(X) vreinterpretq_u64_u32(X)
#define U64TOU32(X) vreinterpretq_u32_u64(X)
#define U8TOU32(X)  vreinterpretq_u32_u8(X)
#define U32TOU8(X)  vreinterpretq_u8_u32(X)
#define U8TOU64(X)  vreinterpretq_u64_u8(X)
#define U64TOU8(X)  vreinterpretq_u8_u64(X)

#define LOAD(IN) U8TOU64( vld1q_u8((uint8_t *)(IN)) )
#define STORE(OUT, IN) vst1q_u8( (uint8_t *)(OUT), U64TOU8(IN) )
#define LOADU(IN) LOAD(IN)
#define STOREU(OUT, IN) STORE(OUT, IN)

#define LOU64(X) vget_low_u64((X))
#define HIU64(X) vget_high_u64((X))
#define COMBU64(X, Y) vcombine_u64((X), (Y))

#define XOR(X, Y) veorq_u64((X), (Y))
#define ADD(X, Y) vaddq_u64((X), (Y))
#define SHL(X, C) vshlq_n_u64((X), (C))
#define SHR(X, C) vshrq_n_u64((X), (C))
#define ROT(X, C) XOR(SHR(X, C), SHL(X, (STORM_W)-(C)))

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

#define DIAGONALIZE(S)                        \
do                                            \
{                                             \
    uint64x2_t T0, T1;                        \
                                              \
    T0 = COMBU64( HIU64(S[2]), LOU64(S[3]) ); \
    T1 = COMBU64( HIU64(S[3]), LOU64(S[2]) ); \
    S[2] = T0;                                \
    S[3] = T1;                                \
                                              \
    T0 = S[4];                                \
    S[4] = S[5];                              \
    S[5] = T0;                                \
                                              \
    T0 = COMBU64( HIU64(S[6]), LOU64(S[7]) ); \
    T1 = COMBU64( HIU64(S[7]), LOU64(S[6]) ); \
    S[6] = T1;                                \
    S[7] = T0;                                \
} while(0)

#define UNDIAGONALIZE(S)                      \
do                                            \
{                                             \
    uint64x2_t T0, T1;                        \
                                              \
    T0 = COMBU64( HIU64(S[3]), LOU64(S[2]) ); \
    T1 = COMBU64( HIU64(S[2]), LOU64(S[3]) ); \
    S[2] = T0;                                \
    S[3] = T1;                                \
                                              \
    T0 = S[4];                                \
    S[4] = S[5];                              \
    S[5] = T0;                                \
                                              \
    T0 = COMBU64( HIU64(S[7]), LOU64(S[6]) ); \
    T1 = COMBU64( HIU64(S[6]), LOU64(S[7]) ); \
    S[6] = T1;                                \
    S[7] = T0;                                \
} while(0)

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
    int i;                       \
    for(i = 0; i < STORM_R; ++i) \
    {                            \
        F(S);                    \
    }                            \
} while(0)

#define PAD(OUT, OUTLEN, IN, INLEN) \
do                                  \
{                                   \
    memset(OUT, 0, OUTLEN);         \
    memcpy(OUT, IN, INLEN);         \
    OUT[INLEN] = 0x01;              \
    OUT[OUTLEN - 1] |= 0x80;        \
} while(0)

#define INIT(K, KEY, IV, IVLEN, TAG)                            \
do {                                                            \
    size_t i;                                                   \
    K[0] = COMBU64(vcreate_u64(0), vcreate_u64(0));             \
    K[1] = COMBU64(vcreate_u64(0), vcreate_u64(0));             \
    for (i = 0; i < IVLEN; ++i)                                 \
    {                                                           \
        K[i] = LOADU(IV + i * 2 * BYTES(STORM_W));              \
    }                                                           \
    K[2] = LOADU(KEY +  0);                                     \
    K[3] = LOADU(KEY + 16);                                     \
    K[4] = COMBU64(vcreate_u64(0), vcreate_u64(0));             \
    K[5] = COMBU64(vcreate_u64(0), vcreate_u64(0));             \
    K[6] = COMBU64(vcreate_u64(STORM_W), vcreate_u64(STORM_R)); \
    K[7] = COMBU64(vcreate_u64(STORM_T), vcreate_u64(TAG));     \
    PERMUTE(K);                                                 \
} while(0)

#define UPDATE(K)                                                                                                       \
do                                                                                                                      \
{                                                                                                                       \
    uint64x2_t T = XOR(ROT(COMBU64( LOU64(K[0]), vcreate_u64(0)), 11), SHL(COMBU64( HIU64(K[2]), vcreate_u64(0)), 13)); \
    K[0] = COMBU64( HIU64(K[0]), LOU64(K[1]));                                                                          \
    K[1] = COMBU64( HIU64(K[1]), LOU64(K[2]));                                                                          \
    K[2] = COMBU64( HIU64(K[2]), LOU64(K[3]));                                                                          \
    K[3] = COMBU64( HIU64(K[3]), LOU64(K[4]));                                                                          \
    K[4] = COMBU64( HIU64(K[4]), LOU64(K[5]));                                                                          \
    K[5] = COMBU64( HIU64(K[5]), LOU64(K[6]));                                                                          \
    K[6] = COMBU64( HIU64(K[6]), LOU64(K[7]));                                                                          \
    K[7] = COMBU64( HIU64(K[7]), LOU64(T   ));                                                                          \
} while(0)

#define ABSORB_BLOCK(S, K, IN)         \
do                                     \
{                                      \
    uint64x2_t B[8];                   \
    B[0] = XOR(K[0], LOADU(IN +   0)); \
    B[1] = XOR(K[1], LOADU(IN +  16)); \
    B[2] = XOR(K[2], LOADU(IN +  32)); \
    B[3] = XOR(K[3], LOADU(IN +  48)); \
    B[4] = XOR(K[4], LOADU(IN +  64)); \
    B[5] = XOR(K[5], LOADU(IN +  80)); \
    B[6] = XOR(K[6], LOADU(IN +  96)); \
    B[7] = XOR(K[7], LOADU(IN + 112)); \
    PERMUTE(B);                        \
    S[0] = XOR(S[0], XOR(K[0], B[0])); \
    S[1] = XOR(S[1], XOR(K[1], B[1])); \
    S[2] = XOR(S[2], XOR(K[2], B[2])); \
    S[3] = XOR(S[3], XOR(K[3], B[3])); \
    S[4] = XOR(S[4], XOR(K[4], B[4])); \
    S[5] = XOR(S[5], XOR(K[5], B[5])); \
    S[6] = XOR(S[6], XOR(K[6], B[6])); \
    S[7] = XOR(S[7], XOR(K[7], B[7])); \
    UPDATE(K);                         \
} while(0)

#define ABSORB_LASTBLOCK(S, K, IN, INLEN)          \
do                                                 \
{                                                  \
    ALIGN(32) unsigned char BLOCK[BYTES(STORM_B)]; \
    PAD(BLOCK, sizeof BLOCK, IN, INLEN);           \
    ABSORB_BLOCK(S, K, BLOCK);                     \
} while(0)

#define ABSORB_FINALISE(S, K, HLEN, MLEN)                            \
do                                                                   \
{                                                                    \
    uint64x2_t B[8];                                                 \
    B[0] = K[0];                                                     \
    B[1] = K[1];                                                     \
    B[2] = K[2];                                                     \
    B[3] = K[3];                                                     \
    B[4] = K[4];                                                     \
    B[5] = K[5];                                                     \
    B[6] = K[6];                                                     \
    B[7] = XOR(K[7], COMBU64(vcreate_u64(HLEN), vcreate_u64(MLEN))); \
    PERMUTE(B);                                                      \
    S[0] = XOR(S[0], XOR(K[0], B[0]));                               \
    S[1] = XOR(S[1], XOR(K[1], B[1]));                               \
    S[2] = XOR(S[2], XOR(K[2], B[2]));                               \
    S[3] = XOR(S[3], XOR(K[3], B[3]));                               \
    S[4] = XOR(S[4], XOR(K[4], B[4]));                               \
    S[5] = XOR(S[5], XOR(K[5], B[5]));                               \
    S[6] = XOR(S[6], XOR(K[6], B[6]));                               \
    S[7] = XOR(S[7], XOR(K[7], B[7]));                               \
    UPDATE(K);                                                       \
} while(0)

#define ENCRYPT_BLOCK(K, BLOCK_NR, OUT, IN)                           \
do                                                                    \
{                                                                     \
    uint64x2_t B[8];                                                  \
    B[0] = K[0];                                                      \
    B[1] = K[1];                                                      \
    B[2] = K[2];                                                      \
    B[3] = K[3];                                                      \
    B[4] = K[4];                                                      \
    B[5] = K[5];                                                      \
    B[6] = K[6];                                                      \
    B[7] = XOR(K[7], COMBU64(vcreate_u64(0), vcreate_u64(BLOCK_NR))); \
    PERMUTE(B);                                                       \
    STOREU(OUT +   0, XOR(B[0], XOR(K[0], LOADU(IN +   0))));         \
    STOREU(OUT +  16, XOR(B[1], XOR(K[1], LOADU(IN +  16))));         \
    STOREU(OUT +  32, XOR(B[2], XOR(K[2], LOADU(IN +  32))));         \
    STOREU(OUT +  48, XOR(B[3], XOR(K[3], LOADU(IN +  48))));         \
    STOREU(OUT +  64, XOR(B[4], XOR(K[4], LOADU(IN +  64))));         \
    STOREU(OUT +  80, XOR(B[5], XOR(K[5], LOADU(IN +  80))));         \
    STOREU(OUT +  96, XOR(B[6], XOR(K[6], LOADU(IN +  96))));         \
    STOREU(OUT + 112, XOR(B[7], XOR(K[7], LOADU(IN + 112))));         \
} while(0)

#define ENCRYPT_LASTBLOCK(K, BLOCK_NR, OUT, IN, INLEN) \
do                                                     \
{                                                      \
    ALIGN(32) unsigned char BLOCK[BYTES(STORM_B)];     \
    memset(BLOCK, 0, BYTES(STORM_B));                  \
    memcpy(BLOCK, IN, INLEN);                          \
    ENCRYPT_BLOCK(K, BLOCK_NR, BLOCK, BLOCK);          \
    memcpy(OUT, BLOCK, INLEN);                         \
} while(0)

#define ABSORB_DATA(S, K, IN, INLEN)                        \
do                                                          \
{                                                           \
    size_t i = 0;                                           \
    size_t l = INLEN;                                       \
    while (l >= BYTES(STORM_B))                             \
    {                                                       \
        ABSORB_BLOCK(S, K, IN + i * BYTES(STORM_B));        \
        i += 1; l -= BYTES(STORM_B);                        \
    }                                                       \
    if (l > 0)                                              \
    {                                                       \
        ABSORB_LASTBLOCK(S, K, IN + i * BYTES(STORM_B), l); \
    }                                                       \
} while(0)

#define ENCRYPT_DATA(K, OUT, IN, INLEN)                                                \
do                                                                                     \
{                                                                                      \
    size_t i = 0;                                                                      \
    size_t l = INLEN;                                                                  \
    while (l >= BYTES(STORM_B))                                                        \
    {                                                                                  \
        ENCRYPT_BLOCK(K, i, OUT + i * BYTES(STORM_B), IN + i * BYTES(STORM_B));        \
        i += 1; l -= BYTES(STORM_B);                                                   \
    }                                                                                  \
    if (l > 0)                                                                         \
    {                                                                                  \
        ENCRYPT_LASTBLOCK(K, i, OUT + i * BYTES(STORM_B), IN + i * BYTES(STORM_B), l); \
    }                                                                                  \
} while(0)

#define DECRYPT_DATA(K, OUT, IN, INLEN) \
do                                      \
{                                       \
    ENCRYPT_DATA(K, OUT, IN, INLEN);    \
} while(0)

static void* (* const volatile burn)(void*, int, size_t) = memset;

typedef enum tag__
{
    ABS_TAG     = 0x00,
    ENC_TAG     = 0x01
} tag_t;

void storm_aead_encrypt(
    unsigned char *c, size_t *clen,
    const unsigned char *h, size_t hlen,
    const unsigned char *m, size_t mlen,
    const unsigned char *nonce,
    const unsigned char *key
    )
{
    uint64x2_t S[8], K[8];

    /* absorb header and message */
    memset(S, 0,  8 * sizeof(uint64x2_t));
    INIT(K, key, nonce, WORDS(STORM_N)/2, ABS_TAG);
    ABSORB_DATA(S, K, h, hlen);
    ABSORB_DATA(S, K, m, mlen);
    ABSORB_FINALISE(S, K, hlen, mlen);

    /* extract tag */
    STOREU(c + mlen, S[0]);
    STOREU(c + mlen + BYTES(STORM_T)/2, S[1]);
    *clen = mlen + BYTES(STORM_T);

    /* encrypt message */
    INIT(K, key, c + mlen, WORDS(STORM_T)/2, ENC_TAG);
    ENCRYPT_DATA(K, c, m, mlen);
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
    uint64x2_t S[8], K[8];
    uint32x4_t T[2];

    if (clen < BYTES(STORM_T)) { return result; }

    /* decrypt message */
    INIT(K, key, c + clen - BYTES(STORM_T), WORDS(STORM_T)/2, ENC_TAG);
    DECRYPT_DATA(K, m, c, clen - BYTES(STORM_T));
    *mlen = clen - BYTES(STORM_T);

    /* absorb header and message */
    memset(S, 0,  8 * sizeof(uint64x2_t));
    INIT(K, key, nonce, WORDS(STORM_N)/2, ABS_TAG);
    ABSORB_DATA(S, K, h, hlen);
    ABSORB_DATA(S, K, m, *mlen);
    ABSORB_FINALISE(S, K, hlen, *mlen);

    /* verify tag */
    T[0] = vceqq_u32( U64TOU32(S[0]), U8TOU32(vld1q_u8((uint8_t *)(c + clen - BYTES(STORM_T)  ))) );
    T[1] = vceqq_u32( U64TOU32(S[1]), U8TOU32(vld1q_u8((uint8_t *)(c + clen - BYTES(STORM_T)/2))) );
    T[0] = vandq_u32(T[0], T[1]);
    result = (0xFFFFFFFFFFFFFFFFULL == (vgetq_lane_u64(U32TOU64(T[0]), 0) & vgetq_lane_u64(U32TOU64(T[0]), 1)) ? 0 : -1);

    /* burn decrypted plaintext on authentication failure */
    if (result != 0) { burn(m, 0, *mlen); }

    return result;
}
