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

/* constants for ROTL */
#define L0 (STORM_W - R0)
#define L1 (STORM_W - R1)
#define L2 (STORM_W - R2)
#define L3 (STORM_W - R3)

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
#define SUB(X, Y) vsubq_u64((X), (Y))
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
    int i;                       \
    for(i = 0; i < STORM_R; ++i) \
    {                            \
        F(S);                    \
    }                            \
} while(0)

#define PERMUTE_INVERSE(S)       \
do                               \
{                                \
    int i;                       \
    for(i = 0; i < STORM_R; ++i) \
    {                            \
        FI(S);                   \
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

#define INIT(K, KEY, NONCE, TAG)                                \
do {                                                            \
    K[0] = LOADU(NONCE + 0);                                    \
    K[1] = COMBU64(vcreate_u64(0), vcreate_u64(0));             \
    K[2] = LOADU(KEY +  0);                                     \
    K[3] = LOADU(KEY + 16);                                     \
    K[4] = COMBU64(vcreate_u64(0), vcreate_u64(0));             \
    K[5] = COMBU64(vcreate_u64(0), vcreate_u64(0));             \
    K[6] = COMBU64(vcreate_u64(STORM_W), vcreate_u64(STORM_R)); \
    K[7] = COMBU64(vcreate_u64(STORM_T), vcreate_u64(TAG));     \
    PERMUTE(K);                                                 \
} while(0)

#define UPDATE(K)                                                                                                     \
do                                                                                                                    \
{                                                                                                                     \
    uint64x2_t T = XOR(ROT(COMBU64( LOU64(K[0]), vcreate_u64(0)), 9), SHR(COMBU64( HIU64(K[4]), vcreate_u64(0)), 7)); \
    K[0] = COMBU64( HIU64(K[0]), LOU64(K[1]));                                                                        \
    K[1] = COMBU64( HIU64(K[1]), LOU64(K[2]));                                                                        \
    K[2] = COMBU64( HIU64(K[2]), LOU64(K[3]));                                                                        \
    K[3] = COMBU64( HIU64(K[3]), LOU64(K[4]));                                                                        \
    K[4] = COMBU64( HIU64(K[4]), LOU64(K[5]));                                                                        \
    K[5] = COMBU64( HIU64(K[5]), LOU64(K[6]));                                                                        \
    K[6] = COMBU64( HIU64(K[6]), LOU64(K[7]));                                                                        \
    K[7] = COMBU64( HIU64(K[7]), LOU64(T   ));                                                                        \
} while(0)

#define ROTL256(K)                                             \
do                                                             \
{                                                              \
    uint64x2_t T;                                              \
    T = K[0]; K[0] = K[2]; K[2] = K[4]; K[4] = K[6]; K[6] = T; \
    T = K[1]; K[1] = K[3]; K[3] = K[5]; K[5] = K[7]; K[7] = T; \
} while(0)

#define ROTL512(K) \
do                 \
{                  \
    ROTL256(K);    \
    ROTL256(K);    \
} while(0)

#define ROTL768(K)                                             \
do                                                             \
{                                                              \
    uint64x2_t T;                                              \
    T = K[6]; K[6] = K[4]; K[4] = K[2]; K[2] = K[0]; K[0] = T; \
    T = K[7]; K[7] = K[5]; K[5] = K[3]; K[3] = K[1]; K[1] = T; \
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
    ROTL256(K);                                    \
    ABSORB_BLOCK(S, K, BLOCK);                     \
} while(0)


#define ENCRYPT_BLOCK(S, K, OUT, IN)    \
do                                      \
{                                       \
    uint64x2_t B[8];                    \
    B[0] = XOR(K[0], LOADU(IN +   0));  \
    B[1] = XOR(K[1], LOADU(IN +  16));  \
    B[2] = XOR(K[2], LOADU(IN +  32));  \
    B[3] = XOR(K[3], LOADU(IN +  48));  \
    B[4] = XOR(K[4], LOADU(IN +  64));  \
    B[5] = XOR(K[5], LOADU(IN +  80));  \
    B[6] = XOR(K[6], LOADU(IN +  96));  \
    B[7] = XOR(K[7], LOADU(IN + 112));  \
    PERMUTE(B);                         \
    STOREU(OUT +   0, XOR(K[0], B[0])); \
    STOREU(OUT +  16, XOR(K[1], B[1])); \
    STOREU(OUT +  32, XOR(K[2], B[2])); \
    STOREU(OUT +  48, XOR(K[3], B[3])); \
    STOREU(OUT +  64, XOR(K[4], B[4])); \
    STOREU(OUT +  80, XOR(K[5], B[5])); \
    STOREU(OUT +  96, XOR(K[6], B[6])); \
    STOREU(OUT + 112, XOR(K[7], B[7])); \
    S[0] = XOR(S[0], LOADU(IN +   0));  \
    S[1] = XOR(S[1], LOADU(IN +  16));  \
    S[2] = XOR(S[2], LOADU(IN +  32));  \
    S[3] = XOR(S[3], LOADU(IN +  48));  \
    S[4] = XOR(S[4], LOADU(IN +  64));  \
    S[5] = XOR(S[5], LOADU(IN +  80));  \
    S[6] = XOR(S[6], LOADU(IN +  96));  \
    S[7] = XOR(S[7], LOADU(IN + 112));  \
    UPDATE(K);                          \
} while(0)



#define ENCRYPT_LASTBLOCK(S, K, OUT, IN, INLEN)                    \
do                                                                 \
{                                                                  \
    uint64x2_t B[8];                                               \
    ALIGN(32) unsigned char BLOCK[BYTES(STORM_B)];                 \
    ROTL768(K);                                                    \
    B[0] = K[0];                                                   \
    B[1] = K[1];                                                   \
    B[2] = K[2];                                                   \
    B[3] = K[3];                                                   \
    B[4] = K[4];                                                   \
    B[5] = K[5];                                                   \
    B[6] = K[6];                                                   \
    B[7] = K[7];                                                   \
    PERMUTE(B);                                                    \
    PAD(BLOCK, sizeof BLOCK, IN, INLEN);                           \
    STOREU(BLOCK +   0, XOR(B[0], XOR(K[0], LOADU(BLOCK +   0)))); \
    STOREU(BLOCK +  16, XOR(B[1], XOR(K[1], LOADU(BLOCK +  16)))); \
    STOREU(BLOCK +  32, XOR(B[2], XOR(K[2], LOADU(BLOCK +  32)))); \
    STOREU(BLOCK +  48, XOR(B[3], XOR(K[3], LOADU(BLOCK +  48)))); \
    STOREU(BLOCK +  64, XOR(B[4], XOR(K[4], LOADU(BLOCK +  64)))); \
    STOREU(BLOCK +  80, XOR(B[5], XOR(K[5], LOADU(BLOCK +  80)))); \
    STOREU(BLOCK +  96, XOR(B[6], XOR(K[6], LOADU(BLOCK +  96)))); \
    STOREU(BLOCK + 112, XOR(B[7], XOR(K[7], LOADU(BLOCK + 112)))); \
    memcpy(OUT, BLOCK, INLEN);                                     \
    PAD(BLOCK, sizeof BLOCK, IN, INLEN);                           \
    S[0] = XOR(S[0], LOADU(BLOCK +   0));                          \
    S[1] = XOR(S[1], LOADU(BLOCK +  16));                          \
    S[2] = XOR(S[2], LOADU(BLOCK +  32));                          \
    S[3] = XOR(S[3], LOADU(BLOCK +  48));                          \
    S[4] = XOR(S[4], LOADU(BLOCK +  64));                          \
    S[5] = XOR(S[5], LOADU(BLOCK +  80));                          \
    S[6] = XOR(S[6], LOADU(BLOCK +  96));                          \
    S[7] = XOR(S[7], LOADU(BLOCK + 112));                          \
} while(0)


#define DECRYPT_BLOCK(S, K, OUT, IN)    \
do                                      \
{                                       \
    uint64x2_t B[8];                    \
    B[0] = XOR(K[0], LOADU(IN +   0));  \
    B[1] = XOR(K[1], LOADU(IN +  16));  \
    B[2] = XOR(K[2], LOADU(IN +  32));  \
    B[3] = XOR(K[3], LOADU(IN +  48));  \
    B[4] = XOR(K[4], LOADU(IN +  64));  \
    B[5] = XOR(K[5], LOADU(IN +  80));  \
    B[6] = XOR(K[6], LOADU(IN +  96));  \
    B[7] = XOR(K[7], LOADU(IN + 112));  \
    PERMUTE_INVERSE(B);                 \
    STOREU(OUT +   0, XOR(B[0], K[0])); \
    STOREU(OUT +  16, XOR(B[1], K[1])); \
    STOREU(OUT +  32, XOR(B[2], K[2])); \
    STOREU(OUT +  48, XOR(B[3], K[3])); \
    STOREU(OUT +  64, XOR(B[4], K[4])); \
    STOREU(OUT +  80, XOR(B[5], K[5])); \
    STOREU(OUT +  96, XOR(B[6], K[6])); \
    STOREU(OUT + 112, XOR(B[7], K[7])); \
    S[0] = XOR(S[0], LOADU(OUT +   0)); \
    S[1] = XOR(S[1], LOADU(OUT +  16)); \
    S[2] = XOR(S[2], LOADU(OUT +  32)); \
    S[3] = XOR(S[3], LOADU(OUT +  48)); \
    S[4] = XOR(S[4], LOADU(OUT +  64)); \
    S[5] = XOR(S[5], LOADU(OUT +  80)); \
    S[6] = XOR(S[6], LOADU(OUT +  96)); \
    S[7] = XOR(S[7], LOADU(OUT + 112)); \
    UPDATE(K);                          \
} while(0)


#define DECRYPT_LASTBLOCK(S, K, OUT, IN, INLEN)                    \
do                                                                 \
{                                                                  \
    uint64x2_t B[8];                                               \
    ALIGN(32) unsigned char BLOCK[BYTES(STORM_B)];                 \
    ROTL768(K);                                                    \
    B[0] = K[0];                                                   \
    B[1] = K[1];                                                   \
    B[2] = K[2];                                                   \
    B[3] = K[3];                                                   \
    B[4] = K[4];                                                   \
    B[5] = K[5];                                                   \
    B[6] = K[6];                                                   \
    B[7] = K[7];                                                   \
    PERMUTE(B);                                                    \
    PAD(BLOCK, sizeof BLOCK, IN, INLEN);                           \
    STOREU(BLOCK +   0, XOR(B[0], XOR(K[0], LOADU(BLOCK +   0)))); \
    STOREU(BLOCK +  16, XOR(B[1], XOR(K[1], LOADU(BLOCK +  16)))); \
    STOREU(BLOCK +  32, XOR(B[2], XOR(K[2], LOADU(BLOCK +  32)))); \
    STOREU(BLOCK +  48, XOR(B[3], XOR(K[3], LOADU(BLOCK +  48)))); \
    STOREU(BLOCK +  64, XOR(B[4], XOR(K[4], LOADU(BLOCK +  64)))); \
    STOREU(BLOCK +  80, XOR(B[5], XOR(K[5], LOADU(BLOCK +  80)))); \
    STOREU(BLOCK +  96, XOR(B[6], XOR(K[6], LOADU(BLOCK +  96)))); \
    STOREU(BLOCK + 112, XOR(B[7], XOR(K[7], LOADU(BLOCK + 112)))); \
    memcpy(OUT, BLOCK, INLEN);                                     \
    PAD(BLOCK, sizeof BLOCK, OUT, INLEN);                          \
    S[0] = XOR(S[0], LOADU(BLOCK +   0));                          \
    S[1] = XOR(S[1], LOADU(BLOCK +  16));                          \
    S[2] = XOR(S[2], LOADU(BLOCK +  32));                          \
    S[3] = XOR(S[3], LOADU(BLOCK +  48));                          \
    S[4] = XOR(S[4], LOADU(BLOCK +  64));                          \
    S[5] = XOR(S[5], LOADU(BLOCK +  80));                          \
    S[6] = XOR(S[6], LOADU(BLOCK +  96));                          \
    S[7] = XOR(S[7], LOADU(BLOCK + 112));                          \
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

#define ENCRYPT_DATA(S, K, OUT, IN, INLEN)                                             \
do                                                                                     \
{                                                                                      \
    size_t i = 0;                                                                      \
    size_t l = INLEN;                                                                  \
    while (l >= BYTES(STORM_B))                                                        \
    {                                                                                  \
        ENCRYPT_BLOCK(S, K, OUT + i * BYTES(STORM_B), IN + i * BYTES(STORM_B));        \
        i += 1; l -= BYTES(STORM_B);                                                   \
    }                                                                                  \
    if (l > 0)                                                                         \
    {                                                                                  \
        ENCRYPT_LASTBLOCK(S, K, OUT + i * BYTES(STORM_B), IN + i * BYTES(STORM_B), l); \
    }                                                                                  \
} while(0)

#define DECRYPT_DATA(S, K, OUT, IN, INLEN)                                             \
do                                                                                     \
{                                                                                      \
    size_t i = 0;                                                                      \
    size_t l = INLEN;                                                                  \
    while (l >= BYTES(STORM_B))                                                        \
    {                                                                                  \
        DECRYPT_BLOCK(S, K, OUT + i * BYTES(STORM_B), IN + i * BYTES(STORM_B));        \
        i += 1; l -= BYTES(STORM_B);                                                   \
    }                                                                                  \
    if (l > 0)                                                                         \
    {                                                                                  \
        DECRYPT_LASTBLOCK(S, K, OUT + i * BYTES(STORM_B), IN + i * BYTES(STORM_B), l); \
    }                                                                                  \
} while(0)

#define FINALISE(SA, SE, K)               \
do {                                      \
    ROTL512(K);                           \
    SE[0] = XOR(SE[0], K[0]);             \
    SE[1] = XOR(SE[1], K[1]);             \
    SE[2] = XOR(SE[2], K[2]);             \
    SE[3] = XOR(SE[3], K[3]);             \
    SE[4] = XOR(SE[4], K[4]);             \
    SE[5] = XOR(SE[5], K[5]);             \
    SE[6] = XOR(SE[6], K[6]);             \
    SE[7] = XOR(SE[7], K[7]);             \
    PERMUTE(SE);                          \
    SA[0] = XOR(SA[0], XOR(SE[0], K[0])); \
    SA[1] = XOR(SA[1], XOR(SE[1], K[1])); \
    SA[2] = XOR(SA[2], XOR(SE[2], K[2])); \
    SA[3] = XOR(SA[3], XOR(SE[3], K[3])); \
    SA[4] = XOR(SA[4], XOR(SE[4], K[4])); \
    SA[5] = XOR(SA[5], XOR(SE[5], K[5])); \
    SA[6] = XOR(SA[6], XOR(SE[6], K[6])); \
    SA[7] = XOR(SA[7], XOR(SE[7], K[7])); \
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
    uint64x2_t SA[8], SE[8], KA[8], KE[8];

    /* init states and keys */
    memset(SA, 0,  8 * sizeof(uint64x2_t));
    memset(SE, 0,  8 * sizeof(uint64x2_t));
    INIT(KA, key, nonce, ABS_TAG);
    INIT(KE, key, nonce, ENC_TAG);

    /* absorb header */
    ABSORB_DATA(SA, KA, h, hlen);

    /* encrypt message */
    ENCRYPT_DATA(SE, KE, c, m, mlen);
    *clen = mlen + BYTES(STORM_T);

    /* finalise */
    FINALISE(SA, SE, KA);

    /* extract tag */
    STOREU(c + mlen, SA[0]);
    STOREU(c + mlen + BYTES(STORM_T)/2, SA[1]);
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
    uint64x2_t SA[8], SE[8], KA[8], KE[8];
    uint32x4_t T[2];

    if (clen < BYTES(STORM_T)) { return result; }

    /* init states and keys */
    memset(SA, 0,  8 * sizeof(uint64x2_t));
    memset(SE, 0,  8 * sizeof(uint64x2_t));
    INIT(KA, key, nonce, ABS_TAG);
    INIT(KE, key, nonce, ENC_TAG);

    /* absorb header */
    ABSORB_DATA(SA, KA, h, hlen);

    /* decrypt message */
    DECRYPT_DATA(SE, KE, m, c, clen - BYTES(STORM_T));
    *mlen = clen - BYTES(STORM_T);

    /* finalise */
    FINALISE(SA, SE, KA);

    /* verify tag */
    T[0] = vceqq_u32( U64TOU32(SA[0]), U8TOU32(vld1q_u8((uint8_t *)(c + clen - BYTES(STORM_T)  ))) );
    T[1] = vceqq_u32( U64TOU32(SA[1]), U8TOU32(vld1q_u8((uint8_t *)(c + clen - BYTES(STORM_T)/2))) );
    T[0] = vandq_u32(T[0], T[1]);
    result = (0xFFFFFFFFFFFFFFFFULL == (vgetq_lane_u64(U32TOU64(T[0]), 0) & vgetq_lane_u64(U32TOU64(T[0]), 1)) ? 0 : -1);

    /* burn decrypted plaintext on authentication failure */
    if (result != 0) { burn(m, 0, *mlen); }

    return result;
}
