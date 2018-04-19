/*
    OPP - MEM AEAD source code package

    :copyright: (c) 2015 by Philipp Jovanovic and Samuel Neves
    :license: Creative Commons CC0 1.0
*/
#include "opp.h"
#include <string.h>
#include <arm_neon.h>

#define OPP_W 64           /* word size */
#define OPP_L 4            /* round number */
#define OPP_T (OPP_W *  4) /* tag size */
#define OPP_N (OPP_W *  2) /* nonce size */
#define OPP_K (OPP_W *  4) /* key size */
#define OPP_B (OPP_W * 16) /* permutation width */

/* constants for ROTR */
#define R0 32
#define R1 24
#define R2 16
#define R3 63

/* constants for ROTL */
#define L0 (OPP_W - R0)
#define L1 (OPP_W - R1)
#define L2 (OPP_W - R2)
#define L3 (OPP_W - R3)

#define BYTES(X) (((X) + 7) / 8)
#define WORDS(X) (((X) + (OPP_W - 1)) / OPP_W)

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

#define LO(X) vget_low_u64((X))
#define HI(X) vget_high_u64((X))
#define COMB(X, Y) vcombine_u64((X), (Y))

#define XOR(X, Y) veorq_u64((X), (Y))
#define ADD(X, Y) vaddq_u64((X), (Y))
#define SUB(X, Y) vsubq_u64((X), (Y))
#define SHL(X, C) vshlq_n_u64((X), (C))
#define SHR(X, C) vshrq_n_u64((X), (C))
#define ROT(X, C) XOR(SHR(X, C), SHL(X, (OPP_W)-(C)))
#define ZERO vdupq_n_u64(0)

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
    S[2] = ROT(S[2],   L3);     S[3] = ROT(S[3],  L3); \
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

#define DIAGONALIZE(S)               \
do                                   \
{                                    \
    uint64x2_t T0, T1;               \
                                     \
    T0 = COMB( HI(S[2]), LO(S[3]) ); \
    T1 = COMB( HI(S[3]), LO(S[2]) ); \
    S[2] = T0;                       \
    S[3] = T1;                       \
                                     \
    T0 = S[4];                       \
    S[4] = S[5];                     \
    S[5] = T0;                       \
                                     \
    T0 = COMB( HI(S[6]), LO(S[7]) ); \
    T1 = COMB( HI(S[7]), LO(S[6]) ); \
    S[6] = T1;                       \
    S[7] = T0;                       \
} while(0)

#define UNDIAGONALIZE(S)             \
do                                   \
{                                    \
    uint64x2_t T0, T1;               \
                                     \
    T0 = COMB( HI(S[3]), LO(S[2]) ); \
    T1 = COMB( HI(S[2]), LO(S[3]) ); \
    S[2] = T0;                       \
    S[3] = T1;                       \
                                     \
    T0 = S[4];                       \
    S[4] = S[5];                     \
    S[5] = T0;                       \
                                     \
    T0 = COMB( HI(S[7]), LO(S[6]) ); \
    T1 = COMB( HI(S[6]), LO(S[7]) ); \
    S[6] = T1;                       \
    S[7] = T0;                       \
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

#define PERMUTE(S)             \
do                             \
{                              \
    int i;                     \
    for(i = 0; i < OPP_L; ++i) \
    {                          \
        F(S);                  \
    }                          \
} while(0)

#define PERMUTE_INVERSE(S)     \
do                             \
{                              \
    int i;                     \
    for(i = 0; i < OPP_L; ++i) \
    {                          \
        FI(S);                 \
    }                          \
} while(0)

#define PAD(OUT, OUTLEN, IN, INLEN) \
do                                  \
{                                   \
    memset(OUT, 0, OUTLEN);         \
    memcpy(OUT, IN, INLEN);         \
    OUT[INLEN] = 0x01;              \
} while(0)

#define INIT_MASK(L, KEY, NONCE)                         \
do                                                       \
{                                                        \
    L[0] = LOADU(NONCE + 0);                             \
    L[1] = ZERO;                                         \
    L[2] = ZERO;                                         \
    L[3] = ZERO;                                         \
    L[4] = ZERO;                                         \
    L[5] = COMB(vcreate_u64(OPP_L), vcreate_u64(OPP_T)); \
    L[6] = LOADU(KEY +  0);                              \
    L[7] = LOADU(KEY + 16);                              \
    PERMUTE(L);                                          \
} while(0)

#define ALPHA(L)                                                                \
do                                                                              \
{                                                                               \
    uint64x2_t T = XOR(ROT(L[0], 11), SHL(COMB(HI(L[2]), vcreate_u64(0)), 13)); \
    L[0] = COMB( HI(L[0]), LO(L[1]) );                                          \
    L[1] = COMB( HI(L[1]), LO(L[2]) );                                          \
    L[2] = COMB( HI(L[2]), LO(L[3]) );                                          \
    L[3] = COMB( HI(L[3]), LO(L[4]) );                                          \
    L[4] = COMB( HI(L[4]), LO(L[5]) );                                          \
    L[5] = COMB( HI(L[5]), LO(L[6]) );                                          \
    L[6] = COMB( HI(L[6]), LO(L[7]) );                                          \
    L[7] = COMB( HI(L[7]), LO(T   ) );                                          \
} while(0)

#define BETA(L)                                                                 \
do                                                                              \
{                                                                               \
    uint64x2_t T = XOR(ROT(L[0], 11), SHL(COMB(HI(L[2]), vcreate_u64(0)), 13)); \
    L[0] = XOR(L[0], COMB( HI(L[0]), LO(L[1]) ));                               \
    L[1] = XOR(L[1], COMB( HI(L[1]), LO(L[2]) ));                               \
    L[2] = XOR(L[2], COMB( HI(L[2]), LO(L[3]) ));                               \
    L[3] = XOR(L[3], COMB( HI(L[3]), LO(L[4]) ));                               \
    L[4] = XOR(L[4], COMB( HI(L[4]), LO(L[5]) ));                               \
    L[5] = XOR(L[5], COMB( HI(L[5]), LO(L[6]) ));                               \
    L[6] = XOR(L[6], COMB( HI(L[6]), LO(L[7]) ));                               \
    L[7] = XOR(L[7], COMB( HI(L[7]), LO(T   ) ));                               \
} while(0)

#define GAMMA(L)                                                          \
do                                                                        \
{                                                                         \
    uint64x2_t T = XOR(ROT(L[0], 11), SHL(COMB(HI(L[2]), LO(L[3])), 13)); \
    L[0] = XOR(L[0], XOR(L[1], COMB( HI(L[0]), LO(L[1]) )));              \
    L[1] = XOR(L[1], XOR(L[2], COMB( HI(L[1]), LO(L[2]) )));              \
    L[2] = XOR(L[2], XOR(L[3], COMB( HI(L[2]), LO(L[3]) )));              \
    L[3] = XOR(L[3], XOR(L[4], COMB( HI(L[3]), LO(L[4]) )));              \
    L[4] = XOR(L[4], XOR(L[5], COMB( HI(L[4]), LO(L[5]) )));              \
    L[5] = XOR(L[5], XOR(L[6], COMB( HI(L[5]), LO(L[6]) )));              \
    L[6] = XOR(L[6], XOR(L[7], COMB( HI(L[6]), LO(L[7]) )));              \
    L[7] = XOR(L[7], XOR(T,    COMB( HI(L[7]), LO(T   ) )));              \
} while(0)


#define ABSORB_BLOCK(S, L, IN)         \
do                                     \
{                                      \
    uint64x2_t B[8];                   \
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
    uint64x2_t B[8];                             \
    ALIGN(32) unsigned char BLOCK[BYTES(OPP_B)]; \
    PAD(BLOCK, sizeof BLOCK, IN, INLEN);         \
    B[0] = XOR(L[0], LOADU(BLOCK +   0));        \
    B[1] = XOR(L[1], LOADU(BLOCK +  16));        \
    B[2] = XOR(L[2], LOADU(BLOCK +  32));        \
    B[3] = XOR(L[3], LOADU(BLOCK +  48));        \
    B[4] = XOR(L[4], LOADU(BLOCK +  64));        \
    B[5] = XOR(L[5], LOADU(BLOCK +  80));        \
    B[6] = XOR(L[6], LOADU(BLOCK +  96));        \
    B[7] = XOR(L[7], LOADU(BLOCK + 112));        \
    PERMUTE(B);                                  \
    S[0] = XOR(S[0], XOR(L[0], B[0]));           \
    S[1] = XOR(S[1], XOR(L[1], B[1]));           \
    S[2] = XOR(S[2], XOR(L[2], B[2]));           \
    S[3] = XOR(S[3], XOR(L[3], B[3]));           \
    S[4] = XOR(S[4], XOR(L[4], B[4]));           \
    S[5] = XOR(S[5], XOR(L[5], B[5]));           \
    S[6] = XOR(S[6], XOR(L[6], B[6]));           \
    S[7] = XOR(S[7], XOR(L[7], B[7]));           \
} while(0)

#define ENCRYPT_BLOCK(S, L, OUT, IN)    \
do                                      \
{                                       \
    uint64x2_t B[8];                    \
    B[0] = XOR(L[0], LOADU(IN +   0));  \
    B[1] = XOR(L[1], LOADU(IN +  16));  \
    B[2] = XOR(L[2], LOADU(IN +  32));  \
    B[3] = XOR(L[3], LOADU(IN +  48));  \
    B[4] = XOR(L[4], LOADU(IN +  64));  \
    B[5] = XOR(L[5], LOADU(IN +  80));  \
    B[6] = XOR(L[6], LOADU(IN +  96));  \
    B[7] = XOR(L[7], LOADU(IN + 112));  \
    PERMUTE(B);                         \
    STOREU(OUT +   0, XOR(B[0], L[0])); \
    STOREU(OUT +  16, XOR(B[1], L[1])); \
    STOREU(OUT +  32, XOR(B[2], L[2])); \
    STOREU(OUT +  48, XOR(B[3], L[3])); \
    STOREU(OUT +  64, XOR(B[4], L[4])); \
    STOREU(OUT +  80, XOR(B[5], L[5])); \
    STOREU(OUT +  96, XOR(B[6], L[6])); \
    STOREU(OUT + 112, XOR(B[7], L[7])); \
    S[0] = XOR(S[0], LOADU(IN +   0));  \
    S[1] = XOR(S[1], LOADU(IN +  16));  \
    S[2] = XOR(S[2], LOADU(IN +  32));  \
    S[3] = XOR(S[3], LOADU(IN +  48));  \
    S[4] = XOR(S[4], LOADU(IN +  64));  \
    S[5] = XOR(S[5], LOADU(IN +  80));  \
    S[6] = XOR(S[6], LOADU(IN +  96));  \
    S[7] = XOR(S[7], LOADU(IN + 112));  \
} while(0)

#define ENCRYPT_LASTBLOCK(S, L, OUT, IN, INLEN)                    \
do                                                                 \
{                                                                  \
    uint64x2_t B[8];                                               \
    ALIGN(32) unsigned char BLOCK[BYTES(OPP_B)];                   \
    B[0] = L[0];                                                   \
    B[1] = L[1];                                                   \
    B[2] = L[2];                                                   \
    B[3] = L[3];                                                   \
    B[4] = L[4];                                                   \
    B[5] = L[5];                                                   \
    B[6] = L[6];                                                   \
    B[7] = L[7];                                                   \
    PERMUTE(B);                                                    \
    PAD(BLOCK, sizeof BLOCK, IN, INLEN);                           \
    STOREU(BLOCK +   0, XOR(B[0], XOR(L[0], LOADU(BLOCK +   0)))); \
    STOREU(BLOCK +  16, XOR(B[1], XOR(L[1], LOADU(BLOCK +  16)))); \
    STOREU(BLOCK +  32, XOR(B[2], XOR(L[2], LOADU(BLOCK +  32)))); \
    STOREU(BLOCK +  48, XOR(B[3], XOR(L[3], LOADU(BLOCK +  48)))); \
    STOREU(BLOCK +  64, XOR(B[4], XOR(L[4], LOADU(BLOCK +  64)))); \
    STOREU(BLOCK +  80, XOR(B[5], XOR(L[5], LOADU(BLOCK +  80)))); \
    STOREU(BLOCK +  96, XOR(B[6], XOR(L[6], LOADU(BLOCK +  96)))); \
    STOREU(BLOCK + 112, XOR(B[7], XOR(L[7], LOADU(BLOCK + 112)))); \
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

#define DECRYPT_BLOCK(S, L, OUT, IN)    \
do                                      \
{                                       \
    uint64x2_t B[8];                    \
    B[0] = XOR(L[0], LOADU(IN +   0));  \
    B[1] = XOR(L[1], LOADU(IN +  16));  \
    B[2] = XOR(L[2], LOADU(IN +  32));  \
    B[3] = XOR(L[3], LOADU(IN +  48));  \
    B[4] = XOR(L[4], LOADU(IN +  64));  \
    B[5] = XOR(L[5], LOADU(IN +  80));  \
    B[6] = XOR(L[6], LOADU(IN +  96));  \
    B[7] = XOR(L[7], LOADU(IN + 112));  \
    PERMUTE_INVERSE(B);                 \
    STOREU(OUT +   0, XOR(B[0], L[0])); \
    STOREU(OUT +  16, XOR(B[1], L[1])); \
    STOREU(OUT +  32, XOR(B[2], L[2])); \
    STOREU(OUT +  48, XOR(B[3], L[3])); \
    STOREU(OUT +  64, XOR(B[4], L[4])); \
    STOREU(OUT +  80, XOR(B[5], L[5])); \
    STOREU(OUT +  96, XOR(B[6], L[6])); \
    STOREU(OUT + 112, XOR(B[7], L[7])); \
    S[0] = XOR(S[0], LOADU(OUT +   0)); \
    S[1] = XOR(S[1], LOADU(OUT +  16)); \
    S[2] = XOR(S[2], LOADU(OUT +  32)); \
    S[3] = XOR(S[3], LOADU(OUT +  48)); \
    S[4] = XOR(S[4], LOADU(OUT +  64)); \
    S[5] = XOR(S[5], LOADU(OUT +  80)); \
    S[6] = XOR(S[6], LOADU(OUT +  96)); \
    S[7] = XOR(S[7], LOADU(OUT + 112)); \
} while(0)

#define DECRYPT_LASTBLOCK(S, L, OUT, IN, INLEN)                    \
do                                                                 \
{                                                                  \
    uint64x2_t B[8];                                               \
    ALIGN(32) unsigned char BLOCK[BYTES(OPP_B)];                   \
    B[0] = L[0];                                                   \
    B[1] = L[1];                                                   \
    B[2] = L[2];                                                   \
    B[3] = L[3];                                                   \
    B[4] = L[4];                                                   \
    B[5] = L[5];                                                   \
    B[6] = L[6];                                                   \
    B[7] = L[7];                                                   \
    PERMUTE(B);                                                    \
    PAD(BLOCK, sizeof BLOCK, IN, INLEN);                           \
    STOREU(BLOCK +   0, XOR(B[0], XOR(L[0], LOADU(BLOCK +   0)))); \
    STOREU(BLOCK +  16, XOR(B[1], XOR(L[1], LOADU(BLOCK +  16)))); \
    STOREU(BLOCK +  32, XOR(B[2], XOR(L[2], LOADU(BLOCK +  32)))); \
    STOREU(BLOCK +  48, XOR(B[3], XOR(L[3], LOADU(BLOCK +  48)))); \
    STOREU(BLOCK +  64, XOR(B[4], XOR(L[4], LOADU(BLOCK +  64)))); \
    STOREU(BLOCK +  80, XOR(B[5], XOR(L[5], LOADU(BLOCK +  80)))); \
    STOREU(BLOCK +  96, XOR(B[6], XOR(L[6], LOADU(BLOCK +  96)))); \
    STOREU(BLOCK + 112, XOR(B[7], XOR(L[7], LOADU(BLOCK + 112)))); \
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

#define ABSORB_DATA(S, L, IN, INLEN)                      \
do                                                        \
{                                                         \
    size_t i = 0;                                         \
    size_t l = INLEN;                                     \
    while (l >= BYTES(OPP_B))                             \
    {                                                     \
        ABSORB_BLOCK(S, L, IN + i * BYTES(OPP_B));        \
        i += 1;                                           \
        l -= BYTES(OPP_B);                                \
        ALPHA(L);                                         \
    }                                                     \
    if (l > 0)                                            \
    {                                                     \
        BETA(L);                                          \
        ABSORB_LASTBLOCK(S, L, IN + i * BYTES(OPP_B), l); \
    }                                                     \
} while(0)

#define ENCRYPT_DATA(S, L, OUT, IN, INLEN)                                         \
do                                                                                 \
{                                                                                  \
    size_t i = 0;                                                                  \
    size_t l = INLEN;                                                              \
    GAMMA(L);                                                                      \
    while (l >= BYTES(OPP_B))                                                      \
    {                                                                              \
        ENCRYPT_BLOCK(S, L, OUT + i * BYTES(OPP_B), IN + i * BYTES(OPP_B));        \
        i += 1;                                                                    \
        l -= BYTES(OPP_B);                                                         \
        ALPHA(L);                                                                  \
    }                                                                              \
    if (l > 0)                                                                     \
    {                                                                              \
        BETA(L);                                                                   \
        ENCRYPT_LASTBLOCK(S, L, OUT + i * BYTES(OPP_B), IN + i * BYTES(OPP_B), l); \
    }                                                                              \
} while(0)

#define DECRYPT_DATA(S, L, OUT, IN, INLEN)                                         \
do                                                                                 \
{                                                                                  \
    size_t i = 0;                                                                  \
    size_t l = INLEN;                                                              \
    GAMMA(L);                                                                      \
    while (l >= BYTES(OPP_B))                                                      \
    {                                                                              \
        DECRYPT_BLOCK(S, L, OUT + i * BYTES(OPP_B), IN + i * BYTES(OPP_B));        \
        i += 1;                                                                    \
        l -= BYTES(OPP_B);                                                         \
        ALPHA(L);                                                                  \
    }                                                                              \
    if (l > 0)                                                                     \
    {                                                                              \
        BETA(L);                                                                   \
        DECRYPT_LASTBLOCK(S, L, OUT + i * BYTES(OPP_B), IN + i * BYTES(OPP_B), l); \
    }                                                                              \
} while(0)

#define FINALISE(SA, SE, L)                                                     \
do                                                                              \
{                                                                               \
    size_t i;                                                                   \
    for(i = 0; i < 2; ++i)                                                      \
    {                                                                           \
       BETA(L);                                                                 \
    }                                                                           \
    SE[0] = XOR(SE[0], L[0]);                                                   \
    SE[1] = XOR(SE[1], L[1]);                                                   \
    SE[2] = XOR(SE[2], L[2]);                                                   \
    SE[3] = XOR(SE[3], L[3]);                                                   \
    SE[4] = XOR(SE[4], L[4]);                                                   \
    SE[5] = XOR(SE[5], L[5]);                                                   \
    SE[6] = XOR(SE[6], L[6]);                                                   \
    SE[7] = XOR(SE[7], L[7]);                                                   \
    PERMUTE(SE);                                                                \
    SA[0] = XOR(SA[0], XOR(SE[0], L[0]));                                       \
    SA[1] = XOR(SA[1], XOR(SE[1], L[1]));                                       \
    SA[2] = XOR(SA[2], XOR(SE[2], L[2]));                                       \
    SA[3] = XOR(SA[3], XOR(SE[3], L[3]));                                       \
    SA[4] = XOR(SA[4], XOR(SE[4], L[4]));                                       \
    SA[5] = XOR(SA[5], XOR(SE[5], L[5]));                                       \
    SA[6] = XOR(SA[6], XOR(SE[6], L[6]));                                       \
    SA[7] = XOR(SA[7], XOR(SE[7], L[7]));                                       \
} while(0)

static void* (* const volatile burn)(void*, int, size_t) = memset;

void crypto_aead_encrypt(
    unsigned char *c, size_t *clen,
    const unsigned char *h, size_t hlen,
    const unsigned char *m, size_t mlen,
    const unsigned char *nonce,
    const unsigned char *key
    )
{
    uint64x2_t SA[8] = {0};
    uint64x2_t SE[8] = {0};
    uint64x2_t LA[8] = {0};
    uint64x2_t LE[8] = {0};

    /* init states and masks */
    INIT_MASK(LA, key, nonce);
    memcpy(LE, LA, 8 * sizeof(uint64x2_t));

    /* absorb header */
    ABSORB_DATA(SA, LA, h, hlen);

    /* encrypt message */
    ENCRYPT_DATA(SE, LE, c, m, mlen);
    *clen = mlen + BYTES(OPP_T);

    /* finalise */
    FINALISE(SA, SE, LE);

    /* extract tag */
    STOREU(c + mlen, SA[0]);
    STOREU(c + mlen + BYTES(OPP_T)/2, SA[1]);
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
    uint64x2_t SA[8] = {0};
    uint64x2_t SE[8] = {0};
    uint64x2_t LA[8] = {0};
    uint64x2_t LE[8] = {0};
    uint32x4_t T[2] = {0};

    if (clen < BYTES(OPP_T)) { return result; }

    /* init states and masks */
    INIT_MASK(LA, key, nonce);
    memcpy(LE, LA, 8 * sizeof(uint64x2_t));

    /* absorb header */
    ABSORB_DATA(SA, LA, h, hlen);

    /* decrypt message */
    DECRYPT_DATA(SE, LE, m, c, clen - BYTES(OPP_T));
    *mlen = clen - BYTES(OPP_T);

    /* finalise */
    FINALISE(SA, SE, LE);

    /* verify tag */
    T[0] = vceqq_u32( U64TOU32(SA[0]), U8TOU32(vld1q_u8((uint8_t *)(c + clen - BYTES(OPP_T)  ))) );
    T[1] = vceqq_u32( U64TOU32(SA[1]), U8TOU32(vld1q_u8((uint8_t *)(c + clen - BYTES(OPP_T)/2))) );
    T[0] = vandq_u32(T[0], T[1]);
    result = (0xFFFFFFFFFFFFFFFFULL == (vgetq_lane_u64(U32TOU64(T[0]), 0) & vgetq_lane_u64(U32TOU64(T[0]), 1)) ? 0 : -1);

    /* burn decrypted plaintext on authentication failure */
    if (result != 0) { burn(m, 0, *mlen); }

    return result;
}
