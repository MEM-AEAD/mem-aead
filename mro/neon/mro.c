#include "mro.h"
#include <string.h>
#include <arm_neon.h>

#define MRO_W 64           /* word size */
#define MRO_L 4            /* round number */
#define MRO_T (MRO_W *  4) /* tag size */
#define MRO_N (MRO_W *  2) /* nonce size */
#define MRO_K (MRO_W *  4) /* key size */
#define MRO_B (MRO_W * 16) /* permutation width */

/* constants for ROTR */
#define R0 32
#define R1 24
#define R2 16
#define R3 63

#define BYTES(X) (((X) + 7) / 8)
#define WORDS(X) (((X) + (MRO_W - 1)) / MRO_W)

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
#define SHL(X, C) vshlq_n_u64((X), (C))
#define SHR(X, C) vshrq_n_u64((X), (C))
#define ROT(X, C) XOR(SHR(X, C), SHL(X, (MRO_W)-(C)))
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

#define PERMUTE(S)             \
do                             \
{                              \
    int i;                     \
    for(i = 0; i < MRO_L; ++i) \
    {                          \
        F(S);                  \
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
    L[5] = COMB(vcreate_u64(MRO_L), vcreate_u64(MRO_T)); \
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
    ALIGN(32) unsigned char BLOCK[BYTES(MRO_B)]; \
    PAD(BLOCK, sizeof BLOCK, IN, INLEN);         \
    ABSORB_BLOCK(S, L, BLOCK);                   \
} while(0)

#define ENCRYPT_BLOCK(L, T, BLOCK_NR, OUT, IN)                     \
do                                                                 \
{                                                                  \
    uint64x2_t B[8];                                               \
    B[0] = XOR(L[0], T[0]);                                        \
    B[1] = XOR(L[1], T[1]);                                        \
    B[2] = L[2];                                                   \
    B[3] = L[3];                                                   \
    B[4] = L[4];                                                   \
    B[5] = L[5];                                                   \
    B[6] = L[6];                                                   \
    B[7] = XOR(L[7], COMB(vcreate_u64(0), vcreate_u64(BLOCK_NR))); \
    PERMUTE(B);                                                    \
    STOREU(OUT +   0, XOR(B[0], XOR(L[0], LOADU(IN +   0))));      \
    STOREU(OUT +  16, XOR(B[1], XOR(L[1], LOADU(IN +  16))));      \
    STOREU(OUT +  32, XOR(B[2], XOR(L[2], LOADU(IN +  32))));      \
    STOREU(OUT +  48, XOR(B[3], XOR(L[3], LOADU(IN +  48))));      \
    STOREU(OUT +  64, XOR(B[4], XOR(L[4], LOADU(IN +  64))));      \
    STOREU(OUT +  80, XOR(B[5], XOR(L[5], LOADU(IN +  80))));      \
    STOREU(OUT +  96, XOR(B[6], XOR(L[6], LOADU(IN +  96))));      \
    STOREU(OUT + 112, XOR(B[7], XOR(L[7], LOADU(IN + 112))));      \
} while(0)

#define ENCRYPT_LASTBLOCK(L, T, BLOCK_NR, OUT, IN, INLEN) \
do                                                        \
{                                                         \
    ALIGN(32) unsigned char BLOCK[BYTES(MRO_B)];          \
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
    while (l >= BYTES(MRO_B))                             \
    {                                                     \
        ABSORB_BLOCK(S, L, IN + i * BYTES(MRO_B));        \
        i += 1; l -= BYTES(MRO_B);                        \
        ALPHA(L);                                         \
    }                                                     \
    if (l > 0)                                            \
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
    while (l >= BYTES(MRO_B))                                                         \
    {                                                                                 \
        ENCRYPT_BLOCK(L, T, i, OUT + i * BYTES(MRO_B), IN + i * BYTES(MRO_B));        \
        i += 1; l -= BYTES(MRO_B);                                                    \
    }                                                                                 \
    if (l > 0)                                                                        \
    {                                                                                 \
        ENCRYPT_LASTBLOCK(L, T, i, OUT + i * BYTES(MRO_B), IN + i * BYTES(MRO_B), l); \
    }                                                                                 \
} while(0)

#define DECRYPT_DATA(L, T, OUT, IN, INLEN) \
do                                         \
{                                          \
    ENCRYPT_DATA(L, T, OUT, IN, INLEN);    \
} while(0)

#define FINALISE(S, L, HLEN, MLEN)                                \
do                                                                \
{                                                                 \
    BETA(L);                                                      \
    BETA(L);                                                      \
    S[0] = XOR(S[0], L[0]);                                       \
    S[1] = XOR(S[1], L[1]);                                       \
    S[2] = XOR(S[2], L[2]);                                       \
    S[3] = XOR(S[3], L[3]);                                       \
    S[4] = XOR(S[4], L[4]);                                       \
    S[5] = XOR(S[5], L[5]);                                       \
    S[6] = XOR(S[6], L[6]);                                       \
    S[7] = XOR(S[7], L[7]);                                       \
    S[7] = XOR(S[7], COMB(vcreate_u64(HLEN), vcreate_u64(MLEN))); \
    PERMUTE(S);                                                   \
    S[0] = XOR(S[0], L[0]);                                       \
    S[1] = XOR(S[1], L[1]);                                       \
    S[2] = XOR(S[2], L[2]);                                       \
    S[3] = XOR(S[3], L[3]);                                       \
    S[4] = XOR(S[4], L[4]);                                       \
    S[5] = XOR(S[5], L[5]);                                       \
    S[6] = XOR(S[6], L[6]);                                       \
    S[7] = XOR(S[7], L[7]);                                       \
} while(0)

static void* (* const volatile burn)(void*, int, size_t) = memset;

typedef enum tag__
{
    ABS_AD     = 0x00,
    ABS_MSG    = 0x01
} tag_t;

void mro_aead_encrypt(
    unsigned char *c, size_t *clen,
    const unsigned char *h, size_t hlen,
    const unsigned char *m, size_t mlen,
    const unsigned char *nonce,
    const unsigned char *key
    )
{
    uint64x2_t S[8] = {0};
    uint64x2_t LA[8] = {0};
    uint64x2_t LE[8] = {0};

    INIT_MASK(LE, key, nonce);

    /* absorb header and message */
    memcpy(LA, LE, 8 * sizeof(uint64x2_t));
    ABSORB_DATA(S, LA, h, hlen, ABS_AD);

    memcpy(LA, LE, 8 * sizeof(uint64x2_t));
    ABSORB_DATA(S, LA, m, mlen, ABS_MSG);

    memcpy(LA, LE, 8 * sizeof(uint64x2_t));
    FINALISE(S, LA, hlen, mlen);

    /* extract tag */
    STOREU(c + mlen, S[0]);
    STOREU(c + mlen + BYTES(MRO_T)/2, S[1]);
    *clen = mlen + BYTES(MRO_T);

    /* encrypt message */
    ENCRYPT_DATA(LE, S, c, m, mlen);
}

int mro_aead_decrypt(
    unsigned char *m, size_t *mlen,
    const unsigned char *h, size_t hlen,
    const unsigned char *c, size_t clen,
    const unsigned char *nonce,
    const unsigned char *key
    )
{
    int result = -1;
    uint64x2_t S[8] = {0};
    uint64x2_t LA[8] = {0};
    uint64x2_t LE[8] = {0};
    uint32x4_t T[2] = {0};

    if (clen < BYTES(MRO_T)) { return result; }

    INIT_MASK(LE, key, nonce);
    memcpy(LA, LE, 8 * sizeof(uint64x2_t));

    *mlen = clen - BYTES(MRO_T);

    /* store received tag temporarily in the first 2 state words */
    S[0] = LOADU(c + *mlen);
    S[1] = LOADU(c + *mlen + BYTES(MRO_T)/2);

    /* decrypt message */
    DECRYPT_DATA(LE, S, m, c, clen - BYTES(MRO_T));

    /* reset state */
    memset(S, 0, 8 * sizeof(uint64x2_t));

    /* absorb header and message */
    memcpy(LE, LA, 8 * sizeof(uint64x2_t));
    ABSORB_DATA(S, LA, h, hlen, ABS_AD);

    memcpy(LA, LE, 8 * sizeof(uint64x2_t));
    ABSORB_DATA(S, LA, m, *mlen, ABS_MSG);

    memcpy(LA, LE, 8 * sizeof(uint64x2_t));
    FINALISE(S, LA, hlen, *mlen);

    /* verify tag */
    T[0] = vceqq_u32( U64TOU32(S[0]), U8TOU32(vld1q_u8((uint8_t *)(c + clen - BYTES(MRO_T)  ))) );
    T[1] = vceqq_u32( U64TOU32(S[1]), U8TOU32(vld1q_u8((uint8_t *)(c + clen - BYTES(MRO_T)/2))) );
    T[0] = vandq_u32(T[0], T[1]);
    result = (0xFFFFFFFFFFFFFFFFULL == (vgetq_lane_u64(U32TOU64(T[0]), 0) & vgetq_lane_u64(U32TOU64(T[0]), 1)) ? 0 : -1);

    /* burn decrypted plaintext on authentication failure */
    if (result != 0) { burn(m, 0, *mlen); }

    return result;
}
