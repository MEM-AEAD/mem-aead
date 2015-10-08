#include "mrs.h"
#include <string.h>
#include <arm_neon.h>

#define MRS_W 64              /* word size */
#define MRS_L 4               /* round number */
#define MRS_T (MRS_W *  4)    /* tag size */
#define MRS_N (MRS_W *  2)    /* nonce size */
#define MRS_K (MRS_W *  4)    /* key size */
#define MRS_B (MRS_W * 16)    /* permutation width */
#define MRS_C (MRS_W *  4)    /* capacity */
#define MRS_R (MRS_B - MRS_C) /* rate */

/* constants for ROTR */
#define R0 32
#define R1 24
#define R2 16
#define R3 63

#define BYTES(X) (((X) + 7) / 8)
#define WORDS(X) (((X) + (MRS_W - 1)) / MRS_W)

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
#define ROT(X, C) XOR(SHR(X, C), SHL(X, (MRS_W)-(C)))

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
    for(i = 0; i < MRS_L; ++i) \
    {                          \
        F(S);                  \
    }                          \
} while(0)

#define PAD(OUT, OUTLEN, IN, INLEN) \
do                                  \
{                                   \
    memset(OUT, 0, OUTLEN);         \
    memcpy(OUT, IN, INLEN);         \
} while(0)

#define ABSORB_BLOCK(S, IN)                                 \
do                                                          \
{                                                           \
    size_t j;                                               \
    PERMUTE(S);                                             \
    for (j = 0; j < WORDS(MRS_B)/2; ++j)                    \
    {                                                       \
        S[j] = XOR(S[j], LOADU(IN + j * 2 * BYTES(MRS_W))); \
    }                                                       \
} while(0)

#define ABSORB_LASTBLOCK(S, IN, INLEN)           \
do                                               \
{                                                \
    ALIGN(32) unsigned char BLOCK[BYTES(MRS_B)]; \
    PAD(BLOCK, sizeof BLOCK, IN, INLEN);         \
    ABSORB_BLOCK(S, BLOCK);                      \
} while(0)

#define ENCRYPT_BLOCK(S, OUT, IN)                           \
do                                                          \
{                                                           \
    size_t j;                                               \
    PERMUTE(S);                                             \
    for (j = 0; j < WORDS(MRS_R)/2; ++j)                    \
    {                                                       \
        S[j] = XOR(S[j], LOADU(IN + j * 2 * BYTES(MRS_W))); \
        STOREU(OUT + j * 2 * BYTES(MRS_W), S[j]);           \
    }                                                       \
} while(0)

#define ENCRYPT_LASTBLOCK(S, OUT, IN, INLEN)     \
do                                               \
{                                                \
    ALIGN(32) unsigned char BLOCK[BYTES(MRS_R)]; \
    PAD(BLOCK, sizeof BLOCK, IN, INLEN);         \
    ENCRYPT_BLOCK(S, BLOCK, BLOCK);              \
    memcpy(OUT, BLOCK, INLEN);                   \
} while(0)

#define DECRYPT_BLOCK(S, OUT, IN)                         \
do                                                        \
{                                                         \
    size_t j;                                             \
    PERMUTE(S);                                           \
    for (j = 0; j < WORDS(MRS_R)/2; ++j)                  \
    {                                                     \
        uint64x2_t T = LOADU(IN + j * 2 * BYTES(MRS_W));  \
        STOREU(OUT + j * 2 * BYTES(MRS_W), XOR(S[j], T)); \
        S[j] = T;                                         \
    }                                                     \
} while(0)

#define DECRYPT_LASTBLOCK(S, OUT, IN, INLEN)                \
do                                                          \
{                                                           \
    size_t j;                                               \
    ALIGN(32) unsigned char BLOCK[BYTES(MRS_R)];            \
    PERMUTE(S);                                             \
    for (j = 0; j < WORDS(MRS_R)/2; ++j)                    \
    {                                                       \
        STOREU(BLOCK + j * 2 * BYTES(MRS_W), S[j]);         \
    }                                                       \
    memcpy(BLOCK, IN, INLEN);                               \
    for (j = 0; j < WORDS(MRS_R)/2; ++j)                    \
    {                                                       \
        uint64x2_t T = LOADU(BLOCK + j * 2 * BYTES(MRS_W)); \
        STOREU(BLOCK + j * 2 * BYTES(MRS_W), XOR(S[j], T)); \
        S[j] = T;                                           \
    }                                                       \
    memcpy(OUT, BLOCK, INLEN);                              \
} while(0)

#define INIT(S, KEY, IV, IVLEN, TAG)                   \
do {                                                   \
    size_t j;                                          \
    memset(S, 0, 8 * sizeof(uint64x2_t));              \
    for (j = 0; j < IVLEN; ++j)                        \
    {                                                  \
        S[j] = LOADU(IV + j * 2 * BYTES(MRS_W));       \
    }                                                  \
    S[4] = COMB(vcreate_u64(0), vcreate_u64(MRS_L));   \
    S[5] = COMB(vcreate_u64(MRS_T), vcreate_u64(TAG)); \
    S[6] = LOADU(KEY + 0 * 2 * BYTES(MRS_W));          \
    S[7] = LOADU(KEY + 1 * 2 * BYTES(MRS_W));          \
} while(0)

#define ABSORB_DATA(S, IN, INLEN)             \
do                                            \
{                                             \
    size_t i = 0;                             \
    size_t l = INLEN;                         \
    while (l >= BYTES(MRS_B))                 \
    {                                         \
        ABSORB_BLOCK(S, IN + i);              \
        i += BYTES(MRS_B); l -= BYTES(MRS_B); \
    }                                         \
    if (l > 0)                                \
    {                                         \
        ABSORB_LASTBLOCK(S, IN + i, l);       \
    }                                         \
} while(0)

#define ENCRYPT_DATA(S, OUT, IN, INLEN)           \
do                                                \
{                                                 \
   size_t i = 0;                                  \
   size_t l = INLEN;                              \
   while (l >= BYTES(MRS_R))                      \
   {                                              \
       ENCRYPT_BLOCK(S, OUT + i, IN + i);         \
       i += BYTES(MRS_R); l -= BYTES(MRS_R);      \
   }                                              \
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
    while (l >= BYTES(MRS_R))                     \
    {                                             \
        DECRYPT_BLOCK(S, OUT + i, IN + i);        \
        i += BYTES(MRS_R); l -= BYTES(MRS_R);     \
    }                                             \
    if (l > 0)                                    \
    {                                             \
        DECRYPT_LASTBLOCK(S, OUT + i, IN + i, l); \
    }                                             \
} while(0)

#define FINALISE(S, HLEN, MLEN)                                   \
do                                                                \
{                                                                 \
    PERMUTE(S);                                                   \
    S[0] = XOR(S[0], COMB(vcreate_u64(HLEN), vcreate_u64(MLEN))); \
    PERMUTE(S);                                                   \
} while(0)

static void* (* const volatile burn)(void*, int, size_t) = memset;

typedef enum tag__
{
    ABS_TAG     = 0x00,
    ENC_TAG     = 0x01
} tag_t;

void crypto_aead_encrypt(
    unsigned char *c, size_t *clen,
    const unsigned char *h, size_t hlen,
    const unsigned char *m, size_t mlen,
    const unsigned char *nonce,
    const unsigned char *key
    )
{
    uint64x2_t S[8];

    /* absorption phase */
    INIT(S, key, nonce, WORDS(MRS_N)/2, ABS_TAG);
    ABSORB_DATA(S, h, hlen);
    ABSORB_DATA(S, m, mlen);
    FINALISE(S, hlen, mlen);

    /* extract tag */
    STORE(c + mlen, S[0]);
    STORE(c + mlen + BYTES(MRS_T)/2, S[1]);
    *clen = mlen + BYTES(MRS_T);

    /* encrypt message */
    INIT(S, key, c + mlen, WORDS(MRS_T)/2, ENC_TAG);
    ENCRYPT_DATA(S, c, m, mlen);
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
    uint64x2_t S[8];
    uint32x4_t T[2];

    if (clen < BYTES(MRS_T)) { return result; }

    /* decrypt message */
    INIT(S, key, c + clen - BYTES(MRS_T), WORDS(MRS_T)/2, ENC_TAG);
    DECRYPT_DATA(S, m, c, clen - BYTES(MRS_T));
    *mlen = clen - BYTES(MRS_T);

    /* absorb header and message */
    INIT(S, key, nonce, WORDS(MRS_N)/2, ABS_TAG);
    ABSORB_DATA(S, h, hlen);
    ABSORB_DATA(S, m, *mlen);
    FINALISE(S, hlen, *mlen);

    /* verify tag */
    T[0] = vceqq_u32( U64TOU32(S[0]), U8TOU32(vld1q_u8((uint8_t *)(c + clen - BYTES(MRS_T)  ))) );
    T[1] = vceqq_u32( U64TOU32(S[1]), U8TOU32(vld1q_u8((uint8_t *)(c + clen - BYTES(MRS_T)/2))) );
    T[0] = vandq_u32(T[0], T[1]);
    result = (0xFFFFFFFFFFFFFFFFULL == (vgetq_lane_u64(U32TOU64(T[0]), 0) & vgetq_lane_u64(U32TOU64(T[0]), 1)) ? 0 : -1);

    /* burn decrypted plaintext on authentication failure */
    if (result != 0) { burn(m, 0, *mlen); }

    return result;
}
