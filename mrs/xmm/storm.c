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
#define STORM_C (STORM_W *  6)   /* capacity */
#define RATE (STORM_B - STORM_C) /* rate */

#define BYTES(x) (((x) + 7) / 8)

#if defined(_MSC_VER)
    #define ALIGN(x) __declspec(align(x))
#else
    #define ALIGN(x) __attribute__((aligned(x)))
#endif

#define XOR(A, B) _mm_xor_si128((A), (B))
#define AND(A, B) _mm_and_si128((A), (B))
#define ADD(A, B) _mm_add_epi64((A), (B))

#define U0 0x901ABF1E4E0D2CA6ULL
#define U1 0x2A501C50B300F172ULL
#define U2 0x0F0D0CE7DD5462FBULL
#define U3 0xF96D0588D0A6E052ULL
#define U4 0x53DFA665CCEB91E9ULL
#define U5 0xC0364181E8CB838FULL
#define U6 0x7A91E62A3D2FD2C2ULL
#define U7 0xDB55CC0737F96AB3ULL
#define U8 0xD79CD8A45EA80C9CULL
#define U9 0xB52B68E611239851ULL

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

#define G(A0, A1, B0, B1, C0, C1, D0, D1)   \
do                                          \
{                                           \
    A0 = ADD(A0, B0);    A1 = ADD(A1, B1);  \
    D0 = XOR(D0, A0);    D1 = XOR(D1, A1);  \
    D0 = ROT(D0, R0);    D1 = ROT(D1, R0);  \
                                            \
    C0 = ADD(C0, D0);    C1 = ADD(C1, D1);  \
    B0 = XOR(B0, C0);    B1 = XOR(B1, C1);  \
    B0 = ROT(B0, R1);    B1 = ROT(B1, R1);  \
                                            \
    A0 = ADD(A0, B0);    A1 = ADD(A1, B1);  \
    D0 = XOR(D0, A0);    D1 = XOR(D1, A1);  \
    D0 = ROT(D0, R2);    D1 = ROT(D1, R2);  \
                                            \
    C0 = ADD(C0, D0);    C1 = ADD(C1, D1);  \
    B0 = XOR(B0, C0);    B1 = XOR(B1, C1);  \
    B0 = ROT(B0, R3);    B1 = ROT(B1, R3);  \
} while(0)

#if defined(__SSSE3__)
#define DIAGONALIZE(A0, A1, B0, B1, C0, C1, D0, D1) \
do                                                  \
{                                                   \
    __m128i t0, t1;                                 \
                                                    \
    t0 = _mm_alignr_epi8(B1, B0, 8);                \
    t1 = _mm_alignr_epi8(B0, B1, 8);                \
    B0 = t0;                                        \
    B1 = t1;                                        \
                                                    \
    t0 = C0;                                        \
    C0 = C1;                                        \
    C1 = t0;                                        \
                                                    \
    t0 = _mm_alignr_epi8(D1, D0, 8);                \
    t1 = _mm_alignr_epi8(D0, D1, 8);                \
    D0 = t1;                                        \
    D1 = t0;                                        \
} while(0)

#define UNDIAGONALIZE(A0, A1, B0, B1, C0, C1, D0, D1) \
do                                                    \
{                                                     \
    __m128i t0, t1;                                   \
                                                      \
    t0 = _mm_alignr_epi8(B0, B1, 8);                  \
    t1 = _mm_alignr_epi8(B1, B0, 8);                  \
    B0 = t0;                                          \
    B1 = t1;                                          \
                                                      \
    t0 = C0;                                          \
    C0 = C1;                                          \
    C1 = t0;                                          \
                                                      \
    t0 = _mm_alignr_epi8(D0, D1, 8);                  \
    t1 = _mm_alignr_epi8(D1, D0, 8);                  \
    D0 = t1;                                          \
    D1 = t0;                                          \
} while(0)

#else

#define DIAGONALIZE(A0, A1, B0, B1, C0, C1, D0, D1)          \
do                                                           \
{                                                            \
    __m128i t0, t1;                                          \
                                                             \
    t0 = D0; t1 = B0;                                        \
    D0 = C0; C0 = C1; C1 = D0;                               \
    D0 = _mm_unpackhi_epi64(D1, _mm_unpacklo_epi64(t0, t0)); \
    D1 = _mm_unpackhi_epi64(t0, _mm_unpacklo_epi64(D1, D1)); \
    B0 = _mm_unpackhi_epi64(B0, _mm_unpacklo_epi64(B1, B1)); \
    B1 = _mm_unpackhi_epi64(B1, _mm_unpacklo_epi64(t1, t1)); \
} while(0)

#define UNDIAGONALIZE(A0, A1, B0, B1, C0, C1, D0, D1)        \
do                                                           \
{                                                            \
    __m128i t0, t1;                                          \
                                                             \
    t0 = C0; C0 = C1; C1 = t0;                               \
    t0 = B0; t1 = D0;                                        \
    B0 = _mm_unpackhi_epi64(B1, _mm_unpacklo_epi64(B0, B0)); \
    B1 = _mm_unpackhi_epi64(t0, _mm_unpacklo_epi64(B1, B1)); \
    D0 = _mm_unpackhi_epi64(D0, _mm_unpacklo_epi64(D1, D1)); \
    D1 = _mm_unpackhi_epi64(D1, _mm_unpacklo_epi64(t1, t1)); \
} while(0)

#endif

#define F(A0, A1, B0, B1, C0, C1, D0, D1)          \
do                                                 \
{                                                  \
    G(A0, A1, B0, B1, C0, C1, D0, D1);             \
    DIAGONALIZE(A0, A1, B0, B1, C0, C1, D0, D1);   \
    G(A0, A1, B0, B1, C0, C1, D0, D1);             \
    UNDIAGONALIZE(A0, A1, B0, B1, C0, C1, D0, D1); \
} while(0)

#define PERMUTE(A0, A1, B0, B1, C0, C1, D0, D1) \
do                                              \
{                                               \
    int i;                                      \
    for(i = 0; i < STORM_R; ++i)                \
    {                                           \
        F(A0, A1, B0, B1, C0, C1, D0, D1);      \
    }                                           \
} while(0)

#define PAD(BLOCK, BLOCKLEN, IN, INLEN) \
do                                      \
{                                       \
    memset(BLOCK, 0, BLOCKLEN);         \
    block_copy(BLOCK, IN, INLEN);       \
    BLOCK[INLEN] = 0x01;                \
    BLOCK[BLOCKLEN - 1] |= 0x80;        \
} while(0)

#define INJECT_DOMAIN_CONSTANT(A0, A1, B0, B1, C0, C1, D0, D1, TAG)       \
do                                                                        \
{                                                                         \
    D1 = XOR(D1, _mm_set_epi64x(TAG, 0));                                 \
} while(0)

#define ABSORB_BLOCK(A0, A1, B0, B1, C0, C1, D0, D1, IN, TAG)    \
do                                                               \
{                                                                \
    INJECT_DOMAIN_CONSTANT(A0, A1, B0, B1, C0, C1, D0, D1, TAG); \
    PERMUTE(A0, A1, B0, B1, C0, C1, D0, D1);                     \
    A0 = XOR(A0, LOADU(IN +  0));                                \
    A1 = XOR(A1, LOADU(IN + 16));                                \
    B0 = XOR(B0, LOADU(IN + 32));                                \
    B1 = XOR(B1, LOADU(IN + 48));                                \
    C0 = XOR(C0, LOADU(IN + 64));                                \
} while(0)

#define ABSORB_LASTBLOCK(A0, A1, B0, B1, C0, C1, D0, D1, IN, INLEN, TAG)    \
do                                                                          \
{                                                                           \
    ALIGN(64) unsigned char lastblock[BYTES(RATE)];                         \
    PAD(lastblock, sizeof lastblock, IN, INLEN);                            \
    ABSORB_BLOCK(A0, A1, B0, B1, C0, C1, D0, D1, lastblock, TAG);           \
} while(0)

#define ENCRYPT_BLOCK(A0, A1, B0, B1, C0, C1, D0, D1, OUT, IN)         \
do                                                                     \
{                                                                      \
    ABSORB_BLOCK(A0, A1, B0, B1, C0, C1, D0, D1, IN, ENC_PAYLOAD_TAG); \
    STOREU(OUT +  0, A0);                                              \
    STOREU(OUT + 16, A1);                                              \
    STOREU(OUT + 32, B0);                                              \
    STOREU(OUT + 48, B1);                                              \
    STOREU(OUT + 64, C0);                                              \
} while(0)

#define ENCRYPT_LASTBLOCK(A0, A1, B0, B1, C0, C1, D0, D1, OUT, IN, INLEN)   \
do                                                                          \
{                                                                           \
    ALIGN(64) unsigned char lastblock[BYTES(RATE)];                         \
    PAD(lastblock, sizeof lastblock, IN, INLEN);                            \
    ENCRYPT_BLOCK(A0, A1, B0, B1, C0, C1, D0, D1, lastblock, lastblock);    \
    block_copy(OUT, lastblock, INLEN);                                      \
} while(0)

#define DECRYPT_BLOCK(A0, A1, B0, B1, C0, C1, D0, D1, OUT, IN)               \
do                                                                           \
{                                                                            \
    __m128i W0, W1, W2, W3, W4;                                              \
    INJECT_DOMAIN_CONSTANT(A0, A1, B0, B1, C0, C1, D0, D1, ENC_PAYLOAD_TAG); \
    PERMUTE(A0, A1, B0, B1, C0, C1, D0, D1);                                 \
    W0 = LOADU(IN +  0); STOREU(OUT +  0, XOR(A0, W0)); A0 = W0;             \
    W1 = LOADU(IN + 16); STOREU(OUT + 16, XOR(A1, W1)); A1 = W1;             \
    W2 = LOADU(IN + 32); STOREU(OUT + 32, XOR(B0, W2)); B0 = W2;             \
    W3 = LOADU(IN + 48); STOREU(OUT + 48, XOR(B1, W3)); B1 = W3;             \
    W4 = LOADU(IN + 64); STOREU(OUT + 64, XOR(C0, W4)); C0 = W4;             \
} while(0)

#define DECRYPT_LASTBLOCK(A0, A1, B0, B1, C0, C1, D0, D1, OUT, IN, INLEN)      \
do                                                                             \
{                                                                              \
    ALIGN(64) unsigned char lastblock[BYTES(RATE)];                            \
    __m128i W0, W1, W2, W3, W4;                                                \
    INJECT_DOMAIN_CONSTANT(A0, A1, B0, B1, C0, C1, D0, D1, ENC_PAYLOAD_TAG);   \
    PERMUTE(A0, A1, B0, B1, C0, C1, D0, D1);                                   \
    STOREU(lastblock +   0, A0);                                               \
    STOREU(lastblock +  16, A1);                                               \
    STOREU(lastblock +  32, B0);                                               \
    STOREU(lastblock +  48, B1);                                               \
    STOREU(lastblock +  64, C0);                                               \
    block_copy(lastblock, IN, INLEN);                                          \
    lastblock[INLEN] ^= 0x01;                                                  \
    lastblock[BYTES(RATE)-1] ^= 0x80;                                          \
    W0 = LOADU(lastblock +  0); STOREU(lastblock +  0, XOR(A0, W0)); A0 = W0;  \
    W1 = LOADU(lastblock + 16); STOREU(lastblock + 16, XOR(A1, W1)); A1 = W1;  \
    W2 = LOADU(lastblock + 32); STOREU(lastblock + 32, XOR(B0, W2)); B0 = W2;  \
    W3 = LOADU(lastblock + 48); STOREU(lastblock + 48, XOR(B1, W3)); B1 = W3;  \
    W4 = LOADU(lastblock + 64); STOREU(lastblock + 64, XOR(C0, W4)); C0 = W4;  \
    block_copy(OUT, lastblock, INLEN);                                         \
} while(0)

#define INITIALISE(A0, A1, B0, B1, C0, C1, D0, D1, N, K0, K1, TAG)          \
do                                                                          \
{                                                                           \
    A0 = _mm_set_epi64x(U1, U0);                                            \
    A1 = K0;                                                                \
    B0 = _mm_set_epi64x(U2, N[0]);                                          \
    B1 = _mm_set_epi64x(N[1], U3);                                          \
    C0 = K1;                                                                \
    C1 = _mm_set_epi64x(U5, U4);                                            \
    D0 = _mm_set_epi64x(U7, U6);                                            \
    D1 = _mm_set_epi64x(U9, U8);                                            \
    A0 = XOR(A0, _mm_set_epi64x(0, STORM_W));                               \
    B0 = XOR(B0, _mm_set_epi64x(STORM_R, 0));                               \
    C1 = XOR(C1, _mm_set_epi64x(0, STORM_T));                               \
    D1 = XOR(D1, _mm_set_epi64x(TAG, 0));                                   \
    PERMUTE(A0, A1, B0, B1, C0, C1, D0, D1);                                \
} while(0)

#define ABSORB_DATA(A0, A1, B0, B1, C0, C1, D0, D1, IN, INLEN, TAG)         \
do                                                                          \
{                                                                           \
    if(INLEN > 0)                                                           \
    {                                                                       \
        size_t i = 0;                                                       \
        size_t l = INLEN;                                                   \
        while(l >= BYTES(RATE))                                             \
        {                                                                   \
            ABSORB_BLOCK(A0, A1, B0, B1, C0, C1, D0, D1, IN + i, TAG);      \
            i += BYTES(RATE); l -= BYTES(RATE);                             \
        }                                                                   \
        ABSORB_LASTBLOCK(A0, A1, B0, B1, C0, C1, D0, D1, IN + i, l, TAG);   \
    }                                                                       \
} while(0)

#define ENCRYPT_DATA(A0, A1, B0, B1, C0, C1, D0, D1, OUT, IN, INLEN)            \
do                                                                              \
{                                                                               \
    if(INLEN > 0)                                                               \
    {                                                                           \
        size_t i = 0;                                                           \
        size_t l = INLEN;                                                       \
        while(l >= BYTES(RATE))                                                 \
        {                                                                       \
            ENCRYPT_BLOCK(A0, A1, B0, B1, C0, C1, D0, D1, OUT + i, IN + i);     \
            i += BYTES(RATE); l -= BYTES(RATE);                                 \
        }                                                                       \
        ENCRYPT_LASTBLOCK(A0, A1, B0, B1, C0, C1, D0, D1, OUT + i, IN + i, l);  \
    }                                                                           \
} while(0)

#define DECRYPT_DATA(A0, A1, B0, B1, C0, C1, D0, D1, OUT, IN, INLEN)            \
do                                                                              \
{                                                                               \
    if(INLEN > 0)                                                               \
    {                                                                           \
        size_t i = 0;                                                           \
        size_t l = INLEN;                                                       \
        while(l >= BYTES(RATE))                                                 \
        {                                                                       \
            DECRYPT_BLOCK(A0, A1, B0, B1, C0, C1, D0, D1, OUT + i, IN + i);     \
            i += BYTES(RATE), l -= BYTES(RATE);                                 \
        }                                                                       \
        DECRYPT_LASTBLOCK(A0, A1, B0, B1, C0, C1, D0, D1, OUT + i, IN + i, l);  \
    }                                                                           \
} while(0)

#define FINALISE(A0, A1, B0, B1, C0, C1, D0, D1)                            \
do                                                                          \
{                                                                           \
    INJECT_DOMAIN_CONSTANT(A0, A1, B0, B1, C0, C1, D0, D1, ABS_FINAL_TAG);  \
    PERMUTE(A0, A1, B0, B1, C0, C1, D0, D1);                                \
    PERMUTE(A0, A1, B0, B1, C0, C1, D0, D1);                                \
} while(0)

/* inlen <= 80 */
static void block_copy(unsigned char *out, const unsigned char *in, const size_t inlen)
{
    size_t i = 0;
    if( inlen & 64 )
    {
        STOREU(out + i +  0, LOADU(in + i +  0));
        STOREU(out + i + 16, LOADU(in + i + 16));
        STOREU(out + i + 32, LOADU(in + i + 32));
        STOREU(out + i + 48, LOADU(in + i + 48));
        i += 64;
    }
    if( inlen & 32 )
    {
        STOREU(out + i +  0, LOADU(in + i +  0));
        STOREU(out + i + 16, LOADU(in + i + 16));
        i += 32;
    }
    if( inlen & 16 )
    {
        STOREU(out + i +  0, LOADU(in + i +  0));
        i += 16;
    }
    if( inlen & 8 )
    {
        memcpy(out + i, in + i, 8);
        i += 8;
    }
    if( inlen & 4 )
    {
        memcpy(out + i, in + i, 4);
        i += 4;
    }
    if( inlen & 2 )
    {
        memcpy(out + i, in + i, 2);
        i += 2;
    }
    if( inlen & 1 )
    {
        memcpy(out + i, in + i, 1);
        i += 1;
    }
}

typedef enum tag__
{
    ABS_INIT_TAG     = 0x00,
    ABS_HEADER_TAG   = 0x01,
    ABS_PAYLOAD_TAG  = 0x02,
    ABS_TRAILER_TAG  = 0x04,
    ABS_FINAL_TAG    = 0x08,
    ENC_INIT_TAG     = 0x10,
    ENC_PAYLOAD_TAG  = 0x12
} tag_t;

void storm_aead_encrypt(
    unsigned char *c, size_t *clen,
    const unsigned char *h, size_t hlen,
    const unsigned char *p, size_t mlen,
    const unsigned char *nonce,
    const unsigned char *key
    )
{
    __m128i A0, A1, B0, B1, C0, C1, D0, D1;
    const __m128i N  = LOADU(nonce);
    const __m128i K0 = LOADU(key +  0);
    const __m128i K1 = LOADU(key + 16);

    *clen = mlen + BYTES(STORM_T);

    /* absorption phase */
    INITIALISE(A0, A1, B0, B1, C0, C1, D0, D1, N, K0, K1, ABS_INIT_TAG);
    ABSORB_DATA(A0, A1, B0, B1, C0, C1, D0, D1, h, hlen, ABS_HEADER_TAG);
    ABSORB_DATA(A0, A1, B0, B1, C0, C1, D0, D1, p, mlen, ABS_PAYLOAD_TAG);
    FINALISE(A0, A1, B0, B1, C0, C1, D0, D1);
    STOREU(c + mlen +                0, A0);
    STOREU(c + mlen + BYTES(STORM_T)/2, A1);

    /* encryption phase */
    INITIALISE(A0, A1, B0, B1, C0, C1, D0, D1, LOADU(c + mlen), K0, K1, ENC_INIT_TAG);
    ABSORB_DATA(A0, A1, B0, B1, C0, C1, D0, D1, c + mlen + BYTES(STORM_T)/2, BYTES(STORM_T)/2, ENC_INIT_TAG);
    ENCRYPT_DATA(A0, A1, B0, B1, C0, C1, D0, D1, c, p, mlen);
}

int storm_aead_decrypt(
    unsigned char *m, size_t *mlen,
    const unsigned char *h, size_t hlen,
    const unsigned char *c, size_t clen,
    const unsigned char *nonce,
    const unsigned char *key
    )
{
    __m128i A0, A1, B0, B1, C0, C1, D0, D1;
    const __m128i N  = LOADU(nonce);
    const __m128i K0 = LOADU(key +  0);
    const __m128i K1 = LOADU(key + 16);

    if (clen < BYTES(STORM_T))
        return -1;

    /* decryption phase */
    INITIALISE(A0, A1, B0, B1, C0, C1, D0, D1, LOADU(c + clen - BYTES(STORM_T)), K0, K1, ENC_INIT_TAG);
    ABSORB_DATA(A0, A1, B0, B1, C0, C1, D0, D1, c + clen - BYTES(STORM_T)/2, BYTES(STORM_T)/2, ENC_INIT_TAG);
    DECRYPT_DATA(A0, A1, B0, B1, C0, C1, D0, D1, m, c, clen - BYTES(STORM_T));
    *mlen = clen - BYTES(STORM_T);

    /* absorption phase */
    INITIALISE(A0, A1, B0, B1, C0, C1, D0, D1, N, K0, K1, ABS_INIT_TAG);
    ABSORB_DATA(A0, A1, B0, B1, C0, C1, D0, D1, h, hlen, ABS_HEADER_TAG);
    ABSORB_DATA(A0, A1, B0, B1, C0, C1, D0, D1, m, *mlen, ABS_PAYLOAD_TAG);
    FINALISE(A0, A1, B0, B1, C0, C1, D0, D1);

    /* verify tag */
    A0 = _mm_cmpeq_epi8(A0, LOADU(c + clen - BYTES(STORM_T)  ));
    A1 = _mm_cmpeq_epi8(A1, LOADU(c + clen - BYTES(STORM_T)/2));
    return (((unsigned long)_mm_movemask_epi8(AND(A0, A1)) + 1) >> 16) - 1;
}
