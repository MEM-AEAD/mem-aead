#include "mrs.h"

#if defined(MRS_DEBUG)
#include <stdio.h>
void print_state(uint64_t S[16]);
void print_bytes(const uint8_t * in, size_t inlen);
#endif

#define MRS_W 64              /* word size */
#define MRS_L 4               /* round number */
#define MRS_T (MRS_W *  4)    /* tag size */
#define MRS_N (MRS_W *  2)    /* nonce size */
#define MRS_K (MRS_W *  4)    /* key size */
#define MRS_B (MRS_W * 16)    /* permutation width */
#define MRS_C (MRS_W *  4)    /* capacity */
#define MRS_R (MRS_B - MRS_C) /* rate */

/* Rotation constants (BLAKE2) */
#define R0 32
#define R1 24
#define R2 16
#define R3 63

/* workaround for C89 compilers */
#if !defined(__cplusplus) && (!defined(__STDC_VERSION__) || __STDC_VERSION__ < 199901L)
  #if   defined(_MSC_VER)
    #define MRS_INLINE __inline
  #elif defined(__GNUC__)
    #define MRS_INLINE __inline__
  #else
    #define MRS_INLINE
  #endif
#else
  #define MRS_INLINE inline
#endif

#define BITS(x) (sizeof(x) * CHAR_BIT)
#define BYTES(x) (((x) + 7) / 8)
#define WORDS(x) (((x) + (MRS_W-1)) / MRS_W)

static MRS_INLINE uint64_t load64(const void * in)
{
#if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
    uint64_t v;
    memcpy(&v, in, sizeof v);
    return v;
#else
    const uint8_t * p = (const uint8_t *)in;
    return ((uint64_t)p[0] <<  0) |
           ((uint64_t)p[1] <<  8) |
           ((uint64_t)p[2] << 16) |
           ((uint64_t)p[3] << 24) |
           ((uint64_t)p[4] << 32) |
           ((uint64_t)p[5] << 40) |
           ((uint64_t)p[6] << 48) |
           ((uint64_t)p[7] << 56);
#endif
}

static MRS_INLINE void store64(void * out, const uint64_t v)
{
#if defined(__BYTE_ORDER__) && (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
    memcpy(out, &v, sizeof v);
#else
    uint8_t * p = (uint8_t *)out;
    p[0] = (uint8_t)(v >>  0);
    p[1] = (uint8_t)(v >>  8);
    p[2] = (uint8_t)(v >> 16);
    p[3] = (uint8_t)(v >> 24);
    p[4] = (uint8_t)(v >> 32);
    p[5] = (uint8_t)(v >> 40);
    p[6] = (uint8_t)(v >> 48);
    p[7] = (uint8_t)(v >> 56);
#endif
}

#define LOAD load64
#define STORE store64

static void* (* const volatile burn)(void*, int, size_t) = memset;

#define ROTR(x, c) ( ((x) >> (c)) | ((x) << (BITS(x) - (c))) )

/* quarter round */
#define G(A, B, C, D)                            \
do                                               \
{                                                \
    (A) += (B); (D) ^= (A); (D) = ROTR((D), R0); \
    (C) += (D); (B) ^= (C); (B) = ROTR((B), R1); \
    (A) += (B); (D) ^= (A); (D) = ROTR((D), R2); \
    (C) += (D); (B) ^= (C); (B) = ROTR((B), R3); \
} while (0)

/* double round */
static MRS_INLINE void F(mrs_word_t S[16])
{
    /* Column step */
    G(S[ 0], S[ 4], S[ 8], S[12]);
    G(S[ 1], S[ 5], S[ 9], S[13]);
    G(S[ 2], S[ 6], S[10], S[14]);
    G(S[ 3], S[ 7], S[11], S[15]);
    /* Diagonal step */
    G(S[ 0], S[ 5], S[10], S[15]);
    G(S[ 1], S[ 6], S[11], S[12]);
    G(S[ 2], S[ 7], S[ 8], S[13]);
    G(S[ 3], S[ 4], S[ 9], S[14]);
}

static MRS_INLINE void mrs_permute(mrs_state_t state)
{
    size_t i;
    mrs_word_t * S = state->S;
    for(i = 0; i < MRS_L; ++i)
    {
        F(S);
    }
}

static MRS_INLINE void mrs_pad(unsigned char * out, const size_t outlen, const unsigned char * in, const size_t inlen)
{
    memset(out, 0, outlen);
    memcpy(out, in, inlen);
    /*out[inlen] = 0x01;
    out[outlen - 1] |= 0x80;*/
}

static MRS_INLINE void mrs_absorb_block(mrs_state_t state, const unsigned char * in)
{
    size_t i;
    mrs_word_t * S = state->S;
    mrs_permute(state);
    for (i = 0; i < WORDS(MRS_B); ++i)
    {
        S[i] ^= LOAD(in + i * BYTES(MRS_W));
    }

#if defined(MRS_DEBUG)
    printf("ABSORBING BLOCK:\n");
    print_bytes(in, BYTES(MRS_B));
    printf("STATE:\n");
    print_state(state->S);
#endif
}

void mrs_absorb_lastblock(mrs_state_t state, const unsigned char * in, size_t inlen)
{
    uint8_t lastblock[BYTES(MRS_B)];
    mrs_pad(lastblock, sizeof lastblock, in, inlen);
    mrs_absorb_block(state, lastblock);
    burn(lastblock, 0, BYTES(MRS_B));
}

static MRS_INLINE void mrs_encrypt_block(mrs_state_t state, unsigned char * out, const unsigned char * in)
{
    size_t i;
    mrs_word_t * S = state->S;
    mrs_permute(state);
    for (i = 0; i < WORDS(MRS_R); ++i)
    {
        S[i] ^= LOAD(in + i * BYTES(MRS_W));
        STORE(out + i * BYTES(MRS_W), S[i]);
    }

#if defined(MRS_DEBUG)
    printf("ENCRYPTING BLOCK:\n");
    print_bytes(in, BYTES(MRS_R));
    printf("STATE:\n");
    print_state(state->S);
#endif
}

static MRS_INLINE void mrs_encrypt_lastblock(mrs_state_t state, unsigned char * out, const unsigned char * in, size_t inlen)
{
    uint8_t lastblock[BYTES(MRS_R)];
    mrs_pad(lastblock, sizeof lastblock, in, inlen);
    mrs_encrypt_block(state, lastblock, lastblock);
    memcpy(out, lastblock, inlen);
    burn(lastblock, 0, BYTES(MRS_R));
}

static MRS_INLINE void mrs_decrypt_block(mrs_state_t state, unsigned char * out, const unsigned char * in)
{
    size_t i;
    mrs_word_t * S = state->S;
    mrs_permute(state);
    for (i = 0; i < WORDS(MRS_R); ++i)
    {
        const mrs_word_t c = LOAD(in + i * BYTES(MRS_W));
        STORE(out + i * BYTES(MRS_W), S[i] ^ c);
        S[i] = c;
    }

#if defined(MRS_DEBUG)
    printf("DECRYPTING BLOCK:\n");
    print_bytes(in, BYTES(MRS_R));
    printf("STATE:\n");
    print_state(state->S);
#endif
}

static MRS_INLINE void mrs_decrypt_lastblock(mrs_state_t state, unsigned char * out, const unsigned char * in, size_t inlen)
{
    size_t i;
    mrs_word_t * S = state->S;
    uint8_t lastblock[BYTES(MRS_R)];
    mrs_permute(state);
    for(i = 0; i < WORDS(MRS_R); ++i)
    {
        STORE(lastblock + i * BYTES(MRS_W), S[i]);
    }

    /* undo padding */
    memcpy(lastblock, in, inlen);
    /*lastblock[inlen] ^= 0x01;
    lastblock[BYTES(MRS_R) - 1] ^= 0x80;*/

    for (i = 0; i < WORDS(MRS_R); ++i)
    {
        const mrs_word_t c = LOAD(lastblock + i * BYTES(MRS_W));
        STORE(lastblock + i * BYTES(MRS_W), S[i] ^ c);
        S[i] = c;
    }
    memcpy(out, lastblock, inlen);

#if defined(MRS_DEBUG)
    printf("DECRYPTING LASTBLOCK:\n");
    print_bytes(lastblock, BYTES(MRS_R));
    printf("STATE:\n");
    print_state(state->S);
#endif

    burn(lastblock, 0, sizeof lastblock);
}

/* low-level interface functions */
void mrs_init(mrs_state_t state, const unsigned char * k, const unsigned char * iv, const size_t ivlen, tag_t tag)
{
    size_t i;
    mrs_word_t * S = state->S;
    memset(state, 0, sizeof(mrs_state_t));
    for(i = 0; i < ivlen; ++i)
    {
        S[i] = LOAD(iv + i * BYTES(MRS_W));
    }
    S[ 9] = MRS_L;
    S[10] = MRS_T;
    S[11] = tag;

    S[12] = LOAD(k + 0 * BYTES(MRS_W));
    S[13] = LOAD(k + 1 * BYTES(MRS_W));
    S[14] = LOAD(k + 2 * BYTES(MRS_W));
    S[15] = LOAD(k + 3 * BYTES(MRS_W));

#if defined(MRS_DEBUG)
    printf("SETUP (%02X):\n", tag);
    print_state(state->S);
#endif
}

void mrs_absorb_data(mrs_state_t state, const unsigned char * in, size_t inlen)
{
    while (inlen >= BYTES(MRS_B))
    {
        mrs_absorb_block(state, in);
        inlen -= BYTES(MRS_B);
        in    += BYTES(MRS_B);
    }
    if(inlen > 0) {
        mrs_absorb_lastblock(state, in, inlen);
    }
}

void mrs_encrypt_data(mrs_state_t state, unsigned char * out, const unsigned char * in, size_t inlen)
{

    while (inlen >= BYTES(MRS_R))
    {
        mrs_encrypt_block(state, out, in);
        inlen -= BYTES(MRS_R);
        in    += BYTES(MRS_R);
        out   += BYTES(MRS_R);
    }
    if(inlen > 0) {
        mrs_encrypt_lastblock(state, out, in, inlen);
    }
}

void mrs_decrypt_data(mrs_state_t state, unsigned char * out, const unsigned char * in, size_t inlen)
{
    while (inlen >= BYTES(MRS_R))
    {
        mrs_decrypt_block(state, out, in);
        inlen -= BYTES(MRS_R);
        in    += BYTES(MRS_R);
        out   += BYTES(MRS_R);
    }
    if(inlen > 0) {
        mrs_decrypt_lastblock(state, out, in, inlen);
    }
}

void mrs_finalise(mrs_state_t state, size_t hlen, size_t mlen, unsigned char * tag)
{
    mrs_word_t * S = state->S;
    uint8_t lastblock[BYTES(MRS_R)];
    size_t i;

    /* finalise state */
    mrs_permute(state);

    S[0] ^= hlen;
    S[1] ^= mlen;

    mrs_permute(state);

    /* extract tag */
    for (i = 0; i < WORDS(MRS_R); ++i)
    {
        STORE(lastblock + i * BYTES(MRS_W), S[i]);
    }
    memcpy(tag, lastblock, BYTES(MRS_T));

#if defined(MRS_DEBUG)
    printf("FINALISED:\n");
    print_state(state->S);
#endif

    burn(lastblock, 0, BYTES(MRS_R));
}

int mrs_verify_tag(const unsigned char * tag1, const unsigned char * tag2)
{
    unsigned acc = 0;
    size_t i;

    for(i = 0; i < BYTES(MRS_T); ++i)
    {
        acc |= tag1[i] ^ tag2[i];
    }

    return (((acc - 1) >> 8) & 1) - 1;
}

/* high level interface functions */
void crypto_aead_encrypt(
    unsigned char *c, size_t *clen,
    const unsigned char *h, size_t hlen,
    const unsigned char *m, size_t mlen,
    const unsigned char *nonce,
    const unsigned char *key
    )
{
    mrs_state_t state;

    /* absorption phase */
    mrs_init(state, key, nonce, WORDS(MRS_N), ABS_TAG);
    mrs_absorb_data(state, h, hlen);
    mrs_absorb_data(state, m, mlen);
    mrs_finalise(state, hlen, mlen, c + mlen);
    *clen = mlen + BYTES(MRS_T);

    /* encryption phase */
    mrs_init(state, key, c + mlen, WORDS(MRS_T), ENC_TAG); /* re-initialise with key and authentication tag */
    mrs_encrypt_data(state, c, m, mlen);
    burn(state, 0, sizeof(mrs_state_t));
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
    unsigned char tag[BYTES(MRS_T)];
    mrs_state_t state;

    if (clen < BYTES(MRS_T)) { return -1; }

    /* decryption phase */
    mrs_init(state, key, c + clen - BYTES(MRS_T), WORDS(MRS_T), ENC_TAG); /* initialise with key and authentication tag */
    mrs_decrypt_data(state, m, c, clen - BYTES(MRS_T));
    *mlen = clen - BYTES(MRS_T);

    /* absorption phase */
    mrs_init(state, key, nonce, WORDS(MRS_N), ABS_TAG);
    mrs_absorb_data(state, h, hlen);
    mrs_absorb_data(state, m, *mlen);
    mrs_finalise(state, hlen, *mlen, tag);

    /* verification phase */
    result = mrs_verify_tag(c + clen - BYTES(MRS_T), tag);

    /* burn decrypted plaintext on authentication failure */
    if(result != 0) { burn(m, 0, *mlen); }

    burn(state, 0, sizeof(mrs_state_t));

    return result;
}
