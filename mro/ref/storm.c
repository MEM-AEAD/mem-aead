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

/* rotation constants (BLAKE2) */
#define R0 32
#define R1 24
#define R2 16
#define R3 63

#define ROTR(x, c) ( ((x) >> (c)) | ((x) << (BITS(x) - (c))) )
#define ROTL(x, c) ( ((x) << (c)) | ((x) >> (BITS(x) - (c))) )

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
static STORM_INLINE void F(storm_word_t S[16])
{
    /* column step */
    G(S[ 0], S[ 4], S[ 8], S[12]);
    G(S[ 1], S[ 5], S[ 9], S[13]);
    G(S[ 2], S[ 6], S[10], S[14]);
    G(S[ 3], S[ 7], S[11], S[15]);
    /* diagonal step */
    G(S[ 0], S[ 5], S[10], S[15]);
    G(S[ 1], S[ 6], S[11], S[12]);
    G(S[ 2], S[ 7], S[ 8], S[13]);
    G(S[ 3], S[ 4], S[ 9], S[14]);
}

static STORM_INLINE void storm_permute(storm_state_t state, size_t R)
{
    storm_word_t * S = state->S;
    size_t i;
    for (i = 0; i < R; ++i)
    {
        F(S);
    }
}

static STORM_INLINE void storm_pad(uint8_t * out, const uint8_t * in, const size_t inlen)
{
    memset(out, 0, BYTES(STORM_B));
    memcpy(out, in, inlen);
    out[inlen] = 0x01;
}

static STORM_INLINE void storm_init_mask(storm_state_t mask, const unsigned char * k, const unsigned char * n)
{
    storm_word_t * L = mask->S;

    L[ 0] = LOAD(n + 0 * BYTES(STORM_W));
    L[ 1] = LOAD(n + 1 * BYTES(STORM_W));
    L[ 2] = 0;
    L[ 3] = 0;

    L[ 4] = 0;
    L[ 5] = 0;
    L[ 6] = 0;
    L[ 7] = 0;

    L[ 8] = 0;
    L[ 9] = 0;
    L[10] = STORM_L;
    L[11] = STORM_T;

    L[12] = LOAD(k + 0 * BYTES(STORM_W));
    L[13] = LOAD(k + 1 * BYTES(STORM_W));
    L[14] = LOAD(k + 2 * BYTES(STORM_W));
    L[15] = LOAD(k + 3 * BYTES(STORM_W));

    /* apply permutation */
    storm_permute(mask, STORM_L);

#if defined(STORM_DEBUG)
    printf("SETUP MASK:\n");
    print_state(mask);
#endif
}

/* phi */
static STORM_INLINE void storm_phi(storm_state_t mask)
{
    size_t i;
    storm_word_t * L = mask->S;
    storm_word_t t = ROTL(L[0], 53) ^ (L[5] << 13);
    for (i = 0; i < WORDS(STORM_B) - 1; ++i)
    {
        L[i] = L[i+1];
    }
    L[15] = t;
}

/* sigma: phi(x) ^ x */
static STORM_INLINE void storm_sigma(storm_state_t mask)
{
    size_t i;
    storm_word_t * L = mask->S;
    storm_word_t t = ROTL(L[0], 53) ^ (L[5] << 13);
    for (i = 0; i < WORDS(STORM_B) - 1; ++i)
    {
        L[i] ^= L[i+1];
    }
    L[15] ^= t;
}

/* lambda: phi^2(x) ^ phi(x) ^ x */
static STORM_INLINE void storm_lambda(storm_state_t mask)
{
    size_t i;
    storm_word_t * L = mask->S;
    storm_word_t t0 = ROTL(L[0], 53) ^ (L[5] << 13);
    storm_word_t t1 = ROTL(L[1], 53) ^ (L[6] << 13);
    for (i = 0; i < WORDS(STORM_B) - 2; ++i)
    {
        L[i] ^= L[i+1] ^ L[i+2];
    }
    L[14] ^= (L[15] ^ t0);
    L[15] ^= (t0 ^ t1);
}

static STORM_INLINE void storm_absorb_block(storm_state_t state, storm_state_t mask, const uint8_t * in)
{
    size_t i;
    storm_state_t block;
    storm_word_t * B = block->S;
    storm_word_t * S = state->S;
    storm_word_t * L = mask->S;

    /* load data and XOR mask */
    memset(block, 0, sizeof(storm_state_t));
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        B[i] = LOAD(in + i * BYTES(STORM_W)) ^ L[i];
    }

    /* apply permutation */
    storm_permute(block, STORM_L);

    /* XOR mask and absorb into S */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        S[i] ^=  B[i] ^ L[i];
    }

#if defined(STORM_DEBUG)
    printf("ABSORBING BLOCK\n");
    printf("IN:\n");
    print_bytes(in, BYTES(STORM_B));
    printf("\nSTATE:\n");
    print_state(state);
    printf("MASK:\n");
    print_state(mask);
#endif
}

static STORM_INLINE void storm_absorb_lastblock(storm_state_t state, storm_state_t mask, const uint8_t * in, size_t inlen)
{
    uint8_t block[BYTES(STORM_B)];
    storm_pad(block, in, inlen);
    storm_absorb_block(state, mask, block);
    burn(block, 0, BYTES(STORM_B));
}

static STORM_INLINE void storm_encrypt_block(storm_state_t mask, storm_state_t tag, size_t block_nr, uint8_t * out, const uint8_t * in)
{
    size_t i;
    storm_state_t block;
    storm_word_t * B = block->S;
    storm_word_t * L = mask->S;
    storm_word_t * T = tag->S;

    /* load mask and XOR authentication tag and block counter */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        B[i] = L[i];
    }
    B[0] ^= T[0];
    B[1] ^= T[1];
    B[2] ^= T[2];
    B[3] ^= T[3];
    B[15] ^= block_nr;

    /* apply permutation */
    storm_permute(block, STORM_L);

    /* encrypt block */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        B[i] ^= LOAD(in + i * BYTES(STORM_W)) ^ L[i];
        STORE(out + i * BYTES(STORM_W), B[i]);
    }

#if defined(STORM_DEBUG)
    printf("ENCRYPTING BLOCK\n");
    printf("IN:\n");
    print_bytes(in, BYTES(STORM_B));
    printf("OUT:\n");
    print_bytes(out, BYTES(STORM_B));
    printf("MASK:\n");
    print_state(mask);
#endif
}

static STORM_INLINE void storm_encrypt_lastblock(storm_state_t mask, storm_state_t tag, size_t block_nr, uint8_t * out, const uint8_t * in, size_t inlen)
{
    uint8_t block[BYTES(STORM_B)];
    memset(block, 0, BYTES(STORM_B));
    memcpy(block, in, inlen);
    storm_encrypt_block(mask, tag, block_nr, block, block);
    memcpy(out, block, inlen);
    burn(block, 0, BYTES(STORM_B));
}

/* low-level interface functions */
void storm_absorb_data(storm_state_t state, storm_state_t mask, const unsigned char * in, size_t inlen, tag_t flag)
{
    if (flag)
    {
        storm_sigma(mask);
    }
    while (inlen >= BYTES(STORM_B))
    {
        storm_absorb_block(state, mask, in);
        inlen -= BYTES(STORM_B);
        in    += BYTES(STORM_B);
        storm_phi(mask);
    }
    if (inlen > 0)
    {
        storm_absorb_lastblock(state, mask, in, inlen);
    }
}

void storm_encrypt_data(storm_state_t mask, storm_state_t tag, unsigned char * out, const unsigned char * in, size_t inlen)
{
    size_t i = 0;
    storm_lambda(mask);
    while (inlen >= BYTES(STORM_B))
    {
        storm_encrypt_block(mask, tag, i, out, in);
        inlen -= BYTES(STORM_B);
        in    += BYTES(STORM_B);
        out   += BYTES(STORM_B);
        i += 1;
    }
    if (inlen > 0)
    {
        storm_encrypt_lastblock(mask, tag, i, out, in, inlen);
    }
}

void storm_decrypt_data(storm_state_t mask, storm_state_t tag, unsigned char * out, const unsigned char * in, size_t inlen)
{
    storm_encrypt_data(mask, tag, out, in, inlen);
}

static STORM_INLINE void storm_finalise(storm_state_t state, storm_state_t mask, size_t hlen, size_t mlen)
{
    size_t i;
    storm_word_t * S = state->S;
    storm_word_t * L = mask->S;

    storm_sigma(mask);
    storm_sigma(mask);

    S[14] ^= hlen;
    S[15] ^= mlen;

    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        S[i] ^= L[i];
    }

    /* apply permutation */
    storm_permute(state, STORM_L);

    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        S[i] ^= L[i];
    }
}

void storm_output_tag(storm_state_t state, unsigned char * tag)
{
    size_t i = 0;
    storm_word_t * S = state->S;
    uint8_t block[BYTES(STORM_T)];

    for (i = 0; i < WORDS(STORM_T); ++i)
    {
        STORE(block + i * BYTES(STORM_W), S[i]);
    }
    memcpy(tag, block, BYTES(STORM_T));
    burn(block, 0, BYTES(STORM_T));

#if defined(STORM_DEBUG)
    printf("EXTRACTING TAG:\n");
    print_state(state);
    printf("TAG:\n");
    print_bytes(tag, BYTES(STORM_T));
#endif
}

int storm_verify_tag(const unsigned char * tag1, const unsigned char * tag2)
{
    unsigned acc = 0;
    size_t i;

    for (i = 0; i < BYTES(STORM_T); ++i)
    {
        acc |= tag1[i] ^ tag2[i];
    }
    return (((acc - 1) >> 8) & 1) - 1;
}


/* high level interface functions */
void storm_aead_encrypt(
    unsigned char *c, size_t *clen,
    const unsigned char *h, size_t hlen,
    const unsigned char *m, size_t mlen,
    const unsigned char *nonce,
    const unsigned char *key
    )
{
    storm_state_t state, la, le;

    memset(state, 0, sizeof(storm_state_t));
    storm_init_mask(le, key, nonce);

    /* absorb header and message */
    memcpy(la, le, sizeof(storm_state_t));
    storm_absorb_data(state, la, h, hlen, ABS_AD);

    memcpy(la, le, sizeof(storm_state_t));
    storm_absorb_data(state, la, m, mlen, ABS_MSG);

    memcpy(la, le, sizeof(storm_state_t));
    storm_finalise(state, la, hlen, mlen);

    /* extract tag */
    storm_output_tag(state, c + mlen);
    *clen = mlen + BYTES(STORM_T);

    /* encrypt message */
    storm_encrypt_data(le, state, c, m, mlen);

    /* empty buffers */
    burn(state, 0, sizeof(storm_state_t));
    burn(la, 0, sizeof(storm_state_t));
    burn(le, 0, sizeof(storm_state_t));
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
    unsigned char tag[BYTES(STORM_T)];
    storm_state_t state, la, le;
    storm_word_t * T = state->S;

    if (clen < BYTES(STORM_T)) { return -1; }

    storm_init_mask(le, key, nonce);
    memcpy(la, le, sizeof(storm_state_t));

    *mlen = clen - BYTES(STORM_T);

    /* set first 4 state words temporarily to received tag */
    T[ 0] = LOAD(c + *mlen + 0 * BYTES(STORM_W));
    T[ 1] = LOAD(c + *mlen + 1 * BYTES(STORM_W));
    T[ 2] = LOAD(c + *mlen + 2 * BYTES(STORM_W));
    T[ 3] = LOAD(c + *mlen + 3 * BYTES(STORM_W));

    /* decrypt message */
    storm_decrypt_data(le, state, m, c, clen - BYTES(STORM_T));

    /* reset state */
    memset(state, 0, sizeof(storm_state_t));

    /* absorb header and message */
    memcpy(le, la, sizeof(storm_state_t));
    storm_absorb_data(state, la, h, hlen, ABS_AD);

    memcpy(la, le, sizeof(storm_state_t));
    storm_absorb_data(state, la, m, *mlen, ABS_MSG);

    memcpy(la, le, sizeof(storm_state_t));
    storm_finalise(state, la, hlen, *mlen);

    /* extract tag */
    storm_output_tag(state, tag);

    /* verify tag */
    result = storm_verify_tag(c + clen - BYTES(STORM_T), tag);

    /* burn decrypted plaintext on authentication failure */
    if (result != 0) { burn(m, 0, *mlen); }

    /* empty buffers */
    burn(state, 0, sizeof(storm_state_t));
    burn(la, 0, sizeof(storm_state_t));
    burn(le, 0, sizeof(storm_state_t));

    return result;
}
