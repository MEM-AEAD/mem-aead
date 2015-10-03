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

/* Rotation constants (BLAKE2) */
#define R0 32
#define R1 24
#define R2 16
#define R3 63

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
static STORM_INLINE void F(storm_word_t S[16])
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

static STORM_INLINE void storm_permute(storm_state_t state)
{
    size_t i;
    storm_word_t * S = state->S;
    for(i = 0; i < STORM_L; ++i)
    {
        F(S);
    }
}

static STORM_INLINE void storm_pad(unsigned char * out, const size_t outlen, const unsigned char * in, const size_t inlen)
{
    memset(out, 0, outlen);
    memcpy(out, in, inlen);
    /*out[inlen] = 0x01;
    out[outlen - 1] |= 0x80;*/
}

static STORM_INLINE void storm_absorb_block(storm_state_t state, const unsigned char * in)
{
    size_t i;
    storm_word_t * S = state->S;
    storm_permute(state);
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        S[i] ^= LOAD(in + i * BYTES(STORM_W));
    }

#if defined(STORM_DEBUG)
    printf("ABSORBING BLOCK:\n");
    print_bytes(in, BYTES(STORM_B));
    printf("STATE:\n");
    print_state(state);
#endif
}

void storm_absorb_lastblock(storm_state_t state, const unsigned char * in, size_t inlen)
{
    uint8_t lastblock[BYTES(STORM_B)];
    storm_pad(lastblock, sizeof lastblock, in, inlen);
    storm_absorb_block(state, lastblock);
    burn(lastblock, 0, BYTES(STORM_B));
}

static STORM_INLINE void storm_encrypt_block(storm_state_t state, unsigned char * out, const unsigned char * in)
{
    size_t i;
    storm_word_t * S = state->S;
    storm_permute(state);
    for (i = 0; i < WORDS(STORM_R); ++i)
    {
        S[i] ^= LOAD(in + i * BYTES(STORM_W));
        STORE(out + i * BYTES(STORM_W), S[i]);
    }

#if defined(STORM_DEBUG)
    printf("ENCRYPTING BLOCK:\n");
    print_bytes(in, BYTES(STORM_R));
    printf("STATE:\n");
    print_state(state);
#endif
}

static STORM_INLINE void storm_encrypt_lastblock(storm_state_t state, unsigned char * out, const unsigned char * in, size_t inlen)
{
    uint8_t lastblock[BYTES(STORM_R)];
    storm_pad(lastblock, sizeof lastblock, in, inlen);
    storm_encrypt_block(state, lastblock, lastblock);
    memcpy(out, lastblock, inlen);
    burn(lastblock, 0, BYTES(STORM_R));
}

static STORM_INLINE void storm_decrypt_block(storm_state_t state, unsigned char * out, const unsigned char * in)
{
    size_t i;
    storm_word_t * S = state->S;
    storm_permute(state);
    for (i = 0; i < WORDS(STORM_R); ++i)
    {
        const storm_word_t c = LOAD(in + i * BYTES(STORM_W));
        STORE(out + i * BYTES(STORM_W), S[i] ^ c);
        S[i] = c;
    }

#if defined(STORM_DEBUG)
    printf("DECRYPTING BLOCK:\n");
    print_bytes(in, BYTES(STORM_R));
    printf("STATE:\n");
    print_state(state);
#endif
}

static STORM_INLINE void storm_decrypt_lastblock(storm_state_t state, unsigned char * out, const unsigned char * in, size_t inlen)
{
    size_t i;
    storm_word_t * S = state->S;
    uint8_t lastblock[BYTES(STORM_R)];
    storm_permute(state);
    for(i = 0; i < WORDS(STORM_R); ++i)
    {
        STORE(lastblock + i * BYTES(STORM_W), S[i]);
    }

    /* undo padding */
    memcpy(lastblock, in, inlen);
    /*lastblock[inlen] ^= 0x01;
    lastblock[BYTES(STORM_R) - 1] ^= 0x80;*/

    for (i = 0; i < WORDS(STORM_R); ++i)
    {
        const storm_word_t c = LOAD(lastblock + i * BYTES(STORM_W));
        STORE(lastblock + i * BYTES(STORM_W), S[i] ^ c);
        S[i] = c;
    }
    memcpy(out, lastblock, inlen);

#if defined(STORM_DEBUG)
    printf("DECRYPTING LASTBLOCK:\n");
    print_bytes(lastblock, BYTES(STORM_R));
    printf("STATE:\n");
    print_state(state);
#endif

    burn(lastblock, 0, sizeof lastblock);
}

/* low-level interface functions */
void storm_init(storm_state_t state, const unsigned char * k, const unsigned char * iv, const size_t ivlen, tag_t tag)
{
    size_t i;
    storm_word_t * S = state->S;
    memset(state, 0, sizeof(storm_state_t));
    for(i = 0; i < ivlen; ++i)
    {
        S[i] = LOAD(iv + i * BYTES(STORM_W));
    }
    S[ 9] = STORM_L;
    S[10] = STORM_T;
    S[11] = tag;

    S[12] = LOAD(k + 0 * BYTES(STORM_W));
    S[13] = LOAD(k + 1 * BYTES(STORM_W));
    S[14] = LOAD(k + 2 * BYTES(STORM_W));
    S[15] = LOAD(k + 3 * BYTES(STORM_W));

#if defined(STORM_DEBUG)
    printf("SETUP (%02X):\n", tag);
    print_state(state);
#endif
}

void storm_absorb_data(storm_state_t state, const unsigned char * in, size_t inlen)
{
    while (inlen >= BYTES(STORM_B))
    {
        storm_absorb_block(state, in);
        inlen -= BYTES(STORM_B);
        in    += BYTES(STORM_B);
    }
    if(inlen > 0) {
        storm_absorb_lastblock(state, in, inlen);
    }
}

void storm_encrypt_data(storm_state_t state, unsigned char * out, const unsigned char * in, size_t inlen)
{

    while (inlen >= BYTES(STORM_R))
    {
        storm_encrypt_block(state, out, in);
        inlen -= BYTES(STORM_R);
        in    += BYTES(STORM_R);
        out   += BYTES(STORM_R);
    }
    if(inlen > 0) {
        storm_encrypt_lastblock(state, out, in, inlen);
    }
}

void storm_decrypt_data(storm_state_t state, unsigned char * out, const unsigned char * in, size_t inlen)
{
    while (inlen >= BYTES(STORM_R))
    {
        storm_decrypt_block(state, out, in);
        inlen -= BYTES(STORM_R);
        in    += BYTES(STORM_R);
        out   += BYTES(STORM_R);
    }
    if(inlen > 0) {
        storm_decrypt_lastblock(state, out, in, inlen);
    }
}

void storm_finalise(storm_state_t state, size_t hlen, size_t mlen, unsigned char * tag)
{
    storm_word_t * S = state->S;
    uint8_t lastblock[BYTES(STORM_R)];
    size_t i;

    /* finalise state */
    storm_permute(state);

    S[0] ^= hlen;
    S[1] ^= mlen;

    storm_permute(state);

    /* extract tag */
    for (i = 0; i < WORDS(STORM_R); ++i)
    {
        STORE(lastblock + i * BYTES(STORM_W), S[i]);
    }
    memcpy(tag, lastblock, BYTES(STORM_T));

#if defined(STORM_DEBUG)
    printf("FINALISED:\n");
    print_state(state);
#endif

    burn(lastblock, 0, BYTES(STORM_R));
}

int storm_verify_tag(const unsigned char * tag1, const unsigned char * tag2)
{
    unsigned acc = 0;
    size_t i;

    for(i = 0; i < BYTES(STORM_T); ++i)
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
    storm_state_t state;

    /* absorption phase */
    storm_init(state, key, nonce, WORDS(STORM_N), ABS_TAG);
    storm_absorb_data(state, h, hlen);
    storm_absorb_data(state, m, mlen);
    storm_finalise(state, hlen, mlen, c + mlen);
    *clen = mlen + BYTES(STORM_T);

    /* encryption phase */
    storm_init(state, key, c + mlen, WORDS(STORM_T), ENC_TAG); /* re-initialise with key and authentication tag */
    storm_encrypt_data(state, c, m, mlen);
    burn(state, 0, sizeof(storm_state_t));
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
    storm_state_t state;

    if (clen < BYTES(STORM_T)) { return -1; }

    /* decryption phase */
    storm_init(state, key, c + clen - BYTES(STORM_T), WORDS(STORM_T), ENC_TAG); /* initialise with key and authentication tag */
    storm_decrypt_data(state, m, c, clen - BYTES(STORM_T));
    *mlen = clen - BYTES(STORM_T);

    /* absorption phase */
    storm_init(state, key, nonce, WORDS(STORM_N), ABS_TAG);
    storm_absorb_data(state, h, hlen);
    storm_absorb_data(state, m, *mlen);
    storm_finalise(state, hlen, *mlen, tag);

    /* verification phase */
    result = storm_verify_tag(c + clen - BYTES(STORM_T), tag);

    /* burn decrypted plaintext on authentication failure */
    if(result != 0) { burn(m, 0, *mlen); }

    burn(state, 0, sizeof(storm_state_t));

    return result;
}
