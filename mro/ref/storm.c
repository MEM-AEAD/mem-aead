/*
   STORM reference source code package - reference C implementations

   Written in 2015 by Philipp Jovanovic <philipp@jovanovic.io>

   To the extent possible under law, the author(s) have dedicated all copyright
   and related and neighboring rights to this software to the public domain
   worldwide. This software is distributed without any warranty.

   You should have received a copy of the CC0 Public Domain Dedication along with
   this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "storm.h"
#include "storm_utils.h"

#if defined(STORM_DEBUG)
#include <inttypes.h>
#include <stdio.h>
#endif

#if STORM_W == 64

    #define LOAD load64
    #define STORE store64
    #define STORM_N (STORM_W *  2)   /* nonce size */
    #define STORM_K (STORM_W *  4)   /* key size */
    #define STORM_B (STORM_W * 16)   /* permutation width */

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

    #if defined(STORM_DEBUG)
        #define FMT "016" PRIX64
    #endif

#else
    #error "Invalid word size!"
#endif

#if defined(STORM_DEBUG)
static void print_state(storm_state_t state)
{
    static const char fmt[] = "%" FMT " "
                              "%" FMT " "
                              "%" FMT " "
                              "%" FMT "\n";
    const storm_word_t * S = state->S;
    printf(fmt, S[ 0],S[ 1],S[ 2],S[ 3]);
    printf(fmt, S[ 4],S[ 5],S[ 6],S[ 7]);
    printf(fmt, S[ 8],S[ 9],S[10],S[11]);
    printf(fmt, S[12],S[13],S[14],S[15]);
    printf("\n");
}

static void print_bytes(const uint8_t * in, size_t inlen)
{
    size_t i;
    for (i = 0; i < inlen; ++i)
    {
        printf("%02X ", in[i]);
        if (i%16 == 15)
        {
            printf("\n");
        }
    }
    printf("\n");
}
#endif

typedef enum tag__
{
    ABS_TAG     = 0x00,
    ENC_TAG     = 0x01
} tag_t;

static STORM_INLINE void storm_permute(storm_state_t state, size_t R)
{
    storm_word_t * S = state->S;
    size_t i;
    for(i = 0; i < R; ++i)
    {
        F(S);
    }
}

static STORM_INLINE void storm_pad(uint8_t * out, const uint8_t * in, const size_t inlen)
{
    memset(out, 0, BYTES(STORM_B));
    memcpy(out, in, inlen);
    out[inlen] = 0x01;
    out[BYTES(STORM_B) - 1] |= 0x80;
}

static STORM_INLINE void storm_init(storm_state_t k, const unsigned char * key, const unsigned char * iv, const size_t ivlen, tag_t tag)
{
    size_t i;
    storm_word_t * K = k->S;

    K[ 0] = 0;
    K[ 1] = 0;
    K[ 2] = 0;
    K[ 3] = 0;

    /* load nonce/tag */
    for(i = 0; i < ivlen; ++i)
    {
        K[i] = LOAD(iv + i * BYTES(STORM_W));
    }

    /* load key */
    K[ 4] = LOAD(key + 0 * BYTES(STORM_W));
    K[ 5] = LOAD(key + 1 * BYTES(STORM_W));
    K[ 6] = LOAD(key + 2 * BYTES(STORM_W));
    K[ 7] = LOAD(key + 3 * BYTES(STORM_W));

    K[ 8] = 0;
    K[ 9] = 0;
    K[10] = 0;
    K[11] = 0;

    /* inject parameters */
    K[12] = STORM_W;
    K[13] = STORM_R;
    K[14] = STORM_T;
    K[15] = tag;

    /* apply permutation */
    storm_permute(k, STORM_R);

#if defined(STORM_DEBUG)
    printf("SETUP KEY (%02X):\n", tag);
    print_state(k);
#endif

}

static STORM_INLINE void storm_update(storm_state_t k)
{
    size_t i;
    storm_word_t * K = k->S;
    storm_word_t t = ROTR(K[0], 11) ^ (K[5] << 13);
    for (i = 0; i < WORDS(STORM_B) - 1; ++i)
    {
        K[i] = K[i+1];
    }
    K[15] = t;
}

static STORM_INLINE void storm_absorb_block(storm_state_t state, storm_state_t k, const uint8_t * in)
{
    size_t i;
    storm_state_t block;
    storm_word_t * B = block->S;
    storm_word_t * S = state->S;
    storm_word_t * K = k->S; /* phi**{i}(K_a) */

    /* load data and XOR key */
    memset(block, 0, sizeof(storm_state_t));
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        B[i] = LOAD(in + i * BYTES(STORM_W)) ^ K[i];
    }

    /* apply permutation */
    storm_permute(block, STORM_R);

    /* XOR key and absorb into S */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        S[i] ^=  B[i] ^ K[i];
    }

    /* update key for the next block */
    storm_update(k);

#if defined(STORM_DEBUG)
    printf("ABSORBING BLOCK\n");
    printf("IN:\n");
    print_bytes(in, BYTES(STORM_B));
    printf("\nSTATE:\n");
    print_state(state);
#endif
}

static STORM_INLINE void storm_absorb_lastblock(storm_state_t state, storm_state_t k, const uint8_t * in, size_t inlen)
{
    uint8_t block[BYTES(STORM_B)];
    storm_pad(block, in, inlen);
    storm_absorb_block(state, k, block);
    burn(block, 0, BYTES(STORM_B));
}

static STORM_INLINE void storm_absorb_finalise(storm_state_t state, storm_state_t k, size_t hlen, size_t mlen)
{
    size_t i;
    storm_state_t block;
    storm_word_t * B = block->S;
    storm_word_t * S = state->S;
    storm_word_t * K = k->S; /* phi**{i}(K_a) */

    /* load key and XOR data lengths */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        B[i] = K[i];
    }
    B[14] ^= hlen;
    B[15] ^= mlen;

    /* apply permutation */
    storm_permute(block, STORM_R);

    /* XOR key and absorb into S */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        S[i] ^=  B[i] ^ K[i];
    }
}

static STORM_INLINE void storm_encrypt_block(const storm_state_t k, size_t block_nr, uint8_t * out, const uint8_t * in)
{
    size_t i;
    storm_state_t block;
    storm_word_t * B = block->S;
    const storm_word_t * K = k->S;

    /* load key and XOR block counter */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        B[i] = K[i];
    }
    B[15] ^= block_nr;

    /* apply permutation */
    storm_permute(block, STORM_R);

    /* encrypt block */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        B[i] ^= LOAD(in + i * BYTES(STORM_W)) ^ K[i];
        STORE(out + i * BYTES(STORM_W), B[i]);
    }

#if defined(STORM_DEBUG)
    printf("ENCRYPTING BLOCK\n");
    printf("IN:\n");
    print_bytes(in, BYTES(STORM_B));
    printf("OUT:\n");
    print_bytes(out, BYTES(STORM_B));
#endif
}

static STORM_INLINE void storm_encrypt_lastblock(const storm_state_t k, size_t block_nr, uint8_t * out, const uint8_t * in, size_t inlen)
{
    uint8_t block[BYTES(STORM_B)];
    memset(block, 0, BYTES(STORM_B));
    memcpy(block, in, inlen);
    storm_encrypt_block(k, block_nr, block, block);
    memcpy(out, block, inlen);
    burn(block, 0, BYTES(STORM_B));
}

/* low-level interface functions */
void storm_absorb_data(storm_state_t state, storm_state_t k, const unsigned char * in, size_t inlen)
{
    while (inlen >= BYTES(STORM_B))
    {
        storm_absorb_block(state, k, in);
        inlen -= BYTES(STORM_B);
        in    += BYTES(STORM_B);
    }
    if (inlen > 0)
    {
        storm_absorb_lastblock(state, k, in, inlen);
    }
}

void storm_encrypt_data(const storm_state_t k, unsigned char * out, const unsigned char * in, size_t inlen)
{
    size_t i = 0;
    while (inlen >= BYTES(STORM_B))
    {
        storm_encrypt_block(k, i, out, in);
        inlen -= BYTES(STORM_B);
        in    += BYTES(STORM_B);
        out   += BYTES(STORM_B);
        i += 1;
    }
    if (inlen > 0)
    {
        storm_encrypt_lastblock(k, i, out, in, inlen);
    }
}

void storm_decrypt_data(const storm_state_t k, unsigned char * out, const unsigned char * in, size_t inlen)
{
    storm_encrypt_data(k, out, in, inlen);
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
    storm_state_t state, k;

    /* absorb header and message */
    memset(state, 0, sizeof(storm_state_t));
    storm_init(k, key, nonce, WORDS(STORM_N), ABS_TAG);
    storm_absorb_data(state, k, h, hlen);
    storm_absorb_data(state, k, m, mlen);
    storm_absorb_finalise(state, k, hlen, mlen);

    /* extract tag */
    storm_output_tag(state, c + mlen);
    *clen = mlen + BYTES(STORM_T);

    /* encrypt message */
    storm_init(k, key, c + mlen, WORDS(STORM_T), ENC_TAG);
    storm_encrypt_data(k, c, m, mlen);

    /* empty buffers */
    burn(state, 0, sizeof(storm_state_t));
    burn(k, 0, sizeof(storm_state_t));
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
    storm_state_t state, k;

    if (clen < BYTES(STORM_T)) { return -1; }

    /* decrypt message */
    storm_init(k, key, c + clen - BYTES(STORM_T), WORDS(STORM_T), ENC_TAG);
    storm_decrypt_data(k, m, c, clen - BYTES(STORM_T));
    *mlen = clen - BYTES(STORM_T);

    /* absorb header and message */
    memset(state, 0, sizeof(storm_state_t));
    storm_init(k, key, nonce, WORDS(STORM_N), ABS_TAG);
    storm_absorb_data(state, k, h, hlen);
    storm_absorb_data(state, k, m, *mlen);
    storm_absorb_finalise(state, k, hlen, *mlen);

    /* extract tag */
    storm_output_tag(state, tag);

    /* verify tag */
    result = storm_verify_tag(c + clen - BYTES(STORM_T), tag);

    /* burn decrypted plaintext on authentication failure */
    if (result != 0) { burn(m, 0, *mlen); }

    /* empty buffers */
    burn(state, 0, sizeof(storm_state_t));
    burn(k, 0, sizeof(storm_state_t));

    return result;
}
