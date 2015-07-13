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

    /*
     * initialisation constants:
     * generated by (u0,...,u15) = F(F(0,...,15));
     */
    /*
    static const storm_word_t storm_u[16] =
    {
        0x901ABF1E4E0D2CA6ULL, 0x2A501C50B300F172ULL, 0x0F0D0CE7DD5462FBULL, 0xF96D0588D0A6E052ULL,
        0x53DFA665CCEB91E9ULL, 0xC0364181E8CB838FULL, 0x7A91E62A3D2FD2C2ULL, 0xDB55CC0737F96AB3ULL,
        0xD79CD8A45EA80C9CULL, 0xB52B68E611239851ULL, 0xCA749CD134276844ULL, 0x5544D09BCC165082ULL,
        0x8A50BFB3B7B174E7ULL, 0xFC3FA901CDFB5FBAULL, 0x2C128F8B5862FBC2ULL, 0xCD989101F81E5AEFULL
    };
    */

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

static STORM_INLINE void storm_permutation(storm_state_t state, size_t R)
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

static STORM_INLINE void storm_load_key(storm_state_t key, const unsigned char * k, const unsigned char * iv, const size_t ivlen, tag_t tag)
{
    size_t i;
    storm_word_t * S = key->S;

    S[ 0] = 0;
    S[ 1] = 0;
    S[ 2] = 0;
    S[ 3] = 0;

    /* load nonce/tag */
    for(i = 0; i < ivlen; ++i)
    {
        S[i] = LOAD(iv + i * BYTES(STORM_W));
    }

    /* load key */
    S[ 4] = LOAD(k + 0 * BYTES(STORM_W));
    S[ 5] = LOAD(k + 1 * BYTES(STORM_W));
    S[ 6] = LOAD(k + 2 * BYTES(STORM_W));
    S[ 7] = LOAD(k + 3 * BYTES(STORM_W));

    S[ 8] = 0;
    S[ 9] = 0;
    S[10] = 0;
    S[11] = 0;

    /* inject parameters */
    S[12] = STORM_W;
    S[13] = STORM_R;
    S[14] = STORM_T;
    S[15] = tag;

    /* apply permutation */
    storm_permutation(key, STORM_R);

#if defined(STORM_DEBUG)
    printf("SETUP KEY (%02X):\n", tag);
    print_state(key);
#endif

}

#if defined(M4)
static STORM_INLINE void storm_update_key(storm_state_t key)
{
    storm_word_t * S = key->S;
    storm_word_t t0 = ROTR(S[0], 9) ^ (S[ 9] >> 7);
    storm_word_t t1 = ROTR(S[1], 9) ^ (S[10] >> 7);
    storm_word_t t2 = ROTR(S[2], 9) ^ (S[11] >> 7);
    storm_word_t t3 = ROTR(S[3], 9) ^ (S[12] >> 7);

    S[ 0] = S[ 4];
    S[ 1] = S[ 5];
    S[ 2] = S[ 6];
    S[ 3] = S[ 7];

    S[ 4] = S[ 8];
    S[ 5] = S[ 9];
    S[ 6] = S[10];
    S[ 7] = S[11];

    S[ 8] = S[12];
    S[ 9] = S[13];
    S[10] = S[14];
    S[11] = S[15];

    S[12] = t0;
    S[13] = t1;
    S[14] = t2;
    S[15] = t3;
}
#elif defined(M256)
static STORM_INLINE void storm_update_key(storm_state_t key)
{
    storm_word_t * S = key->S;
    storm_word_t t0 = S[1];
    storm_word_t t1 = S[2];
    storm_word_t t2 = S[3];
    storm_word_t t3 = ROTR(S[0], 61) ^ (S[ 3] >> 5);

    S[ 0] = t0;
    S[ 1] = t1;
    S[ 2] = t2;
    S[ 3] = t3;

    S[ 4] = t0;
    S[ 5] = t1;
    S[ 6] = t2;
    S[ 7] = t3;

    S[ 8] = t0;
    S[ 9] = t1;
    S[10] = t2;
    S[11] = t3;

    S[12] = t0;
    S[13] = t1;
    S[14] = t2;
    S[15] = t3;
}
#else
static STORM_INLINE void storm_update_key(storm_state_t key)
{
    size_t i;
    storm_word_t * S = key->S;
    storm_word_t t = ROTR(S[0], 9) ^ (S[9] >> 7);
    for (i = 0; i < WORDS(STORM_B) - 1; ++i)
    {
        S[i] = S[i+1];
    }
    S[15] = t;
}
#endif

static STORM_INLINE void storm_absorb_block(storm_state_t state, storm_state_t kx, const uint8_t * in)
{
    size_t i;
    storm_state_t block;
    storm_word_t * BLK = block->S;
    storm_word_t * S = state->S; /* K_a */
    storm_word_t * KX = kx->S; /* phi**{i+1}(K_a) */

    /* update key */
    storm_update_key(kx);

    /* load data and XOR key */
    memset(block, 0, sizeof(storm_state_t));
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        BLK[i] = LOAD(in + i * BYTES(STORM_W)) ^ KX[i];
    }

    /* apply permutation */
    storm_permutation(block, STORM_R);

    /* absorb into K0 */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        S[i] ^=  BLK[i];
    }

#if defined(STORM_DEBUG)
    printf("ABSORBING BLOCK\n");
    printf("IN:\n");
    print_bytes(in, BYTES(STORM_B));
    printf("\nSTATE:\n");
    print_state(state);
#endif
}

static STORM_INLINE void storm_absorb_lastblock(storm_state_t state, storm_state_t kx, const uint8_t * in, size_t inlen)
{
    uint8_t block[BYTES(STORM_B)];
    storm_pad(block, in, inlen);
    storm_absorb_block(state, kx, block);
    burn(block, 0, BYTES(STORM_B));
}

static STORM_INLINE void storm_encrypt_block(const storm_state_t k, size_t block_nr, uint8_t * out, const uint8_t * in)
{
    size_t i;
    storm_state_t block;
    storm_word_t * BLK = block->S;
    const storm_word_t * K = k->S;

    /* Load key and XOR block counter */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        BLK[i] = K[i];
    }
    BLK[15] ^= block_nr;

    /* apply permutation */
    storm_permutation(block, STORM_R);

    /* encrypt block */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        BLK[i] ^= LOAD(in + i * BYTES(STORM_W)) ^ K[i];
        STORE(out + i * BYTES(STORM_W), BLK[i]);
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
void storm_absorb_data(storm_state_t state, storm_state_t kx, const unsigned char * in, size_t inlen)
{
    if(inlen > 0)
    {
        while (inlen >= BYTES(STORM_B))
        {
            storm_absorb_block(state, kx, in);
            inlen -= BYTES(STORM_B);
            in    += BYTES(STORM_B);
        }
        storm_absorb_lastblock(state, kx, in, inlen);
    }
}

void storm_encrypt_data(const storm_state_t k, unsigned char * out, const unsigned char * in, size_t inlen)
{
    if(inlen > 0)
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
        storm_encrypt_lastblock(k, i, out, in, inlen);
    }
}

void storm_decrypt_data(const storm_state_t k, unsigned char * out, const unsigned char * in, size_t inlen)
{
    storm_encrypt_data(k, out, in, inlen);
}

void storm_squeeze_tag(storm_state_t state, unsigned char * tag)
{
    size_t i = 0;
    storm_word_t * S = state->S;
    uint8_t block[BYTES(STORM_B - STORM_K)];

    storm_permutation(state, STORM_R);

    /* squeeze tag (TODO: do it sponge-like) */
    for (i = 0; i < WORDS(STORM_B - STORM_K); ++i)
    {
        STORE(block + i * BYTES(STORM_W), S[i]);
    }
    memcpy(tag, block, BYTES(STORM_T));
    burn(block, 0, BYTES(STORM_B - STORM_K));

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
    storm_word_t * S = state->S;

    /* absorb header and message */
#if defined(STORM_DEBUG)
    storm_load_key(state, key, nonce, WORDS(STORM_N), ABS_TAG); /* K_a */
    storm_update_key(state);
    print_state(state);
#endif
    storm_load_key(state, key, nonce, WORDS(STORM_N), ABS_TAG); /* K_a */
    memcpy(k, state, sizeof(storm_state_t)); /* phi**{i}(K_a) */
    storm_absorb_data(state, k, h, hlen);
    storm_absorb_data(state, k, m, mlen);

    /* add length encoding */
    S[14] ^= hlen;
    S[15] ^= mlen;

    /* extract tag */
    storm_squeeze_tag(state, c + mlen);
    *clen = mlen + BYTES(STORM_T);

    /* encrypt message */
    storm_load_key(k, key, c + mlen, WORDS(STORM_T), ENC_TAG); /* K_e */
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
    storm_word_t * S = state->S;

    if (clen < BYTES(STORM_T))
        return -1;

    /* decrypt message */
    storm_load_key(k, key, c + clen - BYTES(STORM_T), WORDS(STORM_T), ENC_TAG); /* K_e */
    storm_decrypt_data(k, m, c, clen - BYTES(STORM_T));
    *mlen = clen - BYTES(STORM_T);

    /* absorb header and message */
    storm_load_key(state, key, nonce, WORDS(STORM_N), ABS_TAG); /* K_a */
    memcpy(k, state, sizeof(storm_state_t)); /* phi**{i}(K_a) */
    storm_absorb_data(state, k, h, hlen);
    storm_absorb_data(state, k, m, *mlen);

    /* add length encoding */
    S[14] ^= hlen;
    S[15] ^= *mlen;

    /* extract tag */
    storm_squeeze_tag(state, tag);

    /* verify tag */
    result = storm_verify_tag(c + clen - BYTES(STORM_T), tag);

    /* burn decrypted plaintext on authentication failure */
    if(result != 0)
    {
        burn(m, 0, *mlen);
    }

    /* empty buffers */
    burn(state, 0, sizeof(storm_state_t));
    burn(k, 0, sizeof(storm_state_t));

    return result;
}
