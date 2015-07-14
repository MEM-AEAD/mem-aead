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

    /* inverse quarter round */
    #define GI(A, B, C, D)                           \
    do                                               \
    {                                                \
        (B) = ROTL((B), R3); (B) ^= (C); (C) -= (D); \
        (D) = ROTL((D), R2); (D) ^= (A); (A) -= (B); \
        (B) = ROTL((B), R1); (B) ^= (C); (C) -= (D); \
        (D) = ROTL((D), R0); (D) ^= (A); (A) -= (B); \
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

    /* inverse double round */
    static STORM_INLINE void FI(storm_word_t S[16])
    {
        /* Diagonal step */
        GI(S[ 0], S[ 5], S[10], S[15]);
        GI(S[ 1], S[ 6], S[11], S[12]);
        GI(S[ 2], S[ 7], S[ 8], S[13]);
        GI(S[ 3], S[ 4], S[ 9], S[14]);
        /* Column step */
        GI(S[ 0], S[ 4], S[ 8], S[12]);
        GI(S[ 1], S[ 5], S[ 9], S[13]);
        GI(S[ 2], S[ 6], S[10], S[14]);
        GI(S[ 3], S[ 7], S[11], S[15]);
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

static STORM_INLINE void storm_permute(storm_state_t state, size_t R)
{
    storm_word_t * S = state->S;
    size_t i;
    for(i = 0; i < R; ++i)
    {
        F(S);
    }
}

static STORM_INLINE void storm_permute_inverse(storm_state_t state, size_t R)
{
    storm_word_t * S = state->S;
    size_t i;
    for(i = 0; i < R; ++i)
    {
        FI(S);
    }
}

static STORM_INLINE void storm_pad(uint8_t * out, const uint8_t * in, const size_t inlen)
{
    memset(out, 0, BYTES(STORM_B));
    memcpy(out, in, inlen);
    out[inlen] = 0x01;
    out[BYTES(STORM_B) - 1] |= 0x80;
}

static STORM_INLINE void storm_load_key(storm_state_t key, const unsigned char * k, const unsigned char * iv, const size_t ivlen)
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
    S[15] = 0;

    /* apply permutation */
    storm_permute(key, STORM_R);

#if defined(STORM_DEBUG)
    printf("SETUP KEY:\n");
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

static STORM_INLINE void storm_absorb_block(storm_state_t state, storm_state_t x, const uint8_t * in)
{
    size_t i;
    storm_state_t block;
    storm_word_t * BLK = block->S;
    storm_word_t * S = state->S; /* X */
    storm_word_t * X = x->S; /* phi**{i+1}(X) */

    /* update key */
    storm_update_key(x);

    /* load data and XOR key */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        BLK[i] = LOAD(in + i * BYTES(STORM_W)) ^ X[i];
    }

    /* apply permutation */
    storm_permute(block, STORM_R);

    /* absorb into state */
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

static STORM_INLINE void storm_absorb_lastblock(storm_state_t state, storm_state_t x, const uint8_t * in, size_t inlen)
{
    uint8_t lastblock[BYTES(STORM_B)];
    storm_pad(lastblock, in, inlen);
    storm_absorb_block(state, x, lastblock);
    burn(lastblock, 0, BYTES(STORM_B));
}

static STORM_INLINE void storm_encrypt_block(storm_state_t state, storm_state_t x, uint8_t * out, const uint8_t * in)
{
    size_t i;
    storm_state_t block;
    storm_word_t * BLK = block->S;
    storm_word_t * S = state->S;
    storm_word_t * X = x->S;

    /* update key */
    storm_update_key(x);

    /* load message and XOR key */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        BLK[i] = LOAD(in + i * BYTES(STORM_W)) ^ X[i];
    }

    /* apply permutation */
    storm_permute(block, STORM_R);

    /* XOR key to block and store as ciphertext, XOR message to state */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        STORE(out + i * BYTES(STORM_W), BLK[i] ^ X[i]);
        S[i] ^= LOAD(in + i * BYTES(STORM_W));
    }

#if defined(STORM_DEBUG)
    printf("ENCRYPTING BLOCK\n");
    printf("IN:\n");
    print_bytes(in, BYTES(STORM_B));
    printf("OUT:\n");
    print_bytes(out, BYTES(STORM_B));
    printf("STATE:\n");
    print_state(state);
#endif
}

static STORM_INLINE void storm_encrypt_lastblock(storm_state_t state, storm_state_t x, uint8_t * out, const uint8_t * in, size_t inlen)
{
    if (inlen > 0)
    {
        size_t i;
        storm_state_t block;
        storm_word_t * BLK = block->S;
        storm_word_t * S = state->S;
        storm_word_t * X = x->S;
        uint8_t lastblock[BYTES(STORM_B)];

        /* update key */
        storm_update_key(x);

        /* load block with key */
        for (i = 0; i < WORDS(STORM_B); ++i)
        {
            BLK[i] = X[i];
        }
        BLK[15] ^= inlen;

        /* apply permutation */
        storm_permute(block, STORM_R);

        /* XOR padded message to state, XOR padded message and key to block, and extract ciphertext */
        storm_pad(lastblock, in, inlen);
        for (i = 0; i < WORDS(STORM_B); ++i)
        {
            S[i] ^= LOAD(lastblock + i * BYTES(STORM_W));
            STORE(lastblock + i * BYTES(STORM_W), BLK[i] ^ X[i] ^ LOAD(lastblock + i * BYTES(STORM_W)));
        }
        memcpy(out, lastblock, inlen);
        burn(lastblock, 0, BYTES(STORM_B));
#if defined(STORM_DEBUG)
        printf("ENCRYPTING LASTBLOCK\n");
        printf("IN:\n");
        print_bytes(in, inlen);
        printf("OUT:\n");
        print_bytes(out, inlen);
        printf("STATE:\n");
        print_state(state);
#endif
    }
}

static STORM_INLINE void storm_decrypt_block(storm_state_t state, storm_state_t x, uint8_t * out, const uint8_t * in)
{
    size_t i;
    storm_state_t block;
    storm_word_t * BLK = block->S;
    storm_word_t * S = state->S;
    storm_word_t * X = x->S;

    /* update key */
    storm_update_key(x);

    /* load ciphertext and XOR key */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        BLK[i] = LOAD(in + i * BYTES(STORM_W)) ^ X[i];
    }

    /* apply inverse permutation */
    storm_permute_inverse(block, STORM_R);

    /* XOR ciphertext to state, XOR key to block, and extract message */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        STORE(out + i * BYTES(STORM_W), BLK[i] ^ X[i]);
        S[i] ^= LOAD(out + i * BYTES(STORM_W));
    }

#if defined(STORM_DEBUG)
    printf("DECRYPTING BLOCK\n");
    printf("IN:\n");
    print_bytes(in, BYTES(STORM_B));
    printf("OUT:\n");
    print_bytes(out, BYTES(STORM_B));
    printf("STATE:\n");
    print_state(state);
#endif
}

static STORM_INLINE void storm_decrypt_lastblock(storm_state_t state, storm_state_t x, uint8_t * out, const uint8_t * in, size_t inlen)
{
    if (inlen > 0)
    {
        size_t i;
        storm_state_t block;
        storm_word_t * BLK = block->S;
        storm_word_t * S = state->S;
        storm_word_t * X = x->S;
        uint8_t lastblock[BYTES(STORM_B)];

        /* update key */
        storm_update_key(x);

        /* load block with key */
        for (i = 0; i < WORDS(STORM_B); ++i)
        {
            BLK[i] = X[i];
        }
        BLK[15] ^= inlen;

        /* apply permutation */
        storm_permute(block, STORM_R);

        /* XOR padded ciphertext and key to block, store message */
        storm_pad(lastblock, in, inlen);
        for (i = 0; i < WORDS(STORM_B); ++i)
        {
            STORE(lastblock + i * BYTES(STORM_W), BLK[i] ^ X[i] ^ LOAD(lastblock + i * BYTES(STORM_W)));
        }
        memcpy(out, lastblock, inlen);
        /* XOR message to state */
        storm_pad(lastblock, out, inlen);
        for (i = 0; i < WORDS(STORM_B); ++i) {
            S[i] ^= LOAD(lastblock + i * BYTES(STORM_W));
        }
        burn(lastblock, 0, BYTES(STORM_B));
#if defined(STORM_DEBUG)
        printf("DECRYPTING LASTBLOCK\n");
        printf("IN:\n");
        print_bytes(in, inlen);
        printf("OUT:\n");
        print_bytes(out, inlen);
        printf("STATE:\n");
        print_state(state);
#endif
    }
}

/* low-level interface functions */
void storm_absorb_data(storm_state_t state, storm_state_t x, const unsigned char * in, size_t inlen)
{
    if(inlen > 0)
    {
        while (inlen >= BYTES(STORM_B))
        {
            storm_absorb_block(state, x, in);
            inlen -= BYTES(STORM_B);
            in    += BYTES(STORM_B);
        }
        storm_absorb_lastblock(state, x, in, inlen);
    }
}

void storm_encrypt_data(storm_state_t state, storm_state_t x, unsigned char * out, const unsigned char * in, size_t inlen)
{
    if(inlen > 0)
    {
        while (inlen >= BYTES(STORM_B))
        {
            storm_encrypt_block(state, x, out, in);
            inlen -= BYTES(STORM_B);
            in    += BYTES(STORM_B);
            out   += BYTES(STORM_B);
        }
        storm_encrypt_lastblock(state, x, out, in, inlen);
    }
}

void storm_decrypt_data(storm_state_t state, storm_state_t x, unsigned char * out, const unsigned char * in, size_t inlen)
{
    if(inlen > 0)
    {
        while (inlen >= BYTES(STORM_B))
        {
            storm_decrypt_block(state, x, out, in);
            inlen -= BYTES(STORM_B);
            in    += BYTES(STORM_B);
            out   += BYTES(STORM_B);
        }
        storm_decrypt_lastblock(state, x, out, in, inlen);
    }
}

void storm_squeeze_tag(storm_state_t state, unsigned char * tag)
{
    size_t i = 0;
    storm_word_t * S = state->S;
    uint8_t block[BYTES(STORM_B - STORM_K)];

    storm_permute(state, STORM_R);

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
    storm_state_t state, x;

    /* compute X and init state */
    storm_load_key(x, key, nonce, WORDS(STORM_N));
    memcpy(state, x, sizeof(storm_state_t));

    /* absorb header */
    storm_absorb_data(state, x, h, hlen);

    /* encrypt message */
    storm_encrypt_data(state, x, c, m, mlen);
    *clen = mlen + BYTES(STORM_T);

    /* extract tag */
    storm_squeeze_tag(state, c + mlen);

    /* empty buffers */
    burn(state, 0, sizeof(storm_state_t));
    burn(x, 0, sizeof(storm_state_t));
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
    storm_state_t state, x;

    if (clen < BYTES(STORM_T))
        return -1;

    /* compute X and init state */
    storm_load_key(x, key, nonce, WORDS(STORM_N));
    memcpy(state, x, sizeof(storm_state_t));

    /* absorb header */
    storm_absorb_data(state, x, h, hlen);

    /* decrypt message */
    storm_decrypt_data(state, x, m, c, clen - BYTES(STORM_T));
    *mlen = clen - BYTES(STORM_T);

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
    burn(x, 0, sizeof(storm_state_t));

    return result;
}
