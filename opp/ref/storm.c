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

static STORM_INLINE void storm_init(storm_state_t key, const unsigned char * k, const unsigned char * n, tag_t tag)
{
    storm_word_t * K = key->S;

    /* load nonce */
    K[ 0] = LOAD(n + 0 * BYTES(STORM_W));
    K[ 1] = LOAD(n + 1 * BYTES(STORM_W));
    K[ 2] = 0;
    K[ 3] = 0;

    /* load key */
    K[ 4] = LOAD(k + 0 * BYTES(STORM_W));
    K[ 5] = LOAD(k + 1 * BYTES(STORM_W));
    K[ 6] = LOAD(k + 2 * BYTES(STORM_W));
    K[ 7] = LOAD(k + 3 * BYTES(STORM_W));

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
    storm_permute(key, STORM_R);

#if defined(STORM_DEBUG)
    printf("SETUP KEY:\n");
    print_state(key);
#endif

}

#if defined(M4)
static STORM_INLINE void storm_update(storm_state_t key)
{
    storm_word_t * K = key->S;
    storm_word_t t0 = ROTL(K[0], 53) ^ (K[ 5] << 13);
    storm_word_t t1 = ROTL(K[1], 53) ^ (K[ 6] << 13);
    storm_word_t t2 = ROTL(K[2], 53) ^ (K[ 7] << 13);
    storm_word_t t3 = ROTL(K[3], 53) ^ (K[ 8] << 13);

    K[ 0] = K[ 4];
    K[ 1] = K[ 5];
    K[ 2] = K[ 6];
    K[ 3] = K[ 7];

    K[ 4] = K[ 8];
    K[ 5] = K[ 9];
    K[ 6] = K[10];
    K[ 7] = K[11];

    K[ 8] = K[12];
    K[ 9] = K[13];
    K[10] = K[14];
    K[11] = K[15];

    K[12] = t0;
    K[13] = t1;
    K[14] = t2;
    K[15] = t3;
}
#elif defined(M256)
static STORM_INLINE void storm_update(storm_state_t key)
{
    storm_word_t * K = key->S;
    storm_word_t t0 = K[1];
    storm_word_t t1 = K[2];
    storm_word_t t2 = K[3];
    storm_word_t t3 = ROTR(K[0], 61) ^ (K[ 3] >> 5);

    K[ 0] = t0;
    K[ 1] = t1;
    K[ 2] = t2;
    K[ 3] = t3;

    K[ 4] = t0;
    K[ 5] = t1;
    K[ 6] = t2;
    K[ 7] = t3;

    K[ 8] = t0;
    K[ 9] = t1;
    K[10] = t2;
    K[11] = t3;

    K[12] = t0;
    K[13] = t1;
    K[14] = t2;
    K[15] = t3;
}
#else
static STORM_INLINE void storm_update(storm_state_t key)
{
    size_t i;
    storm_word_t * K = key->S;
    storm_word_t t = ROTL(K[0], 53) ^ (K[5] << 13);
    for (i = 0; i < WORDS(STORM_B) - 1; ++i)
    {
        K[i] = K[i+1];
    }
    K[15] = t;
}
#endif

static STORM_INLINE void storm_rotl256(storm_state_t key)
{
    storm_word_t t;
    storm_word_t * K = key->S;
    t = K[ 0]; K[ 0] = K[ 4]; K[ 4] = K[ 8]; K[ 8] = K[12]; K[12] = t;
    t = K[ 1]; K[ 1] = K[ 5]; K[ 5] = K[ 9]; K[ 9] = K[13]; K[13] = t;
    t = K[ 2]; K[ 2] = K[ 6]; K[ 6] = K[10]; K[10] = K[14]; K[14] = t;
    t = K[ 3]; K[ 3] = K[ 7]; K[ 7] = K[11]; K[11] = K[15]; K[15] = t;
}

static STORM_INLINE void storm_rotl512(storm_state_t key)
{
    storm_rotl256(key);
    storm_rotl256(key);
}

static STORM_INLINE void storm_rotl768(storm_state_t key)
{
    storm_word_t t;
    storm_word_t * K = key->S;
    t = K[12]; K[12] = K[ 8]; K[ 8] = K[ 4]; K[ 4] = K[ 0]; K[ 0] = t;
    t = K[13]; K[13] = K[ 9]; K[ 9] = K[ 5]; K[ 5] = K[ 1]; K[ 1] = t;
    t = K[14]; K[14] = K[10]; K[10] = K[ 6]; K[ 6] = K[ 2]; K[ 2] = t;
    t = K[15]; K[15] = K[11]; K[11] = K[ 7]; K[ 7] = K[ 3]; K[ 3] = t;
}

static STORM_INLINE void storm_absorb_block(storm_state_t state, storm_state_t k, const uint8_t * in)
{
    size_t i;
    storm_state_t block;
    storm_word_t * B = block->S;
    storm_word_t * S = state->S;
    storm_word_t * K = k->S;

    /* load data and XOR key */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        B[i] = LOAD(in + i * BYTES(STORM_W)) ^ K[i];
    }

    /* apply permutation */
    storm_permute(block, STORM_R);

    /* XOR key and absorb into state */
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
    uint8_t lastblock[BYTES(STORM_B)];
    storm_pad(lastblock, in, inlen);
    storm_rotl256(k); /* finalise key */
    storm_absorb_block(state, k, lastblock);
    burn(lastblock, 0, BYTES(STORM_B));
}

static STORM_INLINE void storm_encrypt_block(storm_state_t state, storm_state_t k, uint8_t * out, const uint8_t * in)
{
    size_t i;
    storm_state_t block;
    storm_word_t * B = block->S;
    storm_word_t * S = state->S;
    storm_word_t * K = k->S;

    /* load message and XOR key */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        B[i] = LOAD(in + i * BYTES(STORM_W)) ^ K[i];
    }

    /* apply permutation */
    storm_permute(block, STORM_R);

    /* XOR key to block and store as ciphertext, XOR message to state */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        STORE(out + i * BYTES(STORM_W), B[i] ^ K[i]);
        S[i] ^= LOAD(in + i * BYTES(STORM_W));
    }

    /* update key for next block */
    storm_update(k);

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

static STORM_INLINE void storm_encrypt_lastblock(storm_state_t state, storm_state_t k, uint8_t * out, const uint8_t * in, size_t inlen)
{
    size_t i;
    storm_state_t block;
    storm_word_t * B = block->S;
    storm_word_t * S = state->S;
    storm_word_t * K = k->S;
    uint8_t lastblock[BYTES(STORM_B)];

    /* finalise key */
    storm_rotl768(k);

    /* load block with key */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        B[i] = K[i];
    }

    /* apply permutation */
    storm_permute(block, STORM_R);

    /* XOR padded message to state, XOR padded message and key to block, and extract ciphertext */
    storm_pad(lastblock, in, inlen);
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        S[i] ^= LOAD(lastblock + i * BYTES(STORM_W));
        STORE(lastblock + i * BYTES(STORM_W), B[i] ^ K[i] ^ LOAD(lastblock + i * BYTES(STORM_W)));
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

static STORM_INLINE void storm_decrypt_block(storm_state_t state, storm_state_t k, uint8_t * out, const uint8_t * in)
{
    size_t i;
    storm_state_t block;
    storm_word_t * B = block->S;
    storm_word_t * S = state->S;
    storm_word_t * K = k->S;

    /* load ciphertext and XOR key */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        B[i] = LOAD(in + i * BYTES(STORM_W)) ^ K[i];
    }

    /* apply inverse permutation */
    storm_permute_inverse(block, STORM_R);

    /* XOR ciphertext to state, XOR key to block, and extract message */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        STORE(out + i * BYTES(STORM_W), B[i] ^ K[i]);
        S[i] ^= LOAD(out + i * BYTES(STORM_W));
    }

    /* update key for next block */
    storm_update(k);

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

static STORM_INLINE void storm_decrypt_lastblock(storm_state_t state, storm_state_t k, uint8_t * out, const uint8_t * in, size_t inlen)
{
    size_t i;
    storm_state_t block;
    storm_word_t * B = block->S;
    storm_word_t * S = state->S;
    storm_word_t * K = k->S;
    uint8_t lastblock[BYTES(STORM_B)];

    /* finalise key */
    storm_rotl768(k);

    /* load block with key */
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        B[i] = K[i];
    }

    /* apply permutation */
    storm_permute(block, STORM_R);

    /* XOR padded ciphertext and key to block, store message */
    storm_pad(lastblock, in, inlen);
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        STORE(lastblock + i * BYTES(STORM_W), B[i] ^ K[i] ^ LOAD(lastblock + i * BYTES(STORM_W)));
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

/* low-level interface functions */
void storm_absorb_data(storm_state_t state, storm_state_t k, const unsigned char * in, size_t inlen)
{
    while (inlen >= BYTES(STORM_B))
    {
        storm_absorb_block(state, k, in);
        inlen -= BYTES(STORM_B);
        in    += BYTES(STORM_B);
    }
    if(inlen > 0)
    {
        storm_absorb_lastblock(state, k, in, inlen);
    }
}

void storm_encrypt_data(storm_state_t state, storm_state_t k, unsigned char * out, const unsigned char * in, size_t inlen)
{
    while (inlen >= BYTES(STORM_B))
    {
        storm_encrypt_block(state, k, out, in);
        inlen -= BYTES(STORM_B);
        in    += BYTES(STORM_B);
        out   += BYTES(STORM_B);
    }
    if(inlen > 0)
    {
        storm_encrypt_lastblock(state, k, out, in, inlen);
    }
}

void storm_decrypt_data(storm_state_t state, storm_state_t k, unsigned char * out, const unsigned char * in, size_t inlen)
{
    while (inlen >= BYTES(STORM_B))
    {
        storm_decrypt_block(state, k, out, in);
        inlen -= BYTES(STORM_B);
        in    += BYTES(STORM_B);
        out   += BYTES(STORM_B);
    }
    if(inlen > 0)
    {
        storm_decrypt_lastblock(state, k, out, in, inlen);
    }
}

void storm_finalise(storm_state_t sa, storm_state_t se, storm_state_t k, unsigned char * tag)
{
    size_t i = 0;
    storm_word_t * SA = sa->S;
    storm_word_t * SE = se->S;
    storm_word_t * K = k->S;
    uint8_t block[BYTES(STORM_B)];

    /* finalise key */
    storm_rotl512(k);

    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        SE[i] ^= K[i];
    }

    storm_permute(se, STORM_R);

    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        SA[i] ^= SE[i] ^ K[i];
        STORE(block + i * BYTES(STORM_W), SA[i]);
    }

    memcpy(tag, block, BYTES(STORM_T));
    burn(block, 0, BYTES(STORM_B));
#if defined(STORM_DEBUG)
    printf("EXTRACTING TAG:\n");
    print_state(sa);
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
    storm_state_t sa, se, ka, ke;

    /* init states and keys */
    memset(sa, 0, sizeof(storm_state_t));
    memset(se, 0, sizeof(storm_state_t));
    storm_init(ka, key, nonce, ABS_TAG);
    storm_init(ke, key, nonce, ENC_TAG);

    /* absorb header */
    storm_absorb_data(sa, ka, h, hlen);

    /* encrypt message */
    storm_encrypt_data(se, ke, c, m, mlen);
    *clen = mlen + BYTES(STORM_T);

    /* finalise and extract tag */
    storm_finalise(sa, se, ka, c + mlen);

    /* empty buffers */
    burn(sa, 0, sizeof(storm_state_t));
    burn(se, 0, sizeof(storm_state_t));
    burn(ka, 0, sizeof(storm_state_t));
    burn(ke, 0, sizeof(storm_state_t));
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
    storm_state_t sa, se, ka, ke;

    if (clen < BYTES(STORM_T))
        return -1;

    /* init states and keys */
    memset(sa, 0, sizeof(storm_state_t));
    memset(se, 0, sizeof(storm_state_t));
    storm_init(ka, key, nonce, ABS_TAG);
    storm_init(ke, key, nonce, ENC_TAG);

    /* absorb header */
    storm_absorb_data(sa, ka, h, hlen);

    /* decrypt message */
    storm_decrypt_data(se, ke, m, c, clen - BYTES(STORM_T));
    *mlen = clen - BYTES(STORM_T);

    /* finalise and extract tag */
    storm_finalise(sa, se, ka, tag);

    /* verify tag */
    result = storm_verify_tag(c + clen - BYTES(STORM_T), tag);

    /* burn decrypted plaintext on authentication failure */
    if(result != 0)
    {
        burn(m, 0, *mlen);
    }

    /* empty buffers */
    burn(sa, 0, sizeof(storm_state_t));
    burn(se, 0, sizeof(storm_state_t));
    burn(ka, 0, sizeof(storm_state_t));
    burn(ke, 0, sizeof(storm_state_t));

    return result;
}
