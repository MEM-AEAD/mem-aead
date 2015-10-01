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

/* inverse double round */
static STORM_INLINE void FI(storm_word_t S[16])
{
    /* diagonal step */
    GI(S[ 0], S[ 5], S[10], S[15]);
    GI(S[ 1], S[ 6], S[11], S[12]);
    GI(S[ 2], S[ 7], S[ 8], S[13]);
    GI(S[ 3], S[ 4], S[ 9], S[14]);
    /* column step */
    GI(S[ 0], S[ 4], S[ 8], S[12]);
    GI(S[ 1], S[ 5], S[ 9], S[13]);
    GI(S[ 2], S[ 6], S[10], S[14]);
    GI(S[ 3], S[ 7], S[11], S[15]);
}

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
    L[10] = STORM_L * 2; /* we count single not double rounds */
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
    const size_t n = WORDS(STORM_B);
    storm_state_t block;
    storm_word_t * B = block->S;
    storm_word_t * S = state->S;
    storm_word_t * L = mask->S;

    /* load data and XOR mask */
    for (i = 0; i < n; ++i)
    {
        B[i] = LOAD(in + i * BYTES(STORM_W)) ^ L[i];
    }

    /* apply permutation */
    storm_permute(block, STORM_L);

    /* XOR mask and absorb into state */
    for (i = 0; i < n; ++i)
    {
        S[i] ^= B[i] ^ L[i];
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
    size_t i;
    const size_t n = WORDS(STORM_B);
    storm_state_t block;
    storm_word_t * B = block->S;
    storm_word_t * S = state->S;
    storm_word_t * L = mask->S;
    uint8_t lastblock[BYTES(STORM_B)];

    storm_pad(lastblock, in, inlen);

    /* load data and XOR mask */
    for (i = 0; i < n; ++i)
    {
        B[i] = LOAD(lastblock + i * BYTES(STORM_W)) ^ L[i];
    }

    /* apply permutation */
    storm_permute(block, STORM_L);

    /* XOR mask and absorb into state */
    for (i = 0; i < n; ++i)
    {
        S[i] ^= B[i] ^ L[i];
    }

#if defined(STORM_DEBUG)
    printf("ABSORBING LASTBLOCK\n");
    printf("IN:\n");
    print_bytes(in, inlen);
    printf("\nSTATE:\n");
    print_state(state);
    printf("MASK:\n");
    print_state(mask);
#endif
    burn(lastblock, 0, BYTES(STORM_B));
}

static STORM_INLINE void storm_encrypt_block(storm_state_t state, storm_state_t mask, uint8_t * out, const uint8_t * in)
{
    size_t i;
    const size_t n = WORDS(STORM_B);
    storm_state_t block;
    storm_word_t * B = block->S;
    storm_word_t * S = state->S;
    storm_word_t * L = mask->S;

    /* load message and XOR mask */
    for (i = 0; i < n; ++i)
    {
        B[i] = LOAD(in + i * BYTES(STORM_W)) ^ L[i];
    }

    /* apply permutation */
    storm_permute(block, STORM_L);

    /* XOR mask to block and store as ciphertext, XOR message to state */
    for (i = 0; i < n; ++i)
    {
        STORE(out + i * BYTES(STORM_W), B[i] ^ L[i]);
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
    printf("MASK:\n");
    print_state(mask);
#endif
}

static STORM_INLINE void storm_encrypt_lastblock(storm_state_t state, storm_state_t mask, uint8_t * out, const uint8_t * in, size_t inlen)
{
    size_t i;
    const size_t n = WORDS(STORM_B);
    storm_state_t block;
    storm_word_t * B = block->S;
    storm_word_t * S = state->S;
    storm_word_t * L = mask->S;
    uint8_t lastblock[BYTES(STORM_B)];

    /* load block with mask */
    for (i = 0; i < n; ++i)
    {
        B[i] = L[i];
    }

    /* apply permutation */
    storm_permute(block, STORM_L);

    /* XOR padded message to state, XOR padded message and mask to block, and extract ciphertext */
    storm_pad(lastblock, in, inlen);
    for (i = 0; i < WORDS(STORM_B); ++i)
    {
        S[i] ^= LOAD(lastblock + i * BYTES(STORM_W));
        STORE(lastblock + i * BYTES(STORM_W), B[i] ^ L[i] ^ LOAD(lastblock + i * BYTES(STORM_W)));
    }
    memcpy(out, lastblock, inlen);

#if defined(STORM_DEBUG)
    printf("ENCRYPTING LASTBLOCK\n");
    printf("IN:\n");
    print_bytes(in, inlen);
    printf("OUT:\n");
    print_bytes(out, inlen);
    printf("STATE:\n");
    print_state(state);
    printf("MASK:\n");
    print_state(mask);
#endif
    burn(lastblock, 0, BYTES(STORM_B));
}

static STORM_INLINE void storm_decrypt_block(storm_state_t state, storm_state_t mask, uint8_t * out, const uint8_t * in)
{
    size_t i;
    const size_t n = WORDS(STORM_B);
    storm_state_t block;
    storm_word_t * B = block->S;
    storm_word_t * S = state->S;
    storm_word_t * L = mask->S;

    /* load ciphertext and XOR mask */
    for (i = 0; i < n; ++i)
    {
        B[i] = LOAD(in + i * BYTES(STORM_W)) ^ L[i];
    }

    /* apply inverse permutation */
    storm_permute_inverse(block, STORM_L);

    /* XOR ciphertext to state, XOR mask to block, and extract message */
    for (i = 0; i < n; ++i)
    {
        STORE(out + i * BYTES(STORM_W), B[i] ^ L[i]);
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
    printf("MASK:\n");
    print_state(mask);
#endif
}

static STORM_INLINE void storm_decrypt_lastblock(storm_state_t state, storm_state_t mask, uint8_t * out, const uint8_t * in, size_t inlen)
{
    size_t i;
    const size_t n = WORDS(STORM_B);
    storm_state_t block;
    storm_word_t * B = block->S;
    storm_word_t * S = state->S;
    storm_word_t * L = mask->S;
    uint8_t lastblock[BYTES(STORM_B)];

    /* load block with key */
    for (i = 0; i < n; ++i)
    {
        B[i] = L[i];
    }

    /* apply permutation */
    storm_permute(block, STORM_L);

    /* XOR padded ciphertext and key to block, store message */
    storm_pad(lastblock, in, inlen);
    for (i = 0; i < n; ++i)
    {
        STORE(lastblock + i * BYTES(STORM_W), B[i] ^ L[i] ^ LOAD(lastblock + i * BYTES(STORM_W)));
    }
    memcpy(out, lastblock, inlen);

    /* XOR message to state */
    storm_pad(lastblock, out, inlen);
    for (i = 0; i < n; ++i)
    {
        S[i] ^= LOAD(lastblock + i * BYTES(STORM_W));
    }

#if defined(STORM_DEBUG)
    printf("DECRYPTING LASTBLOCK\n");
    printf("IN:\n");
    print_bytes(in, inlen);
    printf("OUT:\n");
    print_bytes(out, inlen);
    printf("STATE:\n");
    print_state(state);
    printf("MASK:\n");
    print_state(mask);
#endif
    burn(lastblock, 0, BYTES(STORM_B));
}

/* low-level interface functions */
void storm_absorb_data(storm_state_t state, storm_state_t mask, const unsigned char * in, size_t inlen)
{
    while (inlen >= BYTES(STORM_B))
    {
        storm_absorb_block(state, mask, in);
        inlen -= BYTES(STORM_B);
        in    += BYTES(STORM_B);
        storm_phi(mask);
    }
    if (inlen > 0)
    {
        storm_sigma(mask);
        storm_absorb_lastblock(state, mask, in, inlen);
    }
}

void storm_encrypt_data(storm_state_t state, storm_state_t mask, unsigned char * out, const unsigned char * in, size_t inlen)
{
    storm_lambda(mask);
    while (inlen >= BYTES(STORM_B))
    {
        storm_encrypt_block(state, mask, out, in);
        inlen -= BYTES(STORM_B);
        in    += BYTES(STORM_B);
        out   += BYTES(STORM_B);
        storm_phi(mask);
    }
    if (inlen > 0)
    {
        storm_sigma(mask);
        storm_encrypt_lastblock(state, mask, out, in, inlen);
    }
}

void storm_decrypt_data(storm_state_t state, storm_state_t mask, unsigned char * out, const unsigned char * in, size_t inlen)
{
    storm_lambda(mask);
    while (inlen >= BYTES(STORM_B))
    {
        storm_decrypt_block(state, mask, out, in);
        inlen -= BYTES(STORM_B);
        in    += BYTES(STORM_B);
        out   += BYTES(STORM_B);
        storm_phi(mask);
    }
    if (inlen > 0)
    {
        storm_sigma(mask);
        storm_decrypt_lastblock(state, mask, out, in, inlen);
    }
}

void storm_finalise(storm_state_t sa, storm_state_t se, storm_state_t mask, unsigned char *tag, size_t hlen, size_t mlen)
{
    size_t i, j;
    const size_t n = WORDS(STORM_B);
    storm_word_t * SA = sa->S;
    storm_word_t * SE = se->S;
    storm_word_t * L = mask->S;
    uint8_t block[BYTES(STORM_B)];

    /* determine how often to update mask depending on hlen and mlen */
    i = BYTES(STORM_B);
    j = 2 + ((mlen % i) + i - 1) / i;
    j = j - ((hlen % i) + i - 1) / i;

    for (i = 0; i < j; ++i)
    {
        storm_sigma(mask);
    }

    for (i = 0; i < n; ++i)
    {
        SE[i] ^= L[i];
    }

    storm_permute(se, STORM_L);

    for (i = 0; i < n; ++i)
    {
        SA[i] ^= SE[i] ^ L[i];
        STORE(block + i * BYTES(STORM_W), SA[i]);
    }
    memcpy(tag, block, BYTES(STORM_T));

#if defined(STORM_DEBUG)
    printf("EXTRACTING TAG:\n");
    print_bytes(tag, BYTES(STORM_T));
    printf("STATE:\n");
    print_state(sa);
    printf("MASK:\n");
    print_state(mask);
#endif
    burn(block, 0, BYTES(STORM_B));
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
    storm_state_t sa, se, la, le;

    /* init checksums and masks */
    memset(sa, 0, sizeof(storm_state_t));
    memset(se, 0, sizeof(storm_state_t));
    storm_init_mask(la, key, nonce);
    memcpy(le, la, sizeof(storm_state_t));

    /* absorb header */
    storm_absorb_data(sa, la, h, hlen);

    /* encrypt message */
    storm_encrypt_data(se, le, c, m, mlen);
    *clen = mlen + BYTES(STORM_T);

    /* finalise and extract tag */
    storm_finalise(sa, se, la, c + mlen, hlen, mlen);

    /* empty buffers */
    burn(sa, 0, sizeof(storm_state_t));
    burn(se, 0, sizeof(storm_state_t));
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
    storm_state_t sa, se, la, le;

    if (clen < BYTES(STORM_T)) { return result; }

    /* init checksums and masks */
    memset(sa, 0, sizeof(storm_state_t));
    memset(se, 0, sizeof(storm_state_t));
    storm_init_mask(la, key, nonce);
    memcpy(le, la, sizeof(storm_state_t));

    /* absorb header */
    storm_absorb_data(sa, la, h, hlen);

    /* decrypt message */
    storm_decrypt_data(se, le, m, c, clen - BYTES(STORM_T));
    *mlen = clen - BYTES(STORM_T);

    /* finalise and extract tag */
    storm_finalise(sa, se, la, tag, hlen, *mlen);

    /* verify tag */
    result = storm_verify_tag(c + clen - BYTES(STORM_T), tag);

    /* burn decrypted plaintext on authentication failure */
    if(result != 0) { burn(m, 0, *mlen); }

    /* empty buffers */
    burn(sa, 0, sizeof(storm_state_t));
    burn(se, 0, sizeof(storm_state_t));
    burn(la, 0, sizeof(storm_state_t));
    burn(le, 0, sizeof(storm_state_t));

    return result;
}
