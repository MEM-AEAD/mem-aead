#include "opp.h"

#if defined(OPP_DEBUG)
#include <stdio.h>
void print_state(uint64_t S[16]);
void print_bytes(const uint8_t * in, size_t inlen);
#endif

#define OPP_W 64           /* word size */
#define OPP_L 4            /* double round number */
#define OPP_T (OPP_W *  4) /* tag size */
#define OPP_N (OPP_W *  2) /* nonce size */
#define OPP_K (OPP_W *  4) /* key size */
#define OPP_B (OPP_W * 16) /* permutation width */

/* Workaround for C89 compilers */
#if !defined(__cplusplus) && (!defined(__STDC_VERSION__) || __STDC_VERSION__ < 199901L)
  #if   defined(_MSC_VER)
    #define OPP_INLINE __inline
  #elif defined(__GNUC__)
    #define OPP_INLINE __inline__
  #else
    #define OPP_INLINE
  #endif
#else
  #define OPP_INLINE inline
#endif

#define BITS(x) (sizeof(x) * CHAR_BIT)
#define BYTES(x) (((x) + 7) / 8)
#define WORDS(x) (((x) + (OPP_W-1)) / OPP_W)

static OPP_INLINE uint64_t load64(const void * in)
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

static OPP_INLINE void store64(void * out, const uint64_t v)
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
static OPP_INLINE void F(opp_word_t S[16])
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
static OPP_INLINE void FI(opp_word_t S[16])
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

static OPP_INLINE void opp_permute(opp_state_t state, size_t R)
{
    opp_word_t * S = state->S;
    size_t i;
    for(i = 0; i < R; ++i)
    {
        F(S);
    }
}

static OPP_INLINE void opp_permute_inverse(opp_state_t state, size_t R)
{
    opp_word_t * S = state->S;
    size_t i;
    for(i = 0; i < R; ++i)
    {
        FI(S);
    }
}

static OPP_INLINE void opp_pad(uint8_t * out, const uint8_t * in, const size_t inlen)
{
    memset(out, 0, BYTES(OPP_B));
    memcpy(out, in, inlen);
    out[inlen] = 0x01;
}

static OPP_INLINE void opp_init_mask(opp_state_t mask, const unsigned char * k, const unsigned char * n)
{
    opp_word_t * L = mask->S;

    L[ 0] = LOAD(n + 0 * BYTES(OPP_W));
    L[ 1] = LOAD(n + 1 * BYTES(OPP_W));
    L[ 2] = 0;
    L[ 3] = 0;

    L[ 4] = 0;
    L[ 5] = 0;
    L[ 6] = 0;
    L[ 7] = 0;

    L[ 8] = 0;
    L[ 9] = 0;
    L[10] = OPP_L;
    L[11] = OPP_T;

    L[12] = LOAD(k + 0 * BYTES(OPP_W));
    L[13] = LOAD(k + 1 * BYTES(OPP_W));
    L[14] = LOAD(k + 2 * BYTES(OPP_W));
    L[15] = LOAD(k + 3 * BYTES(OPP_W));

    /* apply permutation */
    opp_permute(mask, OPP_L);

#if defined(OPP_DEBUG)
    printf("SETUP MASK:\n");
    print_state(mask->S);
#endif

}

/* alpha(x) = phi(x) */
static OPP_INLINE void opp_alpha(opp_state_t mask)
{
    size_t i;
    opp_word_t * L = mask->S;
    opp_word_t t = ROTL(L[0], 53) ^ (L[5] << 13);
    for (i = 0; i < WORDS(OPP_B) - 1; ++i)
    {
        L[i] = L[i+1];
    }
    L[15] = t;
}

/* beta(x) = phi(x) ^ x */
static OPP_INLINE void opp_beta(opp_state_t mask)
{
    size_t i;
    opp_word_t * L = mask->S;
    opp_word_t t = ROTL(L[0], 53) ^ (L[5] << 13);
    for (i = 0; i < WORDS(OPP_B) - 1; ++i)
    {
        L[i] ^= L[i+1];
    }
    L[15] ^= t;
}

/* gamma(x) = phi^2(x) ^ phi(x) ^ x */
static OPP_INLINE void opp_gamma(opp_state_t mask)
{
    size_t i;
    opp_word_t * L = mask->S;
    opp_word_t t0 = ROTL(L[0], 53) ^ (L[5] << 13);
    opp_word_t t1 = ROTL(L[1], 53) ^ (L[6] << 13);
    for (i = 0; i < WORDS(OPP_B) - 2; ++i)
    {
        L[i] ^= L[i+1] ^ L[i+2];
    }
    L[14] ^= (L[15] ^ t0);
    L[15] ^= (t0 ^ t1);
}

static OPP_INLINE void opp_absorb_block(opp_state_t state, opp_state_t mask, const uint8_t * in)
{
    size_t i;
    const size_t n = WORDS(OPP_B);
    opp_state_t block;
    opp_word_t * B = block->S;
    opp_word_t * S = state->S;
    opp_word_t * L = mask->S;

    /* load data and XOR mask */
    for (i = 0; i < n; ++i)
    {
        B[i] = LOAD(in + i * BYTES(OPP_W)) ^ L[i];
    }

    /* apply permutation */
    opp_permute(block, OPP_L);

    /* XOR mask and absorb into state */
    for (i = 0; i < n; ++i)
    {
        S[i] ^= B[i] ^ L[i];
    }

#if defined(OPP_DEBUG)
    printf("ABSORBING BLOCK\n");
    printf("IN:\n");
    print_bytes(in, BYTES(OPP_B));
    printf("\nSTATE:\n");
    print_state(state->S);
    printf("MASK:\n");
    print_state(mask->S);
#endif
}

static OPP_INLINE void opp_absorb_lastblock(opp_state_t state, opp_state_t mask, const uint8_t * in, size_t inlen)
{
    size_t i;
    const size_t n = WORDS(OPP_B);
    opp_state_t block;
    opp_word_t * B = block->S;
    opp_word_t * S = state->S;
    opp_word_t * L = mask->S;
    uint8_t lastblock[BYTES(OPP_B)];

    opp_pad(lastblock, in, inlen);

    /* load data and XOR mask */
    for (i = 0; i < n; ++i)
    {
        B[i] = LOAD(lastblock + i * BYTES(OPP_W)) ^ L[i];
    }

    /* apply permutation */
    opp_permute(block, OPP_L);

    /* XOR mask and absorb into state */
    for (i = 0; i < n; ++i)
    {
        S[i] ^= B[i] ^ L[i];
    }

#if defined(OPP_DEBUG)
    printf("ABSORBING LASTBLOCK\n");
    printf("IN:\n");
    print_bytes(in, inlen);
    printf("\nSTATE:\n");
    print_state(state->S);
    printf("MASK:\n");
    print_state(mask->S);
#endif
    burn(lastblock, 0, BYTES(OPP_B));
}

static OPP_INLINE void opp_encrypt_block(opp_state_t state, opp_state_t mask, uint8_t * out, const uint8_t * in)
{
    size_t i;
    const size_t n = WORDS(OPP_B);
    opp_state_t block;
    opp_word_t * B = block->S;
    opp_word_t * S = state->S;
    opp_word_t * L = mask->S;

    /* load message and XOR mask */
    for (i = 0; i < n; ++i)
    {
        B[i] = LOAD(in + i * BYTES(OPP_W)) ^ L[i];
    }

    /* apply permutation */
    opp_permute(block, OPP_L);

    /* XOR mask to block and store as ciphertext, XOR message to state */
    for (i = 0; i < n; ++i)
    {
        STORE(out + i * BYTES(OPP_W), B[i] ^ L[i]);
        S[i] ^= LOAD(in + i * BYTES(OPP_W));
    }

#if defined(OPP_DEBUG)
    printf("ENCRYPTING BLOCK\n");
    printf("IN:\n");
    print_bytes(in, BYTES(OPP_B));
    printf("OUT:\n");
    print_bytes(out, BYTES(OPP_B));
    printf("STATE:\n");
    print_state(state->S);
    printf("MASK:\n");
    print_state(mask->S);
#endif
}

static OPP_INLINE void opp_encrypt_lastblock(opp_state_t state, opp_state_t mask, uint8_t * out, const uint8_t * in, size_t inlen)
{
    size_t i;
    const size_t n = WORDS(OPP_B);
    opp_state_t block;
    opp_word_t * B = block->S;
    opp_word_t * S = state->S;
    opp_word_t * L = mask->S;
    uint8_t lastblock[BYTES(OPP_B)];

    /* load block with mask */
    for (i = 0; i < n; ++i)
    {
        B[i] = L[i];
    }

    /* apply permutation */
    opp_permute(block, OPP_L);

    /* XOR padded message to state, XOR padded message and mask to block, and extract ciphertext */
    opp_pad(lastblock, in, inlen);
    for (i = 0; i < WORDS(OPP_B); ++i)
    {
        S[i] ^= LOAD(lastblock + i * BYTES(OPP_W));
        STORE(lastblock + i * BYTES(OPP_W), B[i] ^ L[i] ^ LOAD(lastblock + i * BYTES(OPP_W)));
    }
    memcpy(out, lastblock, inlen);

#if defined(OPP_DEBUG)
    printf("ENCRYPTING LASTBLOCK\n");
    printf("IN:\n");
    print_bytes(in, inlen);
    printf("OUT:\n");
    print_bytes(out, inlen);
    printf("STATE:\n");
    print_state(state->S);
    printf("MASK:\n");
    print_state(mask->S);
#endif
    burn(lastblock, 0, BYTES(OPP_B));
}

static OPP_INLINE void opp_decrypt_block(opp_state_t state, opp_state_t mask, uint8_t * out, const uint8_t * in)
{
    size_t i;
    const size_t n = WORDS(OPP_B);
    opp_state_t block;
    opp_word_t * B = block->S;
    opp_word_t * S = state->S;
    opp_word_t * L = mask->S;

    /* load ciphertext and XOR mask */
    for (i = 0; i < n; ++i)
    {
        B[i] = LOAD(in + i * BYTES(OPP_W)) ^ L[i];
    }

    /* apply inverse permutation */
    opp_permute_inverse(block, OPP_L);

    /* XOR ciphertext to state, XOR mask to block, and extract message */
    for (i = 0; i < n; ++i)
    {
        STORE(out + i * BYTES(OPP_W), B[i] ^ L[i]);
        S[i] ^= LOAD(out + i * BYTES(OPP_W));
    }

#if defined(OPP_DEBUG)
    printf("DECRYPTING BLOCK\n");
    printf("IN:\n");
    print_bytes(in, BYTES(OPP_B));
    printf("OUT:\n");
    print_bytes(out, BYTES(OPP_B));
    printf("STATE:\n");
    print_state(state->S);
    printf("MASK:\n");
    print_state(mask->S);
#endif
}

static OPP_INLINE void opp_decrypt_lastblock(opp_state_t state, opp_state_t mask, uint8_t * out, const uint8_t * in, size_t inlen)
{
    size_t i;
    const size_t n = WORDS(OPP_B);
    opp_state_t block;
    opp_word_t * B = block->S;
    opp_word_t * S = state->S;
    opp_word_t * L = mask->S;
    uint8_t lastblock[BYTES(OPP_B)];

    /* load block with key */
    for (i = 0; i < n; ++i)
    {
        B[i] = L[i];
    }

    /* apply permutation */
    opp_permute(block, OPP_L);

    /* XOR padded ciphertext and key to block, store message */
    opp_pad(lastblock, in, inlen);
    for (i = 0; i < n; ++i)
    {
        STORE(lastblock + i * BYTES(OPP_W), B[i] ^ L[i] ^ LOAD(lastblock + i * BYTES(OPP_W)));
    }
    memcpy(out, lastblock, inlen);

    /* XOR message to state */
    opp_pad(lastblock, out, inlen);
    for (i = 0; i < n; ++i)
    {
        S[i] ^= LOAD(lastblock + i * BYTES(OPP_W));
    }

#if defined(OPP_DEBUG)
    printf("DECRYPTING LASTBLOCK\n");
    printf("IN:\n");
    print_bytes(in, inlen);
    printf("OUT:\n");
    print_bytes(out, inlen);
    printf("STATE:\n");
    print_state(state->S);
    printf("MASK:\n");
    print_state(mask->S);
#endif
    burn(lastblock, 0, BYTES(OPP_B));
}

/* low-level interface functions */
void opp_absorb_data(opp_state_t state, opp_state_t mask, const unsigned char * in, size_t inlen)
{
    while (inlen >= BYTES(OPP_B))
    {
        opp_absorb_block(state, mask, in);
        inlen -= BYTES(OPP_B);
        in    += BYTES(OPP_B);
        opp_alpha(mask);
    }
    if (inlen > 0)
    {
        opp_beta(mask);
        opp_absorb_lastblock(state, mask, in, inlen);
    }
}

void opp_encrypt_data(opp_state_t state, opp_state_t mask, unsigned char * out, const unsigned char * in, size_t inlen)
{
    opp_gamma(mask);
    while (inlen >= BYTES(OPP_B))
    {
        opp_encrypt_block(state, mask, out, in);
        inlen -= BYTES(OPP_B);
        in    += BYTES(OPP_B);
        out   += BYTES(OPP_B);
        opp_alpha(mask);
    }
    if (inlen > 0)
    {
        opp_beta(mask);
        opp_encrypt_lastblock(state, mask, out, in, inlen);
    }
}

void opp_decrypt_data(opp_state_t state, opp_state_t mask, unsigned char * out, const unsigned char * in, size_t inlen)
{
    opp_gamma(mask);
    while (inlen >= BYTES(OPP_B))
    {
        opp_decrypt_block(state, mask, out, in);
        inlen -= BYTES(OPP_B);
        in    += BYTES(OPP_B);
        out   += BYTES(OPP_B);
        opp_alpha(mask);
    }
    if (inlen > 0)
    {
        opp_beta(mask);
        opp_decrypt_lastblock(state, mask, out, in, inlen);
    }
}

void opp_finalise(opp_state_t sa, opp_state_t se, opp_state_t mask, unsigned char *tag, size_t hlen, size_t mlen)
{
    size_t i, j;
    const size_t n = WORDS(OPP_B);
    opp_word_t * SA = sa->S;
    opp_word_t * SE = se->S;
    opp_word_t * L = mask->S;
    uint8_t block[BYTES(OPP_B)];

    /* determine how often to update the mask depending on hlen and mlen */
    i = BYTES(OPP_B);
    j = 2 + ( ( mlen % i ) + i - 1 ) / i - ( ( hlen % i ) + i - 1 ) / i;

    for (i = 0; i < j; ++i)
    {
        opp_beta(mask);
    }

    for (i = 0; i < n; ++i)
    {
        SE[i] ^= L[i];
    }

    opp_permute(se, OPP_L);

    for (i = 0; i < n; ++i)
    {
        SA[i] ^= SE[i] ^ L[i];
        STORE(block + i * BYTES(OPP_W), SA[i]);
    }
    memcpy(tag, block, BYTES(OPP_T));

#if defined(OPP_DEBUG)
    printf("EXTRACTING TAG:\n");
    print_bytes(tag, BYTES(OPP_T));
    printf("STATE:\n");
    print_state(sa->S);
    printf("MASK:\n");
    print_state(mask->S);
#endif
    burn(block, 0, BYTES(OPP_B));
}

int opp_verify_tag(const unsigned char * tag1, const unsigned char * tag2)
{
    unsigned acc = 0;
    size_t i;

    for(i = 0; i < BYTES(OPP_T); ++i)
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
    opp_state_t sa, se, la, le;

    /* init checksums and masks */
    memset(sa, 0, sizeof(opp_state_t));
    memset(se, 0, sizeof(opp_state_t));
    opp_init_mask(la, key, nonce);
    memcpy(le, la, sizeof(opp_state_t));

    /* absorb header */
    opp_absorb_data(sa, la, h, hlen);

    /* encrypt message */
    opp_encrypt_data(se, le, c, m, mlen);
    *clen = mlen + BYTES(OPP_T);

    /* finalise and extract tag */
    opp_finalise(sa, se, la, c + mlen, hlen, mlen);

    /* empty buffers */
    burn(sa, 0, sizeof(opp_state_t));
    burn(se, 0, sizeof(opp_state_t));
    burn(la, 0, sizeof(opp_state_t));
    burn(le, 0, sizeof(opp_state_t));
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
    unsigned char tag[BYTES(OPP_T)];
    opp_state_t sa, se, la, le;

    if (clen < BYTES(OPP_T)) { return result; }

    /* init checksums and masks */
    memset(sa, 0, sizeof(opp_state_t));
    memset(se, 0, sizeof(opp_state_t));
    opp_init_mask(la, key, nonce);
    memcpy(le, la, sizeof(opp_state_t));

    /* absorb header */
    opp_absorb_data(sa, la, h, hlen);

    /* decrypt message */
    opp_decrypt_data(se, le, m, c, clen - BYTES(OPP_T));
    *mlen = clen - BYTES(OPP_T);

    /* finalise and extract tag */
    opp_finalise(sa, se, la, tag, hlen, *mlen);

    /* verify tag */
    result = opp_verify_tag(c + clen - BYTES(OPP_T), tag);

    /* burn decrypted plaintext on authentication failure */
    if(result != 0) { burn(m, 0, *mlen); }

    /* empty buffers */
    burn(sa, 0, sizeof(opp_state_t));
    burn(se, 0, sizeof(opp_state_t));
    burn(la, 0, sizeof(opp_state_t));
    burn(le, 0, sizeof(opp_state_t));

    return result;
}
