#include "mro.h"

#if defined(MRO_DEBUG)
#include "../debug.h"
#endif

#define MRO_W 64           /* word size */
#define MRO_L 4            /* round number */
#define MRO_T (MRO_W *  4) /* tag size */
#define MRO_N (MRO_W *  2) /* nonce size */
#define MRO_K (MRO_W *  4) /* key size */
#define MRO_B (MRO_W * 16) /* permutation width */

/* workaround for C89 compilers */
#if !defined(__cplusplus) && (!defined(__STDC_VERSION__) || __STDC_VERSION__ < 199901L)
  #if   defined(_MSC_VER)
    #define MRO_INLINE __inline
  #elif defined(__GNUC__)
    #define MRO_INLINE __inline__
  #else
    #define MRO_INLINE
  #endif
#else
  #define MRO_INLINE inline
#endif

#define BITS(x) (sizeof(x) * CHAR_BIT)
#define BYTES(x) (((x) + 7) / 8)
#define WORDS(x) (((x) + (MRO_W-1)) / MRO_W)

static MRO_INLINE uint64_t load64(const void * in)
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

static MRO_INLINE void store64(void * out, const uint64_t v)
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

/* double round */
static MRO_INLINE void F(mro_word_t S[16])
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

static MRO_INLINE void mro_permute(mro_state_t state, size_t R)
{
    mro_word_t * S = state->S;
    size_t i;
    for (i = 0; i < R; ++i)
    {
        F(S);
    }
}

static MRO_INLINE void mro_pad(uint8_t * out, const uint8_t * in, const size_t inlen)
{
    memset(out, 0, BYTES(MRO_B));
    memcpy(out, in, inlen);
    out[inlen] = 0x01;
}

static MRO_INLINE void mro_init_mask(mro_state_t mask, const unsigned char * k, const unsigned char * n)
{
    mro_word_t * L = mask->S;

    L[ 0] = LOAD(n + 0 * BYTES(MRO_W));
    L[ 1] = LOAD(n + 1 * BYTES(MRO_W));
    L[ 2] = 0;
    L[ 3] = 0;

    L[ 4] = 0;
    L[ 5] = 0;
    L[ 6] = 0;
    L[ 7] = 0;

    L[ 8] = 0;
    L[ 9] = 0;
    L[10] = MRO_L;
    L[11] = MRO_T;

    L[12] = LOAD(k + 0 * BYTES(MRO_W));
    L[13] = LOAD(k + 1 * BYTES(MRO_W));
    L[14] = LOAD(k + 2 * BYTES(MRO_W));
    L[15] = LOAD(k + 3 * BYTES(MRO_W));

    /* apply permutation */
    mro_permute(mask, MRO_L);

#if defined(MRO_DEBUG)
    printf("SETUP MASK:\n");
    print_state(mask);
#endif
}

/* alpha(x) = phi(x) */
static MRO_INLINE void mro_alpha(mro_state_t mask)
{
    size_t i;
    mro_word_t * L = mask->S;
    mro_word_t t = ROTL(L[0], 53) ^ (L[5] << 13);
    for (i = 0; i < WORDS(MRO_B) - 1; ++i)
    {
        L[i] = L[i+1];
    }
    L[15] = t;
}

/* beta(x) = phi(x) ^ x */
static MRO_INLINE void mro_beta(mro_state_t mask)
{
    size_t i;
    mro_word_t * L = mask->S;
    mro_word_t t = ROTL(L[0], 53) ^ (L[5] << 13);
    for (i = 0; i < WORDS(MRO_B) - 1; ++i)
    {
        L[i] ^= L[i+1];
    }
    L[15] ^= t;
}

/* gamma(x) = phi^2(x) ^ phi(x) ^ x */
static MRO_INLINE void mro_gamma(mro_state_t mask)
{
    size_t i;
    mro_word_t * L = mask->S;
    mro_word_t t0 = ROTL(L[0], 53) ^ (L[5] << 13);
    mro_word_t t1 = ROTL(L[1], 53) ^ (L[6] << 13);
    for (i = 0; i < WORDS(MRO_B) - 2; ++i)
    {
        L[i] ^= L[i+1] ^ L[i+2];
    }
    L[14] ^= (L[15] ^ t0);
    L[15] ^= (t0 ^ t1);
}

static MRO_INLINE void mro_absorb_block(mro_state_t state, mro_state_t mask, const uint8_t * in)
{
    size_t i;
    mro_state_t block;
    mro_word_t * B = block->S;
    mro_word_t * S = state->S;
    mro_word_t * L = mask->S;

    /* load data and XOR mask */
    memset(block, 0, sizeof(mro_state_t));
    for (i = 0; i < WORDS(MRO_B); ++i)
    {
        B[i] = LOAD(in + i * BYTES(MRO_W)) ^ L[i];
    }

    /* apply permutation */
    mro_permute(block, MRO_L);

    /* XOR mask and absorb into S */
    for (i = 0; i < WORDS(MRO_B); ++i)
    {
        S[i] ^=  B[i] ^ L[i];
    }

#if defined(MRO_DEBUG)
    printf("ABSORBING BLOCK\n");
    printf("IN:\n");
    print_bytes(in, BYTES(MRO_B));
    printf("\nSTATE:\n");
    print_state(state);
    printf("MASK:\n");
    print_state(mask);
#endif
}

static MRO_INLINE void mro_absorb_lastblock(mro_state_t state, mro_state_t mask, const uint8_t * in, size_t inlen)
{
    uint8_t block[BYTES(MRO_B)];
    mro_pad(block, in, inlen);
    mro_absorb_block(state, mask, block);
    burn(block, 0, BYTES(MRO_B));
}

static MRO_INLINE void mro_encrypt_block(mro_state_t mask, mro_state_t tag, size_t block_nr, uint8_t * out, const uint8_t * in)
{
    size_t i;
    mro_state_t block;
    mro_word_t * B = block->S;
    mro_word_t * L = mask->S;
    mro_word_t * T = tag->S;

    /* load mask and XOR authentication tag and block counter */
    for (i = 0; i < WORDS(MRO_B); ++i)
    {
        B[i] = L[i];
    }
    B[0] ^= T[0];
    B[1] ^= T[1];
    B[2] ^= T[2];
    B[3] ^= T[3];
    B[15] ^= block_nr;

    /* apply permutation */
    mro_permute(block, MRO_L);

    /* encrypt block */
    for (i = 0; i < WORDS(MRO_B); ++i)
    {
        B[i] ^= LOAD(in + i * BYTES(MRO_W)) ^ L[i];
        STORE(out + i * BYTES(MRO_W), B[i]);
    }

#if defined(MRO_DEBUG)
    printf("ENCRYPTING BLOCK\n");
    printf("IN:\n");
    print_bytes(in, BYTES(MRO_B));
    printf("OUT:\n");
    print_bytes(out, BYTES(MRO_B));
    printf("MASK:\n");
    print_state(mask);
#endif
}

static MRO_INLINE void mro_encrypt_lastblock(mro_state_t mask, mro_state_t tag, size_t block_nr, uint8_t * out, const uint8_t * in, size_t inlen)
{
    uint8_t block[BYTES(MRO_B)];
    memset(block, 0, BYTES(MRO_B));
    memcpy(block, in, inlen);
    mro_encrypt_block(mask, tag, block_nr, block, block);
    memcpy(out, block, inlen);
    burn(block, 0, BYTES(MRO_B));
}

/* low-level interface functions */
void mro_absorb_data(mro_state_t state, mro_state_t mask, const unsigned char * in, size_t inlen, tag_t flag)
{
    if (flag)
    {
        mro_beta(mask);
    }
    while (inlen >= BYTES(MRO_B))
    {
        mro_absorb_block(state, mask, in);
        inlen -= BYTES(MRO_B);
        in    += BYTES(MRO_B);
        mro_alpha(mask);
    }
    if (inlen > 0)
    {
        mro_absorb_lastblock(state, mask, in, inlen);
    }
}

void mro_encrypt_data(mro_state_t mask, mro_state_t tag, unsigned char * out, const unsigned char * in, size_t inlen)
{
    size_t i = 0;
    mro_gamma(mask);
    while (inlen >= BYTES(MRO_B))
    {
        mro_encrypt_block(mask, tag, i, out, in);
        inlen -= BYTES(MRO_B);
        in    += BYTES(MRO_B);
        out   += BYTES(MRO_B);
        i += 1;
    }
    if (inlen > 0)
    {
        mro_encrypt_lastblock(mask, tag, i, out, in, inlen);
    }
}

void mro_decrypt_data(mro_state_t mask, mro_state_t tag, unsigned char * out, const unsigned char * in, size_t inlen)
{
    mro_encrypt_data(mask, tag, out, in, inlen);
}

static MRO_INLINE void mro_finalise(mro_state_t state, mro_state_t mask, size_t hlen, size_t mlen)
{
    size_t i;
    mro_word_t * S = state->S;
    mro_word_t * L = mask->S;

    mro_beta(mask);
    mro_beta(mask);

    S[14] ^= hlen;
    S[15] ^= mlen;

    for (i = 0; i < WORDS(MRO_B); ++i)
    {
        S[i] ^= L[i];
    }

    /* apply permutation */
    mro_permute(state, MRO_L);

    for (i = 0; i < WORDS(MRO_B); ++i)
    {
        S[i] ^= L[i];
    }
}

void mro_output_tag(mro_state_t state, unsigned char * tag)
{
    size_t i = 0;
    mro_word_t * S = state->S;
    uint8_t block[BYTES(MRO_T)];

    for (i = 0; i < WORDS(MRO_T); ++i)
    {
        STORE(block + i * BYTES(MRO_W), S[i]);
    }
    memcpy(tag, block, BYTES(MRO_T));
    burn(block, 0, BYTES(MRO_T));

#if defined(MRO_DEBUG)
    printf("EXTRACTING TAG:\n");
    print_state(state);
    printf("TAG:\n");
    print_bytes(tag, BYTES(MRO_T));
#endif
}

int mro_verify_tag(const unsigned char * tag1, const unsigned char * tag2)
{
    unsigned acc = 0;
    size_t i;

    for (i = 0; i < BYTES(MRO_T); ++i)
    {
        acc |= tag1[i] ^ tag2[i];
    }
    return (((acc - 1) >> 8) & 1) - 1;
}


/* high level interface functions */
void mro_aead_encrypt(
    unsigned char *c, size_t *clen,
    const unsigned char *h, size_t hlen,
    const unsigned char *m, size_t mlen,
    const unsigned char *nonce,
    const unsigned char *key
    )
{
    mro_state_t state, la, le;

    memset(state, 0, sizeof(mro_state_t));
    mro_init_mask(le, key, nonce);

    /* absorb header and message */
    memcpy(la, le, sizeof(mro_state_t));
    mro_absorb_data(state, la, h, hlen, ABS_AD);

    memcpy(la, le, sizeof(mro_state_t));
    mro_absorb_data(state, la, m, mlen, ABS_MSG);

    memcpy(la, le, sizeof(mro_state_t));
    mro_finalise(state, la, hlen, mlen);

    /* extract tag */
    mro_output_tag(state, c + mlen);
    *clen = mlen + BYTES(MRO_T);

    /* encrypt message */
    mro_encrypt_data(le, state, c, m, mlen);

    /* empty buffers */
    burn(state, 0, sizeof(mro_state_t));
    burn(la, 0, sizeof(mro_state_t));
    burn(le, 0, sizeof(mro_state_t));
}

int mro_aead_decrypt(
    unsigned char *m, size_t *mlen,
    const unsigned char *h, size_t hlen,
    const unsigned char *c, size_t clen,
    const unsigned char *nonce,
    const unsigned char *key
    )
{
    int result = -1;
    unsigned char tag[BYTES(MRO_T)];
    mro_state_t state, la, le;
    mro_word_t * T = state->S;

    if (clen < BYTES(MRO_T)) { return -1; }

    mro_init_mask(le, key, nonce);
    memcpy(la, le, sizeof(mro_state_t));

    *mlen = clen - BYTES(MRO_T);

    /* store received tag temporarily in the first 4 state words */
    T[ 0] = LOAD(c + *mlen + 0 * BYTES(MRO_W));
    T[ 1] = LOAD(c + *mlen + 1 * BYTES(MRO_W));
    T[ 2] = LOAD(c + *mlen + 2 * BYTES(MRO_W));
    T[ 3] = LOAD(c + *mlen + 3 * BYTES(MRO_W));

    /* decrypt message */
    mro_decrypt_data(le, state, m, c, clen - BYTES(MRO_T));

    /* reset state */
    memset(state, 0, sizeof(mro_state_t));

    /* absorb header and message */
    memcpy(le, la, sizeof(mro_state_t));
    mro_absorb_data(state, la, h, hlen, ABS_AD);

    memcpy(la, le, sizeof(mro_state_t));
    mro_absorb_data(state, la, m, *mlen, ABS_MSG);

    memcpy(la, le, sizeof(mro_state_t));
    mro_finalise(state, la, hlen, *mlen);

    /* extract tag */
    mro_output_tag(state, tag);

    /* verify tag */
    result = mro_verify_tag(c + clen - BYTES(MRO_T), tag);

    /* burn decrypted plaintext on authentication failure */
    if (result != 0) { burn(m, 0, *mlen); }

    /* empty buffers */
    burn(state, 0, sizeof(mro_state_t));
    burn(la, 0, sizeof(mro_state_t));
    burn(le, 0, sizeof(mro_state_t));

    return result;
}
