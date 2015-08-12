#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <immintrin.h>

#define STORM_W 64
#define STORM_R 4
#define STORM_T (STORM_W *  4)
#define STORM_B (STORM_W * 16)
#define STORM_C (STORM_W *  4)

#include "v0.h"
#include "v1.h"
#include "v2.h"
#include "v4.h"

#define STORM_PAD(out, in, inlen) do { \
  memset(out, 0, BYTES(STORM_B));      \
  memcpy(out, in, inlen);              \
  out[inlen] = 0x01;                   \
  out[BYTES(STORM_B) - 1] |= 0x80;     \
} while(0)

static void storm_init_abs(__m256i B[4], const uint8_t * k, const uint8_t * n) {
  B[0] = _mm256_castsi128_si256(LOADU128(n));
  B[1] = LOADU256(k);
  B[2] = _mm256_setzero_si256();
  B[3] = _mm256_set_epi64x(0, STORM_T, STORM_R, STORM_W);
  V1_PERMUTE_F(B);
}

static void storm_init_enc(__m256i B[4], const uint8_t * k, const uint8_t * n) {
  B[0] = LOADU256(n);
  B[1] = LOADU256(k);
  B[2] = _mm256_setzero_si256();
  B[3] = _mm256_set_epi64x(1, STORM_T, STORM_R, STORM_W);
  V1_PERMUTE_F(B);
}

static void storm_compute_tag(uint8_t * tag,
  const uint8_t * h, size_t hlen,
  const uint8_t * m, size_t mlen,
  const unsigned char * k,
  const unsigned char * n)
{
  const size_t hlen_ = hlen;
  const size_t mlen_ = mlen;
  __m256i B[4];
  uint64_t L[16+4] = {0};

  storm_init_abs(B, k, n);
  memcpy(L, B, BYTES(STORM_B));
  V1_ZERO_BLOCK(B);

  while(hlen >= 4 * BYTES(STORM_B)) {
    __m256i H[16];
    V4_MASK_UPDATE_1(L);
    V4_LOAD_BLOCK(H, h);
    V4_BLOCKCIPHER_F(H, L);
    V4_ACCUMULATE(B, H);
    V4_MASK_UPDATE_2(L);
    h    += 4 * BYTES(STORM_B);
    hlen -= 4 * BYTES(STORM_B);
  }

  /* TODO: V2 */

  while(hlen >= BYTES(STORM_B)) {
    __m256i H[4];
    V1_MASK_UPDATE_1(L);

    V1_LOAD_BLOCK(H, h);
    V1_BLOCKCIPHER_F(H, L);
    V1_ACCUMULATE(B, H);

    V1_MASK_UPDATE_2(L);
    h    += BYTES(STORM_B);
    hlen -= BYTES(STORM_B);
  }

  if(hlen > 0) {
    uint8_t lastblock[BYTES(STORM_B)];
    __m256i H[4];
    STORM_PAD(lastblock, h, hlen);
    V1_MASK_UPDATE_1(L);
    V1_LOAD_BLOCK(H, lastblock);
    V1_BLOCKCIPHER_F(H, L);
    V1_ACCUMULATE(B, H);
    V1_MASK_UPDATE_2(L);
  }

  while(mlen >= 4 * BYTES(STORM_B)) {
    __m256i M[16];
    V4_MASK_UPDATE_1(L);
    V4_LOAD_BLOCK(M, m);
    V4_BLOCKCIPHER_F(M, L);
    V4_ACCUMULATE(B, M);
    V4_MASK_UPDATE_2(L);
    m    += 4 * BYTES(STORM_B);
    mlen -= 4 * BYTES(STORM_B);
  }

  /* TODO: V2 */

  while(mlen >= BYTES(STORM_B)) {
    __m256i M[4];
    V1_MASK_UPDATE_1(L);

    V1_LOAD_BLOCK(M, m);
    V1_BLOCKCIPHER_F(M, L);
    V1_ACCUMULATE(B, M);

    V1_MASK_UPDATE_2(L);
    m    += BYTES(STORM_B);
    mlen -= BYTES(STORM_B);
  }

  if(mlen > 0) {
    uint8_t lastblock[BYTES(STORM_B)];
    __m256i M[4];
    STORM_PAD(lastblock, m, mlen);
    V1_MASK_UPDATE_1(L);
    V1_LOAD_BLOCK(M, lastblock);
    V1_BLOCKCIPHER_F(M, L);
    V1_ACCUMULATE(B, M);
    V1_MASK_UPDATE_2(L);
  }

  {
    __m256i M[4];
    M[0] = _mm256_setzero_si256();
    M[1] = _mm256_setzero_si256();
    M[2] = _mm256_setzero_si256();
    M[3] = _mm256_set_epi64x(mlen_, hlen_, 0, 0);
    V1_BLOCKCIPHER_F(M, L);
    V1_ACCUMULATE(B, M);
  }
  memcpy(tag, B, BYTES(STORM_T));
}


static void storm_encrypt_data(
  uint8_t * c,
  const unsigned char * m,
  size_t mlen,
  const uint8_t * k,
  const uint8_t * n)
{
  __m256i counter = _mm256_set_epi64x(3, 2, 1, 0);
  __m256i K[4];
  __m256i K_[16];
  int i;

  storm_init_enc(K, k, n);

  if(mlen >= 4 * BYTES(STORM_B)) {
    K_[ 0] = _mm256_set1_epi64x(_mm256_extract_epi64(K[0], 0));
    K_[ 1] = _mm256_set1_epi64x(_mm256_extract_epi64(K[0], 1));
    K_[ 2] = _mm256_set1_epi64x(_mm256_extract_epi64(K[0], 2));
    K_[ 3] = _mm256_set1_epi64x(_mm256_extract_epi64(K[0], 3));

    K_[ 4] = _mm256_set1_epi64x(_mm256_extract_epi64(K[1], 0));
    K_[ 5] = _mm256_set1_epi64x(_mm256_extract_epi64(K[1], 1));
    K_[ 6] = _mm256_set1_epi64x(_mm256_extract_epi64(K[1], 2));
    K_[ 7] = _mm256_set1_epi64x(_mm256_extract_epi64(K[1], 3));

    K_[ 8] = _mm256_set1_epi64x(_mm256_extract_epi64(K[2], 0));
    K_[ 9] = _mm256_set1_epi64x(_mm256_extract_epi64(K[2], 1));
    K_[10] = _mm256_set1_epi64x(_mm256_extract_epi64(K[2], 2));
    K_[11] = _mm256_set1_epi64x(_mm256_extract_epi64(K[2], 3));

    K_[12] = _mm256_set1_epi64x(_mm256_extract_epi64(K[3], 0));
    K_[13] = _mm256_set1_epi64x(_mm256_extract_epi64(K[3], 1));
    K_[14] = _mm256_set1_epi64x(_mm256_extract_epi64(K[3], 2));
    K_[15] = _mm256_set1_epi64x(_mm256_extract_epi64(K[3], 3));
  }

  while(mlen >= 4 * BYTES(STORM_B)) {
    __m256i B[16], C[16];
    int i;
    memcpy(C, K_, 4 * BYTES(STORM_B));
    C[15] = XOR256(C[15], counter);

    V4_PERMUTE_F(C);

    for(i = 0; i < 16; ++i)
      C[i] = XOR256(C[i], K_[i]);
    V4_TRANSPOSE_B(B, C);
    for(i = 0; i < 16; ++i)
      STOREU256(&c[i * 32], XOR256(B[i], LOADU256(&m[i * 32])));

    counter = ADD256(counter, _mm256_set_epi64x(4, 4, 4, 4));
    c    += 4 * BYTES(STORM_B);
    m    += 4 * BYTES(STORM_B);
    mlen -= 4 * BYTES(STORM_B);
  }

  counter = _mm256_set_epi64x(_mm256_extract_epi64(counter, 0), 0, 0, 0);
  while(mlen >= BYTES(STORM_B)) {
    __m256i B[4];

    V1_COPY_BLOCK(B, K);
    B[3] = XOR256(B[3], counter);
    V1_PERMUTE_F(B);
    for(i = 0; i < 4; ++i) B[i] = XOR256(B[i], K[i]);

    for(i = 0; i < 4; ++i)
      STOREU256(&c[32 * i], XOR256(B[i], LOADU256(&m[32 * i])));

    counter = ADD256(counter, _mm256_set_epi64x(1, 0, 0, 0));
    c    += BYTES(STORM_B);
    m    += BYTES(STORM_B);
    mlen -= BYTES(STORM_B);
  }

  if(mlen > 0) {
    uint8_t lastblock[BYTES(STORM_B)];
    __m256i B[4];
    memcpy(lastblock, m, mlen);
    for(i = 0; i < 4; ++i) B[i] = K[i];
    B[3] = XOR256(B[3], counter);
    V1_PERMUTE_F(B);
    for(i = 0; i < 4; ++i) B[i] = XOR256(B[i], K[i]);
    for(i = 0; i < 4; ++i)
      B[i] = XOR256(B[i], LOADU256(&lastblock[32 * i]));
    memcpy(c, B, mlen);
  }
}

static void storm_decrypt_data(
  uint8_t * m,
  const unsigned char * c,
  size_t clen,
  const uint8_t * k,
  const uint8_t * n)
{
  storm_encrypt_data(m, c, clen, k, n);
}

/* high level interface functions */
void storm_aead_encrypt(
    unsigned char *c, size_t *clen,
    const unsigned char *h, size_t hlen,
    const unsigned char *m, size_t mlen,
    const unsigned char *n,
    const unsigned char *k
    )
{
  uint8_t t[BYTES(STORM_T)];
  storm_compute_tag(t, h, hlen, m, mlen, k, n);
  storm_encrypt_data(c, m, mlen, k, t);
  memcpy(c + mlen, t, BYTES(STORM_T));
  *clen = mlen + BYTES(STORM_T);
}

int storm_aead_decrypt(
    unsigned char *m, size_t *mlen,
    const unsigned char *h, size_t hlen,
    const unsigned char *c, size_t clen,
    const unsigned char *n,
    const unsigned char *k)
{
  uint8_t t[BYTES(STORM_T)];
  __m256i v;

  if (clen < BYTES(STORM_T))
    return -1;

  *mlen = clen - BYTES(STORM_T);
  storm_decrypt_data(m, c, *mlen, k, c + *mlen);
  storm_compute_tag(t, h, hlen, m, *mlen, k, n);

  v = _mm256_cmpeq_epi8(LOADU256(t), LOADU256(c + clen - BYTES(STORM_T)));
  return (( (_mm256_movemask_epi8(v) & 0xFFFFFFFFULL) + 1) >> 32) - 1;
}

