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
} while(0)

static void storm_init_abs(__m256i B[4], const uint8_t * k, const uint8_t * n) {
  B[0] = _mm256_castsi128_si256(LOADU128(n));
  B[1] = _mm256_setzero_si256();
  B[2] = _mm256_set_epi64x(0, STORM_T, STORM_R, STORM_W);
  B[3] = LOADU256(k);
}

static void storm_init_enc(__m256i B[4], const uint8_t * k, const uint8_t * n) {
  B[0] = LOADU256(n);
  B[1] = _mm256_setzero_si256();
  B[2] = _mm256_set_epi64x(1, STORM_T, STORM_R, STORM_W);
  B[3] = LOADU256(k);
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

  storm_init_abs(B, k, n);

  while(hlen >= BYTES(STORM_B)) {
    __m256i H[4];
    V1_LOAD_BLOCK(H, h);
    V1_PERMUTE_F(B);
    V1_ACCUMULATE(B, H);
    h    += BYTES(STORM_B);
    hlen -= BYTES(STORM_B);
  }

  if(hlen > 0) {
    __m256i H[4];
    STORM_PAD(H, h, hlen);
    V1_PERMUTE_F(B);
    V1_ACCUMULATE(B, H);
  }

  while(mlen >= BYTES(STORM_B)) {
    __m256i M[4];
    V1_LOAD_BLOCK(M, m);
    V1_PERMUTE_F(B);
    V1_ACCUMULATE(B, M);
    m    += BYTES(STORM_B);
    mlen -= BYTES(STORM_B);
  }

  if(mlen > 0) {
    __m256i M[4];
    STORM_PAD(M, m, mlen);
    V1_PERMUTE_F(B);
    V1_ACCUMULATE(B, M);
  }

  V1_PERMUTE_F(B);
  B[0] = XOR256(B[0], _mm256_set_epi64x(0, 0, mlen_, hlen_));
  V1_PERMUTE_F(B);
  memcpy(tag, B, BYTES(STORM_T));
}


static void storm_encrypt_data(
  uint8_t * c,
  const unsigned char * m,
  size_t mlen,
  const uint8_t * k,
  const uint8_t * n)
{
  __m256i B[4];
  storm_init_enc(B, k, n);

  while(mlen >= BYTES(STORM_B - STORM_C)) {
    __m256i M[3];
    V1_LOAD_BLOCK_N(M, m, 3);
    V1_PERMUTE_F(B);
    V1_ACCUMULATE_N(B, M, 3);
    V1_STORE_BLOCK_N(c, B, 3);
    c    += BYTES(STORM_B - STORM_C);
    m    += BYTES(STORM_B - STORM_C);
    mlen -= BYTES(STORM_B - STORM_C);
  }

  if(mlen > 0) {
    __m256i M[4];
    STORM_PAD(M, m, mlen);
    V1_PERMUTE_F(B);
    V1_ACCUMULATE_N(B, M, 3);
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
  __m256i B[4];
  storm_init_enc(B, k, n);

  while(clen >= BYTES(STORM_B - STORM_C)) {
    __m256i C[3];
    V1_LOAD_BLOCK_N(C, c, 3);
    V1_PERMUTE_F(B);
    V1_ACCUMULATE_N(B, C, 3);
    V1_STORE_BLOCK_N(m, B, 3);
    V1_COPY_BLOCK_N(B, C, 3);
    m    += BYTES(STORM_B - STORM_C);
    c    += BYTES(STORM_B - STORM_C);
    clen -= BYTES(STORM_B - STORM_C);
  }

  if(clen > 0) {
    __m256i C[4];
    STORM_PAD(C, c, clen);
    V1_PERMUTE_F(B);
    V1_ACCUMULATE_N(B, C, 3);
    memcpy(m, B, clen);
  }
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

