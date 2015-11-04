/*
    MRS - MEM AEAD source code package

    :copyright: (c) 2015 by Philipp Jovanovic and Samuel Neves
    :license: Creative Commons CC0 1.0
*/
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <immintrin.h>

#define MRS_W 64
#define MRS_L 4
#define MRS_T (MRS_W *  4)
#define MRS_B (MRS_W * 16)
#define MRS_C (MRS_W *  4)

#include "v0.h"
#include "v1.h"

#define MRS_PAD(out, in, inlen) do { \
  memset(out, 0, BYTES(MRS_B));      \
  memcpy(out, in, inlen);            \
} while(0)

static void mrs_init_abs(__m256i B[4], const uint8_t * k, const uint8_t * n) {
  B[0] = _mm256_castsi128_si256(LOADU128(n));
  B[1] = _mm256_setzero_si256();
  B[2] = _mm256_set_epi64x(0, MRS_T, MRS_L, 0);
  B[3] = LOADU256(k);
}

static void mrs_init_enc(__m256i B[4], const uint8_t * k, const uint8_t * n) {
  B[0] = LOADU256(n);
  B[1] = _mm256_setzero_si256();
  B[2] = _mm256_set_epi64x(1, MRS_T, MRS_L, 0);
  B[3] = LOADU256(k);
}

static void mrs_compute_tag(uint8_t * tag,
  const uint8_t * h, size_t hlen,
  const uint8_t * m, size_t mlen,
  const unsigned char * k,
  const unsigned char * n)
{
  const size_t hlen_ = hlen;
  const size_t mlen_ = mlen;
  __m256i B[4];

  mrs_init_abs(B, k, n);

  while(hlen >= BYTES(MRS_B)) {
    __m256i H[4];
    V1_LOAD_BLOCK(H, h);
    V1_PERMUTE_F(B);
    V1_ACCUMULATE(B, H);
    h    += BYTES(MRS_B);
    hlen -= BYTES(MRS_B);
  }

  if(hlen > 0) {
    __m256i H[4];
    MRS_PAD(H, h, hlen);
    V1_PERMUTE_F(B);
    V1_ACCUMULATE(B, H);
  }

  while(mlen >= BYTES(MRS_B)) {
    __m256i M[4];
    V1_LOAD_BLOCK(M, m);
    V1_PERMUTE_F(B);
    V1_ACCUMULATE(B, M);
    m    += BYTES(MRS_B);
    mlen -= BYTES(MRS_B);
  }

  if(mlen > 0) {
    __m256i M[4];
    MRS_PAD(M, m, mlen);
    V1_PERMUTE_F(B);
    V1_ACCUMULATE(B, M);
  }

  V1_PERMUTE_F(B);
  B[0] = XOR256(B[0], _mm256_set_epi64x(0, 0, mlen_, hlen_));
  V1_PERMUTE_F(B);
  memcpy(tag, B, BYTES(MRS_T));
}


static void mrs_encrypt_data(
  uint8_t * c,
  const unsigned char * m,
  size_t mlen,
  const uint8_t * k,
  const uint8_t * n)
{
  __m256i B[4];
  mrs_init_enc(B, k, n);

  while(mlen >= BYTES(MRS_B - MRS_C)) {
    __m256i M[3];
    V1_LOAD_BLOCK_N(M, m, 3);
    V1_PERMUTE_F(B);
    V1_ACCUMULATE_N(B, M, 3);
    V1_STORE_BLOCK_N(c, B, 3);
    c    += BYTES(MRS_B - MRS_C);
    m    += BYTES(MRS_B - MRS_C);
    mlen -= BYTES(MRS_B - MRS_C);
  }

  if(mlen > 0) {
    __m256i M[4];
    MRS_PAD(M, m, mlen);
    V1_PERMUTE_F(B);
    V1_ACCUMULATE_N(B, M, 3);
    memcpy(c, B, mlen);
  }
}

static void mrs_decrypt_data(
  uint8_t * m,
  const unsigned char * c,
  size_t clen,
  const uint8_t * k,
  const uint8_t * n)
{
  __m256i B[4];
  mrs_init_enc(B, k, n);

  while(clen >= BYTES(MRS_B - MRS_C)) {
    __m256i C[3];
    V1_LOAD_BLOCK_N(C, c, 3);
    V1_PERMUTE_F(B);
    V1_ACCUMULATE_N(B, C, 3);
    V1_STORE_BLOCK_N(m, B, 3);
    V1_COPY_BLOCK_N(B, C, 3);
    m    += BYTES(MRS_B - MRS_C);
    c    += BYTES(MRS_B - MRS_C);
    clen -= BYTES(MRS_B - MRS_C);
  }

  if(clen > 0) {
    __m256i C[4];
    MRS_PAD(C, c, clen);
    V1_PERMUTE_F(B);
    V1_ACCUMULATE_N(B, C, 3);
    memcpy(m, B, clen);
  }
}

/* high level interface functions */
void crypto_aead_encrypt(
    unsigned char *c, size_t *clen,
    const unsigned char *h, size_t hlen,
    const unsigned char *m, size_t mlen,
    const unsigned char *n,
    const unsigned char *k
    )
{
  uint8_t t[BYTES(MRS_T)];
  mrs_compute_tag(t, h, hlen, m, mlen, k, n);
  mrs_encrypt_data(c, m, mlen, k, t);
  memcpy(c + mlen, t, BYTES(MRS_T));
  *clen = mlen + BYTES(MRS_T);
}

int crypto_aead_decrypt(
    unsigned char *m, size_t *mlen,
    const unsigned char *h, size_t hlen,
    const unsigned char *c, size_t clen,
    const unsigned char *n,
    const unsigned char *k)
{
  uint8_t t[BYTES(MRS_T)];
  __m256i v;

  if (clen < BYTES(MRS_T))
    return -1;

  *mlen = clen - BYTES(MRS_T);
  mrs_decrypt_data(m, c, *mlen, k, c + *mlen);
  mrs_compute_tag(t, h, hlen, m, *mlen, k, n);

  v = _mm256_cmpeq_epi8(LOADU256(t), LOADU256(c + clen - BYTES(MRS_T)));
  return (( (_mm256_movemask_epi8(v) & 0xFFFFFFFFULL) + 1) >> 32) - 1;
}

