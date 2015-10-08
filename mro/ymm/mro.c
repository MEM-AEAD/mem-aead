#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <immintrin.h>

#define MRO_W 64
#define MRO_L 4
#define MRO_T (MRO_W *  4)
#define MRO_B (MRO_W * 16)
#define MRO_C (MRO_W *  4)

#include "v0.h"
#include "v1.h"
#include "v2.h"
#include "v4.h"

#define MRO_PAD(out, in, inlen) do { \
  memset(out, 0, BYTES(MRO_B));      \
  memcpy(out, in, inlen);            \
  out[inlen] = 0x01;                 \
} while(0)

static void mro_kdf(uint64_t Ka[16+4], uint64_t Ke[16+4], const uint8_t * k, const uint8_t * n) {
  __m256i B[4];
  B[0] = _mm256_castsi128_si256(LOADU128(n));
  B[1] = _mm256_setzero_si256();
  B[2] = _mm256_set_epi64x(MRO_T, MRO_L, 0, 0);
  B[3] = LOADU256(k);

  V1_PERMUTE_F(B);

  memcpy(Ka, B, sizeof(B));
  memcpy(Ke, B, sizeof(B));
  V1_GAMMA_UPDATE(Ke);
}

static void mro_compute_tag(void * tag,
  const uint8_t * h, size_t hlen,
  const uint8_t * m, size_t mlen,
  uint64_t L[16+4])
{
  const size_t hlen_ = hlen;
  const size_t mlen_ = mlen;
  __m256i B[4];
  uint64_t L_[16];

  memcpy(L_, L, BYTES(MRO_B)); /* we'll need it later */
  /* uint64_t L[16+4] = {0}; */

  /* mro_init_abs(B, k, n); */
  /* memcpy(L, B, BYTES(MRO_B)); */
  V1_ZERO_BLOCK(B);

  /* Absorb AD */
  while(hlen >= 4 * BYTES(MRO_B)) {
    __m256i H[16];
    V4_ALPHA_UPDATE_1(L);
    V4_LOAD_BLOCK(H, h);
    V4_BLOCKCIPHER_F(H, L);
    V4_ACCUMULATE(B, H);
    V4_ALPHA_UPDATE_2(L);
    h    += 4 * BYTES(MRO_B);
    hlen -= 4 * BYTES(MRO_B);
  }

  /* TODO: V2 */

  while(hlen >= BYTES(MRO_B)) {
    __m256i H[4];
    V1_ALPHA_UPDATE_1(L);

    V1_LOAD_BLOCK(H, h);
    V1_BLOCKCIPHER_F(H, L);
    V1_ACCUMULATE(B, H);

    V1_ALPHA_UPDATE_2(L);
    h    += BYTES(MRO_B);
    hlen -= BYTES(MRO_B);
  }

  if(hlen > 0) {
    uint8_t lastblock[BYTES(MRO_B)];
    __m256i H[4];
    MRO_PAD(lastblock, h, hlen);
    V1_ALPHA_UPDATE_1(L);
    V1_LOAD_BLOCK(H, lastblock);
    V1_BLOCKCIPHER_F(H, L);
    V1_ACCUMULATE(B, H);
    V1_ALPHA_UPDATE_2(L);
  }

  /* Absorb Message */
  memcpy(L, L_, BYTES(MRO_B));
  V1_BETA_UPDATE(L);
  while(mlen >= 4 * BYTES(MRO_B)) {
    __m256i M[16];
    V4_ALPHA_UPDATE_1(L);
    V4_LOAD_BLOCK(M, m);
    V4_BLOCKCIPHER_F(M, L);
    V4_ACCUMULATE(B, M);
    V4_ALPHA_UPDATE_2(L);
    m    += 4 * BYTES(MRO_B);
    mlen -= 4 * BYTES(MRO_B);
  }

  /* TODO: V2 */

  while(mlen >= BYTES(MRO_B)) {
    __m256i M[4];
    V1_ALPHA_UPDATE_1(L);

    V1_LOAD_BLOCK(M, m);
    V1_BLOCKCIPHER_F(M, L);
    V1_ACCUMULATE(B, M);

    V1_ALPHA_UPDATE_2(L);
    m    += BYTES(MRO_B);
    mlen -= BYTES(MRO_B);
  }

  if(mlen > 0) {
    uint8_t lastblock[BYTES(MRO_B)];
    __m256i M[4];
    MRO_PAD(lastblock, m, mlen);
    V1_ALPHA_UPDATE_1(L);
    V1_LOAD_BLOCK(M, lastblock);
    V1_BLOCKCIPHER_F(M, L);
    V1_ACCUMULATE(B, M);
    V1_ALPHA_UPDATE_2(L);
  }

  {
    memcpy(L, L_, BYTES(MRO_B));
    V1_BETA2_UPDATE(L);
    B[3] = XOR256(B[3], _mm256_set_epi64x(mlen_, hlen_, 0, 0));
    V1_BLOCKCIPHER_F(B, L);
  }
  memcpy(tag, B, BYTES(MRO_T));
}


static void mro_encrypt_data(
  uint8_t * c,
  const unsigned char * m,
  size_t mlen,
  uint64_t L[16+4],
  uint64_t T[4])
{
  __m256i counter = _mm256_set_epi64x(3, 2, 1, 0);
  __m256i K[4];
  __m256i K_[16];
  int i;
  const __m256i T0 = _mm256_set1_epi64x(T[0]);
  const __m256i T1 = _mm256_set1_epi64x(T[1]);
  const __m256i T2 = _mm256_set1_epi64x(T[2]);
  const __m256i T3 = _mm256_set1_epi64x(T[3]);

  memcpy(K, L, 16 * BYTES(MRO_W));

  if(mlen >= 4 * BYTES(MRO_B)) {
    int i;
    for(i = 0; i < 16; ++i) K_[i] = _mm256_set1_epi64x(L[i]);
/*
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
    */
  }

  while(mlen >= 4 * BYTES(MRO_B)) {
    __m256i B[16], C[16];
    int i;
    memcpy(C, K_, 4 * BYTES(MRO_B));

    C[ 0] = XOR256(C[ 0], T0);
    C[ 1] = XOR256(C[ 1], T1);
    C[ 2] = XOR256(C[ 2], T2);
    C[ 3] = XOR256(C[ 3], T3);
    C[15] = XOR256(C[15], counter);

    V4_PERMUTE_F(C);

    for(i = 0; i < 16; ++i)
      C[i] = XOR256(C[i], K_[i]);
    V4_TRANSPOSE_B(B, C);
    for(i = 0; i < 16; ++i)
      STOREU256(&c[i * 32], XOR256(B[i], LOADU256(&m[i * 32])));

    counter = ADD256(counter, _mm256_set_epi64x(4, 4, 4, 4));
    c    += 4 * BYTES(MRO_B);
    m    += 4 * BYTES(MRO_B);
    mlen -= 4 * BYTES(MRO_B);
  }

  counter = _mm256_set_epi64x(_mm256_extract_epi64(counter, 0), 0, 0, 0);
  while(mlen >= BYTES(MRO_B)) {
    __m256i B[4];

    V1_COPY_BLOCK(B, K);
    B[0] = XOR256(B[0], LOADU256(T));
    B[3] = XOR256(B[3], counter);
    V1_PERMUTE_F(B);
    for(i = 0; i < 4; ++i) B[i] = XOR256(B[i], K[i]);

    for(i = 0; i < 4; ++i)
      STOREU256(&c[32 * i], XOR256(B[i], LOADU256(&m[32 * i])));

    counter = ADD256(counter, _mm256_set_epi64x(1, 0, 0, 0));
    c    += BYTES(MRO_B);
    m    += BYTES(MRO_B);
    mlen -= BYTES(MRO_B);
  }

  if(mlen > 0) {
    uint8_t lastblock[BYTES(MRO_B)];
    __m256i B[4];
    memcpy(lastblock, m, mlen);
    for(i = 0; i < 4; ++i) B[i] = K[i];
    B[0] = XOR256(B[0], LOADU256(T));
    B[3] = XOR256(B[3], counter);
    V1_PERMUTE_F(B);
    for(i = 0; i < 4; ++i) B[i] = XOR256(B[i], K[i]);
    for(i = 0; i < 4; ++i)
      B[i] = XOR256(B[i], LOADU256(&lastblock[32 * i]));
    memcpy(c, B, mlen);
  }
}

static void mro_decrypt_data(
  uint8_t * m,
  const unsigned char * c,
  size_t clen,
  uint64_t L[16+4],
  uint64_t T[4])
{
  mro_encrypt_data(m, c, clen, L, T);
}

/* high level interface functions */
void mro_aead_encrypt(
    unsigned char *c, size_t *clen,
    const unsigned char *h, size_t hlen,
    const unsigned char *m, size_t mlen,
    const unsigned char *n,
    const unsigned char *k
    )
{
  uint64_t La[16+4], Le[16+4];
  uint64_t t[4];
  mro_kdf(La, Le, k, n);
  mro_compute_tag(t, h, hlen, m, mlen, La);
  mro_encrypt_data(c, m, mlen, Le, t);
  memcpy(c + mlen, t, BYTES(MRO_T));
  *clen = mlen + BYTES(MRO_T);
}

int mro_aead_decrypt(
    unsigned char *m, size_t *mlen,
    const unsigned char *h, size_t hlen,
    const unsigned char *c, size_t clen,
    const unsigned char *n,
    const unsigned char *k)
{
  uint64_t La[16+4], Le[16+4];
  uint64_t t[4];
  __m256i v;

  if (clen < BYTES(MRO_T))
    return -1;

  mro_kdf(La, Le, k, n);

  *mlen = clen - BYTES(MRO_T);
  mro_decrypt_data(m, c, *mlen, Le, (uint64_t *)(c + *mlen)); /* unaligned */
  mro_compute_tag(t, h, hlen, m, *mlen, La);

  v = _mm256_cmpeq_epi8(LOADU256(t), LOADU256(c + clen - BYTES(MRO_T)));
  return (( (_mm256_movemask_epi8(v) & 0xFFFFFFFFULL) + 1) >> 32) - 1;
}

