#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <immintrin.h>

#define OPP_W 64
#define OPP_R 4
#define OPP_T (OPP_W * 4)
#define OPP_B 1024

#include "v0.h"
#include "v1.h"
#include "v2.h"
#include "v4.h"

/* this is x86, so we can just memcpy */
/*
static uint64_t load64(const void * in) {
  uint64_t x;
  memcpy(&x, in, sizeof x);
  return x;
}
*/

static void opp_pad(unsigned char * out, const void * in, size_t inlen) {
  memset(out, 0, BYTES(OPP_B));
  memcpy(out, in, inlen);
  out[inlen] = 0x01;
}

static void opp_kdf(uint64_t * Ka, uint64_t * Ke, const uint8_t * k, const uint8_t * n) {
  __m256i B[4];
  B[0] = _mm256_castsi128_si256(LOADU128(n));
  B[1] = _mm256_setzero_si256();
  B[2] = _mm256_set_epi64x(OPP_T, OPP_R, 0, 0);
  B[3] = LOADU256(k);

  V1_PERMUTE_F(B);

  memcpy(Ka, B, sizeof(B));
  memcpy(Ke, B, sizeof(B));
  V1_GAMMA_UPDATE(Ke);
}

static void opp_hash_data(__m256i T[4], const uint8_t * h, size_t hlen, uint64_t L[16+4]) {
  while(hlen >= 4 * BYTES(OPP_B)) {
    __m256i B[16];

    V4_ALPHA_UPDATE_1(L);
    V4_LOAD_BLOCK(B, h);
    V4_BLOCKCIPHER_F(B, L);
    V4_ACCUMULATE(T, B);
    V4_ALPHA_UPDATE_2(L);
    h    += 4 * BYTES(OPP_B);
    hlen -= 4 * BYTES(OPP_B);
  }

  /* TODO: V2 */

  while(hlen >= BYTES(OPP_B)) {
    __m256i B[4];

    V1_LOAD_BLOCK(B, h);
    V1_BLOCKCIPHER_F(B, L);
    V1_ACCUMULATE(T, B);

    V1_ALPHA_UPDATE(L);
    h    += BYTES(OPP_B);
    hlen -= BYTES(OPP_B);
  }

  if(hlen > 0) {
    uint8_t lastblock[BYTES(OPP_B)];
    __m256i B[4];
    V1_BETA_UPDATE(L);
    opp_pad(lastblock, h, hlen);
    V1_LOAD_BLOCK(B, lastblock);
    V1_BLOCKCIPHER_F(B, L);
    V1_ACCUMULATE(T, B);
  }
}


static void opp_encrypt_data(__m256i T[4], uint8_t * c, const uint8_t * m, size_t mlen, uint64_t L[16+4]) {
  while(mlen >= 4 * BYTES(OPP_B)) {
    __m256i B[16];
    V4_ALPHA_UPDATE_1(L);
    V4_LOAD_BLOCK(B, m);
    V4_ACCUMULATE(T, B);
    V4_BLOCKCIPHER_F(B, L);
    V4_STORE_BLOCK(c, B);
    V4_ALPHA_UPDATE_2(L);
    c    += 4 * BYTES(OPP_B);
    m    += 4 * BYTES(OPP_B);
    mlen -= 4 * BYTES(OPP_B);
  }

  /* TODO: V2 */

  while(mlen >= BYTES(OPP_B)) {
    __m256i B[4];

    V1_LOAD_BLOCK(B, m);
    V1_ACCUMULATE(T, B);
    V1_BLOCKCIPHER_F(B, L);
    V1_STORE_BLOCK(c, B);

    V1_ALPHA_UPDATE(L);
    c    += BYTES(OPP_B);
    m    += BYTES(OPP_B);
    mlen -= BYTES(OPP_B);
  }

  if(mlen > 0) { /* handle partial final block */
    uint8_t lastblock[BYTES(OPP_B)];
    __m256i B[4];
    int i;
    V1_BETA_UPDATE(L);
    opp_pad(lastblock, m, mlen);
    V1_ZERO_BLOCK(B);
    V1_BLOCKCIPHER_F(B, L);
    for(i = 0; i < 4; ++i) { /* lastblock xor B and T xor last block */
      const __m256i M_i = LOADU256(&lastblock[32 * i]);
      T[i] = XOR256(T[i], M_i);
      STOREU256(&lastblock[32 * i], XOR256(B[i], M_i));
    }
    memcpy(c, lastblock, mlen);
  }
}

static void opp_decrypt_data(__m256i T[4], uint8_t * m, const uint8_t * c, size_t clen, uint64_t L[16+4]) {
  while(clen >= 4 * BYTES(OPP_B)) {
    __m256i B[16];
    V4_ALPHA_UPDATE_1(L);
    V4_LOAD_BLOCK(B, c);
    V4_BLOCKCIPHER_B(B, L);
    V4_ACCUMULATE(T, B);
    V4_STORE_BLOCK(m, B);
    V4_ALPHA_UPDATE_2(L);
    m    += 4 * BYTES(OPP_B);
    c    += 4 * BYTES(OPP_B);
    clen -= 4 * BYTES(OPP_B);
  }

  /* TODO: V2 */

  while(clen >= BYTES(OPP_B)) {
    __m256i B[4];

    V1_LOAD_BLOCK(B, c);
    V1_BLOCKCIPHER_B(B, L);
    V1_ACCUMULATE(T, B);
    V1_STORE_BLOCK(m, B);

    V1_ALPHA_UPDATE(L);
    m    += BYTES(OPP_B);
    c    += BYTES(OPP_B);
    clen -= BYTES(OPP_B);
  }

  if(clen > 0) { /* handle partial final block */
    uint8_t lastblock[BYTES(OPP_B)];
    __m256i B[4];
    int i;
    V1_BETA_UPDATE(L);
    opp_pad(lastblock, c, clen);
    V1_ZERO_BLOCK(B);
    V1_BLOCKCIPHER_F(B, L);
    for(i = 0; i < 4; ++i) { /* lastblock xor B */
      const __m256i C_i = LOADU256(&lastblock[32 * i]);
      STOREU256(&lastblock[32 * i], XOR256(B[i], C_i));
    }
    memcpy(m, lastblock, clen);
    opp_pad(lastblock, m, clen);
    for(i = 0; i < 4; ++i) { /* T xor last block */
      T[i] = XOR256(T[i], LOADU256(&lastblock[32 * i]));
    }
  }
}

static void opp_tag(__m256i * Te, const __m256i * Ta, uint64_t * L, size_t hlen, size_t mlen) {
  const int mode =   (hlen % BYTES(OPP_B) != 0) +
                   2*(mlen % BYTES(OPP_B) != 0);
  switch(mode) {
  case 2: /* hlen == 128, mlen != 128: i2 = 0 -> 3 */
    V1_BETA_UPDATE(L);
  case 0: /* hlen == 128, mlen == 128: i2 = 0 -> 2 */
  case 3: /* hlen != 128, mlen != 128: i2 = 1 -> 3 */
    V1_BETA2_UPDATE(L);
    break;
  case 1: /* hlen != 128, mlen == 128: i2 = 1 -> 2 */
    V1_BETA_UPDATE(L);
    break;
  }
  V1_BLOCKCIPHER_F(Te, L);
  V1_ACCUMULATE(Te, Ta);
}

#if defined(OPP_DEBUG)
static void print_mask(uint64_t * L) {
  int i;
  for(i = 0; i < 16; ++i) {
    printf("%016lX%c", L[i], i % 4 == 3 ? '\n' : ' ');
  }
  printf("\n");
}

static void print_state(__m256i * B) {
  uint64_t L[16];
  int i;
  for(i = 0; i < 4; ++i) {
    STOREU256(&L[4*i], B[i]);
  }
  print_mask(L);
}
#endif

/* high level interface functions */
void crypto_aead_encrypt(
    unsigned char *c, size_t *clen,
    const unsigned char *h, size_t hlen,
    const unsigned char *m, size_t mlen,
    const unsigned char *n,
    const unsigned char *k
    )
{
  __m256i Ta[4] = {0};
  __m256i Te[4] = {0};
  uint64_t Ka[16+4];
  uint64_t Ke[16+4];

  opp_kdf(Ka, Ke, k, n);

#if defined(OPP_DEBUG)
  print_mask(Ka);
  print_mask(Ke);
#endif

  opp_hash_data(Ta, h, hlen, Ka);
  opp_encrypt_data(Te, c, m, mlen, Ke);
  opp_tag(Te, Ta, Ka, hlen, mlen);

#if defined(OPP_DEBUG)
  print_state(Te);
#endif

  *clen = mlen + BYTES(OPP_T);
  STOREU256(c + mlen, Te[0]);

#if defined(DEBUG)
  {
    int i;
    for(i = 0; i < *clen; ++i)
      printf("%02X ", c[i]);
    printf("\n");
  }
#endif
}

int crypto_aead_decrypt(
    unsigned char *m, size_t *mlen,
    const unsigned char *h, size_t hlen,
    const unsigned char *c, size_t clen,
    const unsigned char *n,
    const unsigned char *k)
{
  __m256i Ta[4] = {0};
  __m256i Te[4] = {0};
  uint64_t Ka[16+4];
  uint64_t Ke[16+4];

  if (clen < BYTES(OPP_T))
    return -1;
  *mlen = clen - BYTES(OPP_T);

  opp_kdf(Ka, Ke, k, n);

  opp_hash_data(Ta, h, hlen, Ka);
  opp_decrypt_data(Te, m, c, clen - BYTES(OPP_T), Ke);
  opp_tag(Te, Ta, Ka, hlen, *mlen);

  Te[0] = _mm256_cmpeq_epi8(Te[0], LOADU256(c + clen - BYTES(OPP_T)));
  return (( (_mm256_movemask_epi8(Te[0]) & 0xFFFFFFFFULL) + 1) >> 32) - 1;
}

