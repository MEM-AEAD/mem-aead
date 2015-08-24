#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <immintrin.h>

#define STORM_W 64
#define STORM_R 4
#define STORM_T (STORM_W * 4)
#define STORM_B 1024

#include "v0.h"
#include "v1.h"
#include "v2.h"
#include "v4.h"

#define STORM_PAD(out, in, inlen) do { \
  memset(out, 0, BYTES(STORM_B));      \
  memcpy(out, in, inlen);              \
  out[inlen] = 0x01;                   \
} while(0)

static void storm_kdf(uint64_t * Ka, uint64_t * Ke, const uint8_t * k, const uint8_t * n) {
  __m128i B[16] = {0};

  B[ 0] = _mm_set1_epi64x(load64(n + 0));
  B[ 1] = _mm_set1_epi64x(load64(n + 8));

  B[ 4] = _mm_set1_epi64x(load64(k +  0));
  B[ 5] = _mm_set1_epi64x(load64(k +  8));
  B[ 6] = _mm_set1_epi64x(load64(k + 16));
  B[ 7] = _mm_set1_epi64x(load64(k + 24));

  B[12] = _mm_set1_epi64x(STORM_W);
  B[13] = _mm_set1_epi64x(STORM_R);
  B[14] = _mm_set1_epi64x(STORM_T);
  B[15] = _mm_set_epi64x(1, 0);

  V2_PERMUTE_F(B);

#define V2_UNPACK(i) do { \
  __m128i t0, t1, t2, t3; \
  t0 = _mm_unpacklo_epi64(B[4*i+0], B[4*i+1]); \
  t1 = _mm_unpackhi_epi64(B[4*i+0], B[4*i+1]); \
  t2 = _mm_unpacklo_epi64(B[4*i+2], B[4*i+3]); \
  t3 = _mm_unpackhi_epi64(B[4*i+2], B[4*i+3]); \
  STOREU256(&Ka[i*4], _mm256_permute2x128_si256(_mm256_castsi128_si256(t0), _mm256_castsi128_si256(t2), 0x20)); \
  STOREU256(&Ke[i*4], _mm256_permute2x128_si256(_mm256_castsi128_si256(t1), _mm256_castsi128_si256(t3), 0x20)); \
} while(0)
  V2_UNPACK(0);
  V2_UNPACK(1);
  V2_UNPACK(2);
  V2_UNPACK(3);
#undef V2_UNPACK
}

static void storm_hash_data(__m256i T[4], const uint8_t * h, size_t hlen, uint64_t L[16+4]) {
  while(hlen >= 4 * BYTES(STORM_B)) {
    __m256i B[16];

    V4_MASK_UPDATE_1(L);
    V4_LOAD_BLOCK(B, h);
    V4_BLOCKCIPHER_F(B, L);
    V4_ACCUMULATE(T, B);
    V4_MASK_UPDATE_2(L);
    h    += 4 * BYTES(STORM_B);
    hlen -= 4 * BYTES(STORM_B);
  }

  /* TODO: V2 */

  while(hlen >= BYTES(STORM_B)) {
    __m256i B[4];
    V1_MASK_UPDATE(L);

    V1_LOAD_BLOCK(B, h);
    V1_BLOCKCIPHER_F(B, L);
    V1_ACCUMULATE(T, B);

    h    += BYTES(STORM_B);
    hlen -= BYTES(STORM_B);
  }

  if(hlen > 0) {
    uint8_t lastblock[BYTES(STORM_B)];
    __m256i B[4];
    V1_MASK_UPDATE(L);
    STORM_PAD(lastblock, h, hlen);
    V1_LOAD_BLOCK(B, lastblock);
    V1_BLOCKCIPHER_ROTATED_F(B, L, 256);
    V1_ACCUMULATE(T, B);
  }
}


static void storm_encrypt_data(__m256i T[4], uint8_t * c, const uint8_t * m, size_t mlen, uint64_t L[16+4]) {
  while(mlen >= 4 * BYTES(STORM_B)) {
    __m256i B[16];
    V4_MASK_UPDATE_1(L);
    V4_LOAD_BLOCK(B, m);
    V4_ACCUMULATE(T, B);
    V4_BLOCKCIPHER_F(B, L);
    V4_STORE_BLOCK(c, B);
    V4_MASK_UPDATE_2(L);
    c    += 4 * BYTES(STORM_B);
    m    += 4 * BYTES(STORM_B);
    mlen -= 4 * BYTES(STORM_B);
  }

  /* TODO: V2 */

  while(mlen >= BYTES(STORM_B)) {
    __m256i B[4];
    V1_MASK_UPDATE(L);

    V1_LOAD_BLOCK(B, m);
    V1_ACCUMULATE(T, B);
    V1_BLOCKCIPHER_F(B, L);
    V1_STORE_BLOCK(c, B);

    c    += BYTES(STORM_B);
    m    += BYTES(STORM_B);
    mlen -= BYTES(STORM_B);
  }

  if(mlen > 0) { /* handle partial final block */
    uint8_t lastblock[BYTES(STORM_B)];
    __m256i B[4];
    int i;
    V1_MASK_UPDATE(L);
    STORM_PAD(lastblock, m, mlen);
    V1_ZERO_BLOCK(B);
    V1_BLOCKCIPHER_ROTATED_F(B, L, 768);
    for(i = 0; i < 4; ++i) { /* lastblock xor B and T xor last block */
      const __m256i M_i = LOADU256(&lastblock[32 * i]);
      T[i] = XOR256(T[i], M_i);
      STOREU256(&lastblock[32 * i], XOR256(B[i], M_i));
    }
    memcpy(c, lastblock, mlen);
  }
}

static void storm_decrypt_data(__m256i T[4], uint8_t * m, const uint8_t * c, size_t clen, uint64_t L[16+4]) {
  while(clen >= 4 * BYTES(STORM_B)) {
    __m256i B[16];
    V4_MASK_UPDATE_1(L);
    V4_LOAD_BLOCK(B, c);
    V4_BLOCKCIPHER_B(B, L);
    V4_ACCUMULATE(T, B);
    V4_STORE_BLOCK(m, B);
    V4_MASK_UPDATE_2(L);
    m    += 4 * BYTES(STORM_B);
    c    += 4 * BYTES(STORM_B);
    clen -= 4 * BYTES(STORM_B);
  }

  /* TODO: V2 */

  while(clen >= BYTES(STORM_B)) {
    __m256i B[4];
    V1_MASK_UPDATE(L);

    V1_LOAD_BLOCK(B, c);
    V1_BLOCKCIPHER_B(B, L);
    V1_ACCUMULATE(T, B);
    V1_STORE_BLOCK(m, B);

    m    += BYTES(STORM_B);
    c    += BYTES(STORM_B);
    clen -= BYTES(STORM_B);
  }

  if(clen > 0) { /* handle partial final block */
    uint8_t lastblock[BYTES(STORM_B)];
    __m256i B[4];
    int i;

    V1_MASK_UPDATE(L);

    STORM_PAD(lastblock, c, clen);
    V1_ZERO_BLOCK(B);
    V1_BLOCKCIPHER_ROTATED_F(B, L, 768);
    for(i = 0; i < 4; ++i) { /* lastblock xor B */
      const __m256i C_i = LOADU256(&lastblock[32 * i]);
      STOREU256(&lastblock[32 * i], XOR256(B[i], C_i));
    }
    memcpy(m, lastblock, clen);
    STORM_PAD(lastblock, m, clen);
    for(i = 0; i < 4; ++i) { /* T xor last block */
      T[i] = XOR256(T[i], LOADU256(&lastblock[32 * i]));
    }
  }
}

static void storm_tag(__m256i * Te, const __m256i * Ta, const uint64_t * L) {
  V1_BLOCKCIPHER_ROTATED_F(Te, L, 512);
  V1_ACCUMULATE(Te, Ta);
}

#if defined(STORM_DEBUG)
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
void storm_aead_encrypt(
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

  storm_kdf(Ka, Ke, k, n);

#if defined(STORM_DEBUG)
  print_mask(Ka);
  print_mask(Ke);
#endif

  storm_hash_data(Ta, h, hlen, Ka);
  storm_encrypt_data(Te, c, m, mlen, Ke);
  storm_tag(Te, Ta, Ka);

#if defined(STORM_DEBUG)
  print_state(Te);
#endif

  *clen = mlen + BYTES(STORM_T);
  STOREU256(c + mlen, Te[0]);

#if defined(DEBUG)
  {
    int i;
    for(i = 0; i < 32; ++i)
      printf("%02X ", c[mlen + i]);
    printf("\n");
  }
#endif
}

int storm_aead_decrypt(
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

  if (clen < BYTES(STORM_T))
    return -1;

  storm_kdf(Ka, Ke, k, n);

  storm_hash_data(Ta, h, hlen, Ka);
  storm_decrypt_data(Te, m, c, clen - BYTES(STORM_T), Ke);
  storm_tag(Te, Ta, Ka);

  *mlen = clen - BYTES(STORM_T);
  Te[0] = _mm256_cmpeq_epi8(Te[0], LOADU256(c + clen - BYTES(STORM_T)));
  return (( (_mm256_movemask_epi8(Te[0]) & 0xFFFFFFFFULL) + 1) >> 32) - 1;
}

