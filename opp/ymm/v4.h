#ifndef STORM_OPP_YMM_V4_H
#define STORM_OPP_YMM_V4_H

#include "v0.h"

#define V4_G_F(A, B, C, D) do {                          \
  A = ADD256(A, B); D = XOR256(D, A); D = ROT256(D, 32); \
  C = ADD256(C, D); B = XOR256(B, C); B = ROT256(B, 24); \
  A = ADD256(A, B); D = XOR256(D, A); D = ROT256(D, 16); \
  C = ADD256(C, D); B = XOR256(B, C); B = ROT256(B, 63); \
} while(0)

#define V4_G_B(A, B, C, D) do {                          \
  B = ROT256(B,  1); B = XOR256(B, C); C = SUB256(C, D); \
  D = ROT256(D, 48); D = XOR256(D, A); A = SUB256(A, B); \
  B = ROT256(B, 40); B = XOR256(B, C); C = SUB256(C, D); \
  D = ROT256(D, 32); D = XOR256(D, A); A = SUB256(A, B); \
} while(0)

#define V4_PERMUTE_F(B) do {            \
  int i;                                \
  for(i = 0; i < STORM_R; ++i) {        \
    /* Column step */                   \
    V4_G_F(B[ 0], B[ 4], B[ 8], B[12]); \
    V4_G_F(B[ 1], B[ 5], B[ 9], B[13]); \
    V4_G_F(B[ 2], B[ 6], B[10], B[14]); \
    V4_G_F(B[ 3], B[ 7], B[11], B[15]); \
    /* Diagonal step */                 \
    V4_G_F(B[ 0], B[ 5], B[10], B[15]); \
    V4_G_F(B[ 1], B[ 6], B[11], B[12]); \
    V4_G_F(B[ 2], B[ 7], B[ 8], B[13]); \
    V4_G_F(B[ 3], B[ 4], B[ 9], B[14]); \
  }                                     \
} while(0)

#define V4_PERMUTE_B(B) do {            \
  int i;                                \
  for(i = 0; i < STORM_R; ++i) {        \
    /* Diagonal step */                 \
    V4_G_B(B[ 0], B[ 5], B[10], B[15]); \
    V4_G_B(B[ 1], B[ 6], B[11], B[12]); \
    V4_G_B(B[ 2], B[ 7], B[ 8], B[13]); \
    V4_G_B(B[ 3], B[ 4], B[ 9], B[14]); \
      /* Column step */                 \
    V4_G_B(B[ 0], B[ 4], B[ 8], B[12]); \
    V4_G_B(B[ 1], B[ 5], B[ 9], B[13]); \
    V4_G_B(B[ 2], B[ 6], B[10], B[14]); \
    V4_G_B(B[ 3], B[ 7], B[11], B[15]); \
  }                                     \
} while(0)

#define V4_TRANSPOSE_F(C, B) do {                   \
  __m256i t0, t1, t2, t3;                           \
  t0 = _mm256_unpacklo_epi64(B[ 0], B[ 4]);         \
  t1 = _mm256_unpackhi_epi64(B[ 0], B[ 4]);         \
  t2 = _mm256_unpacklo_epi64(B[ 8], B[12]);         \
  t3 = _mm256_unpackhi_epi64(B[ 8], B[12]);         \
  C[ 0] = _mm256_permute2x128_si256(t0, t2, 0x20);  \
  C[ 2] = _mm256_permute2x128_si256(t0, t2, 0x31);  \
  C[ 1] = _mm256_permute2x128_si256(t1, t3, 0x20);  \
  C[ 3] = _mm256_permute2x128_si256(t1, t3, 0x31);  \
  t0 = _mm256_unpacklo_epi64(B[ 1], B[ 5]);         \
  t1 = _mm256_unpackhi_epi64(B[ 1], B[ 5]);         \
  t2 = _mm256_unpacklo_epi64(B[ 9], B[13]);         \
  t3 = _mm256_unpackhi_epi64(B[ 9], B[13]);         \
  C[ 4] = _mm256_permute2x128_si256(t0, t2, 0x20);  \
  C[ 6] = _mm256_permute2x128_si256(t0, t2, 0x31);  \
  C[ 5] = _mm256_permute2x128_si256(t1, t3, 0x20);  \
  C[ 7] = _mm256_permute2x128_si256(t1, t3, 0x31);  \
  t0 = _mm256_unpacklo_epi64(B[ 2], B[ 6]);         \
  t1 = _mm256_unpackhi_epi64(B[ 2], B[ 6]);         \
  t2 = _mm256_unpacklo_epi64(B[10], B[14]);         \
  t3 = _mm256_unpackhi_epi64(B[10], B[14]);         \
  C[ 8] = _mm256_permute2x128_si256(t0, t2, 0x20);  \
  C[10] = _mm256_permute2x128_si256(t0, t2, 0x31);  \
  C[ 9] = _mm256_permute2x128_si256(t1, t3, 0x20);  \
  C[11] = _mm256_permute2x128_si256(t1, t3, 0x31);  \
  t0 = _mm256_unpacklo_epi64(B[ 3], B[ 7]);         \
  t1 = _mm256_unpackhi_epi64(B[ 3], B[ 7]);         \
  t2 = _mm256_unpacklo_epi64(B[11], B[15]);         \
  t3 = _mm256_unpackhi_epi64(B[11], B[15]);         \
  C[12] = _mm256_permute2x128_si256(t0, t2, 0x20);  \
  C[14] = _mm256_permute2x128_si256(t0, t2, 0x31);  \
  C[13] = _mm256_permute2x128_si256(t1, t3, 0x20);  \
  C[15] = _mm256_permute2x128_si256(t1, t3, 0x31);  \
} while(0)


#define V4_TRANSPOSE_B(C, B) do {                   \
  __m256i t0, t1, t2, t3;                           \
  t0 = _mm256_unpacklo_epi64(B[ 0], B[ 1]);         \
  t1 = _mm256_unpackhi_epi64(B[ 0], B[ 1]);         \
  t2 = _mm256_unpacklo_epi64(B[ 2], B[ 3]);         \
  t3 = _mm256_unpackhi_epi64(B[ 2], B[ 3]);         \
  C[ 0] = _mm256_permute2x128_si256(t0, t2, 0x20);  \
  C[ 8] = _mm256_permute2x128_si256(t0, t2, 0x31);  \
  C[ 4] = _mm256_permute2x128_si256(t1, t3, 0x20);  \
  C[12] = _mm256_permute2x128_si256(t1, t3, 0x31);  \
  t0 = _mm256_unpacklo_epi64(B[ 4], B[ 5]);         \
  t1 = _mm256_unpackhi_epi64(B[ 4], B[ 5]);         \
  t2 = _mm256_unpacklo_epi64(B[ 6], B[ 7]);         \
  t3 = _mm256_unpackhi_epi64(B[ 6], B[ 7]);         \
  C[ 1] = _mm256_permute2x128_si256(t0, t2, 0x20);  \
  C[ 9] = _mm256_permute2x128_si256(t0, t2, 0x31);  \
  C[ 5] = _mm256_permute2x128_si256(t1, t3, 0x20);  \
  C[13] = _mm256_permute2x128_si256(t1, t3, 0x31);  \
  t0 = _mm256_unpacklo_epi64(B[ 8], B[ 9]);         \
  t1 = _mm256_unpackhi_epi64(B[ 8], B[ 9]);         \
  t2 = _mm256_unpacklo_epi64(B[10], B[11]);         \
  t3 = _mm256_unpackhi_epi64(B[10], B[11]);         \
  C[ 2] = _mm256_permute2x128_si256(t0, t2, 0x20);  \
  C[10] = _mm256_permute2x128_si256(t0, t2, 0x31);  \
  C[ 6] = _mm256_permute2x128_si256(t1, t3, 0x20);  \
  C[14] = _mm256_permute2x128_si256(t1, t3, 0x31);  \
  t0 = _mm256_unpacklo_epi64(B[12], B[13]);         \
  t1 = _mm256_unpackhi_epi64(B[12], B[13]);         \
  t2 = _mm256_unpacklo_epi64(B[14], B[15]);         \
  t3 = _mm256_unpackhi_epi64(B[14], B[15]);         \
  C[ 3] = _mm256_permute2x128_si256(t0, t2, 0x20);  \
  C[11] = _mm256_permute2x128_si256(t0, t2, 0x31);  \
  C[ 7] = _mm256_permute2x128_si256(t1, t3, 0x20);  \
  C[15] = _mm256_permute2x128_si256(t1, t3, 0x31);  \
} while(0)

#define V4_XOR_MASK(B, L) do {              \
  int i;                                    \
  for(i = 0; i < 16; ++i) {                 \
    B[i] = XOR256(B[i], LOADU256(&L[i+0])); \
  }                                         \
} while(0)

#define V4_LOAD_BLOCK(B, m) do { \
  int i;                         \
  for(i = 0; i < 16; ++i) {      \
    B[i] = LOADU256(&m[32 * i]); \
  }                              \
} while(0)

#define V4_STORE_BLOCK(c, B) do {\
  int i;                         \
  for(i = 0; i < 16; ++i) {      \
    STOREU256(&c[32 * i], B[i]); \
  }                              \
} while(0)

#define V4_ACCUMULATE(T, B) do {    \
  int i;                            \
  for(i = 0; i < 4; ++i) {          \
    T[i] = XOR256(T[i], B[i +  0]); \
    T[i] = XOR256(T[i], B[i +  4]); \
    T[i] = XOR256(T[i], B[i +  8]); \
    T[i] = XOR256(T[i], B[i + 12]); \
  }                                 \
} while(0)

#define V4_PHI_UPDATE_1(L) do {                           \
  STOREU256(&L[16], XOR256(ROT256(LOADU256(&L[0]), 11),   \
                           SHL256(LOADU256(&L[5]), 13))); \
} while(0)

#define V4_PHI_UPDATE_2(L) do {        \
  STOREU256(&L[ 0], LOADU256(&L[ 4])); \
  STOREU256(&L[ 4], LOADU256(&L[ 8])); \
  STOREU256(&L[ 8], LOADU256(&L[12])); \
  STOREU256(&L[12], LOADU256(&L[16])); \
} while(0)

#define V4_PHI_UPDATE(L) do {                             \
  STOREU256(&L[ 0], LOADU256(&L[ 4]));                    \
  STOREU256(&L[ 4], LOADU256(&L[ 8]));                    \
  STOREU256(&L[ 8], LOADU256(&L[12]));                    \
  STOREU256(&L[12], LOADU256(&L[16]));                    \
  STOREU256(&L[16], XOR256(ROT256(LOADU256(&L[0]), 11),   \
                           SHL256(LOADU256(&L[5]), 13))); \
} while(0)


#define V4_BLOCKCIPHER_F(B, L) do { \
  __m256i C[16];                    \
  V4_TRANSPOSE_F(C, B);             \
  V4_XOR_MASK(C, L);                \
  V4_PERMUTE_F(C);                  \
  V4_XOR_MASK(C, L);                \
  V4_TRANSPOSE_B(B, C);             \
} while(0)

#define V4_BLOCKCIPHER_B(B, L) do { \
  __m256i C[16];                    \
  V4_TRANSPOSE_F(C, B);             \
  V4_XOR_MASK(C, L);                \
  V4_PERMUTE_B(C);                  \
  V4_XOR_MASK(C, L);                \
  V4_TRANSPOSE_B(B, C);             \
} while(0)

#endif

