#ifndef MRO_YMM_V1_H
#define MRO_YMM_V1_H

#include "v0.h"

#define V1_DIAG_F(A, B, C, D) do {                       \
  D = _mm256_permute4x64_epi64(D, _MM_SHUFFLE(2,1,0,3)); \
  C = _mm256_permute4x64_epi64(C, _MM_SHUFFLE(1,0,3,2)); \
  B = _mm256_permute4x64_epi64(B, _MM_SHUFFLE(0,3,2,1)); \
} while(0)

#define V1_DIAG_B(A, B, C, D) do {                       \
  D = _mm256_permute4x64_epi64(D, _MM_SHUFFLE(0,3,2,1)); \
  C = _mm256_permute4x64_epi64(C, _MM_SHUFFLE(1,0,3,2)); \
  B = _mm256_permute4x64_epi64(B, _MM_SHUFFLE(2,1,0,3)); \
} while(0)

#define V1_G_F(A, B, C, D) do {                          \
  A = ADD256(A, B); D = XOR256(D, A); D = ROT256(D, 32); \
  C = ADD256(C, D); B = XOR256(B, C); B = ROT256(B, 24); \
  A = ADD256(A, B); D = XOR256(D, A); D = ROT256(D, 16); \
  C = ADD256(C, D); B = XOR256(B, C); B = ROT256(B, 63); \
} while(0);

#define V1_G_B(A, B, C, D) do {                          \
  B = ROT256(B,  1); B = XOR256(B, C); C = SUB256(C, D); \
  D = ROT256(D, 48); D = XOR256(D, A); A = SUB256(A, B); \
  B = ROT256(B, 40); B = XOR256(B, C); C = SUB256(C, D); \
  D = ROT256(D, 32); D = XOR256(D, A); A = SUB256(A, B); \
} while(0)

#define V1_PERMUTE_F(B) do {           \
  int i;                               \
  for(i = 0; i < MRO_L; ++i) {       \
    V1_G_F(B[0], B[1], B[2], B[3]);    \
    V1_DIAG_F(B[0], B[1], B[2], B[3]); \
    V1_G_F(B[0], B[1], B[2], B[3]);    \
    V1_DIAG_B(B[0], B[1], B[2], B[3]); \
  }                                    \
} while(0)

#define V1_PERMUTE_B(B) do {           \
  int i;                               \
  for(i = 0; i < MRO_L; ++i) {       \
    V1_DIAG_F(B[0], B[1], B[2], B[3]); \
    V1_G_B(B[0], B[1], B[2], B[3]);    \
    V1_DIAG_B(B[0], B[1], B[2], B[3]); \
    V1_G_B(B[0], B[1], B[2], B[3]);    \
  }                                    \
} while(0)


#define V1_ALPHA_UPDATE_1(L) do {             \
  L[16] = ROT64(L[0], 11) ^ (L[5] << 13);   \
} while(0)

#define V1_ALPHA_UPDATE_2(L) do {     \
  int i;                            \
  for(i = 0; i < 16; ++i) {         \
    L[i] = L[i+1];                  \
  }                                 \
} while(0)

#if 0
#define V1_ALPHA_UPDATE(L) do {                        \
  const uint64_t t = ROT64(L[0], 11) ^ (L[5] << 13); \
  /* int i;  */                                      \
  memmove(&L[0], &L[1], 15 * sizeof(L[0]));          \
  /* This triggers codegen bug on GCC! */            \
  /* for(i = 0; i < 15; ++i) L[i] = L[i+1]; */       \
  /* Also acceptable:                     */         \
  /* STOREU256(&L[ 0], LOADU256(&L[ 1])); */         \
  /* STOREU256(&L[ 4], LOADU256(&L[ 5])); */         \
  /* STOREU256(&L[ 8], LOADU256(&L[ 9])); */         \
  /* STOREU256(&L[12], LOADU256(&L[13])); */         \
  L[15] = t;                                         \
} while(0)
#endif

#define V1_ALPHA_UPDATE(L) do {             \
  L[16] = ROT64(L[0], 11) ^ (L[5] << 13); \
  STOREU256(&L[ 0], LOADU256(&L[ 1]));    \
  STOREU256(&L[ 4], LOADU256(&L[ 5]));    \
  STOREU256(&L[ 8], LOADU256(&L[ 9]));    \
  STOREU256(&L[12], LOADU256(&L[13]));    \
} while(0)

#define V1_BETA_UPDATE(L) do {                                  \
  L[16] = ROT64(L[0], 11) ^ (L[5] << 13);                        \
  STOREU256(&L[ 0], XOR256(LOADU256(&L[ 0]), LOADU256(&L[ 1]))); \
  STOREU256(&L[ 4], XOR256(LOADU256(&L[ 4]), LOADU256(&L[ 5]))); \
  STOREU256(&L[ 8], XOR256(LOADU256(&L[ 8]), LOADU256(&L[ 9]))); \
  STOREU256(&L[12], XOR256(LOADU256(&L[12]), LOADU256(&L[13]))); \
} while(0)

#define V1_BETA2_UPDATE(L) do {                                 \
  L[16] = ROT64(L[0], 11) ^ (L[5] << 13);                        \
  L[17] = ROT64(L[1], 11) ^ (L[6] << 13);                        \
  STOREU256(&L[ 0], XOR256(LOADU256(&L[ 0]), LOADU256(&L[ 2]))); \
  STOREU256(&L[ 4], XOR256(LOADU256(&L[ 4]), LOADU256(&L[ 6]))); \
  STOREU256(&L[ 8], XOR256(LOADU256(&L[ 8]), LOADU256(&L[10]))); \
  STOREU256(&L[12], XOR256(LOADU256(&L[12]), LOADU256(&L[14]))); \
} while(0)

#define V1_GAMMA_UPDATE(L) do {                                                           \
  L[16] = ROT64(L[0], 11) ^ (L[5] << 13);                                                  \
  L[17] = ROT64(L[1], 11) ^ (L[6] << 13);                                                  \
  STOREU256(&L[ 0], XOR256(XOR256(LOADU256(&L[ 0]), LOADU256(&L[ 1])), LOADU256(&L[ 2]))); \
  STOREU256(&L[ 4], XOR256(XOR256(LOADU256(&L[ 4]), LOADU256(&L[ 5])), LOADU256(&L[ 6]))); \
  STOREU256(&L[ 8], XOR256(XOR256(LOADU256(&L[ 8]), LOADU256(&L[ 9])), LOADU256(&L[10]))); \
  STOREU256(&L[12], XOR256(XOR256(LOADU256(&L[12]), LOADU256(&L[13])), LOADU256(&L[14]))); \
} while(0)

#if 0
#define V1_BETA_UPDATE(L) do {                       \
  int i;                                              \
  L[16] = ROT64(L[0], 11) ^ (L[5] << 13);             \
  for(i = 0; i < 16; ++i) {                           \
    L[i] ^= L[i+1];                                   \
  }                                                   \
} while(0)

#define V1_BETA2_UPDATE(L) do {                      \
  int i;                                              \
  L[16] = ROT64(L[0], 11) ^ (L[5] << 13);             \
  L[17] = ROT64(L[1], 11) ^ (L[6] << 13);             \
  for(i = 0; i < 16; ++i) {                           \
    L[i] ^= L[i+2];                                   \
  }                                                   \
} while(0)

#define V1_BETA3_UPDATE(L) do {                      \
  int i;                                              \
  L[16] = ROT64(L[0], 11) ^ (L[5] << 13);             \
  L[17] = ROT64(L[1], 11) ^ (L[6] << 13);             \
  L[18] = ROT64(L[2], 11) ^ (L[7] << 13);             \
  L[19] = ROT64(L[3], 11) ^ (L[8] << 13);             \
  for(i = 0; i < 16; ++i) {                           \
    L[i] ^= L[i+1] ^ L[i+2] ^ L[i+3];                 \
  }                                                   \
} while(0)

#define V1_GAMMA_UPDATE(L) do {                      \
  int i;                                              \
  L[16] = ROT64(L[0], 11) ^ (L[5] << 13);             \
  L[17] = ROT64(L[1], 11) ^ (L[6] << 13);             \
  for(i = 0; i < 16; ++i) {                           \
    L[i] ^= L[i+1] ^ L[i+2];                          \
  }                                                   \
} while(0)
#endif

#define V1_ZERO_BLOCK(B) do {      \
  int i;                           \
  for(i = 0; i < 4; ++i) {         \
  	B[i] = _mm256_setzero_si256(); \
  }                                \
} while(0)

#define V1_LOAD_BLOCK(B, m) do {  \
  int i;                          \
  for(i = 0; i < 4; ++i) {        \
    B[i] = LOADU256(&m[32 * i]);  \
  }                               \
} while(0)

#define V1_STORE_BLOCK(c, B) do { \
  int i;                          \
  for(i = 0; i < 4; ++i) {        \
    STOREU256(&c[32 * i], B[i]);  \
  }                               \
} while(0)

#define V1_COPY_BLOCK_N(C, B, n) do { \
  int i;                              \
  for(i = 0; i < n; ++i) {            \
    C[i] = B[i];                      \
  }                                   \
} while(0)

#define V1_COPY_BLOCK(C, B) V1_COPY_BLOCK_N(C, B, 4)

#define V1_XOR_MASK(B, L) do {               \
  int i;                                    \
  for(i = 0; i < 4; ++i) {                  \
    B[i] = XOR256(B[i], LOADU256(&L[4*i])); \
  }                                         \
} while(0)

#define V1_ACCUMULATE(T, B) do { \
  int i;                         \
  for(i = 0; i < 4; ++i) {       \
    T[i] = XOR256(B[i], T[i]);   \
  }                              \
} while(0)

#define V1_BLOCKCIPHER_F(B, L) do { \
  V1_XOR_MASK(B, L);                \
  V1_PERMUTE_F(B);                  \
  V1_XOR_MASK(B, L);                \
} while(0)

#define V1_BLOCKCIPHER_B(B, L) do { \
  V1_XOR_MASK(B, L);                \
  V1_PERMUTE_B(B);                  \
  V1_XOR_MASK(B, L);                \
} while(0)

#endif
