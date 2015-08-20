#ifndef STORM_OPP_YMM_V1_H
#define STORM_OPP_YMM_V1_H

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
  for(i = 0; i < STORM_R; ++i) {       \
    V1_G_F(B[0], B[1], B[2], B[3]);    \
    V1_DIAG_F(B[0], B[1], B[2], B[3]); \
    V1_G_F(B[0], B[1], B[2], B[3]);    \
    V1_DIAG_B(B[0], B[1], B[2], B[3]); \
  }                                    \
} while(0)

#define V1_PERMUTE_B(B) do {           \
  int i;                               \
  for(i = 0; i < STORM_R; ++i) {       \
    V1_DIAG_F(B[0], B[1], B[2], B[3]); \
    V1_G_B(B[0], B[1], B[2], B[3]);    \
    V1_DIAG_B(B[0], B[1], B[2], B[3]); \
    V1_G_B(B[0], B[1], B[2], B[3]);    \
  }                                    \
} while(0)


#define V1_MASK_UPDATE_1(L) do {            \
  L[16] = ROT64(L[0], 11) ^ (L[5] << 13);   \
} while(0)

#define V1_MASK_UPDATE_2(L) do {    \
  int i;                            \
  for(i = 0; i < 16; ++i) {         \
    L[i] = L[i+1];                  \
  }                                 \
} while(0)

#define V1_MASK_UPDATE(L) do {                       \
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

#define V1_MASK_ROT_768(L) do {        \
	const __m256i t = LOADU256(&L[ 0]);  \
	STOREU256(&L[ 0], LOADU256(&L[ 4])); \
	STOREU256(&L[ 4], LOADU256(&L[ 8])); \
	STOREU256(&L[ 8], LOADU256(&L[12])); \
	STOREU256(&L[12], t);                \
} while(0)

#define V1_MASK_ROT_256(L) do {        \
	const __m256i t = LOADU256(&L[12]);  \
  STOREU256(&L[12], LOADU256(&L[ 8])); \
  STOREU256(&L[ 8], LOADU256(&L[ 4])); \
  STOREU256(&L[ 4], LOADU256(&L[ 0])); \
  STOREU256(&L[ 0], t);                \
} while(0)

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

#define V1_XOR_MASK(B, L) do {              \
  int i;                                    \
  for(i = 0; i < 4; ++i) {                  \
    B[i] = XOR256(B[i], LOADU256(&L[4*i])); \
  }                                         \
} while(0)

#define V1_XOR_ROTATED_MASK(B, L, R) do {                         \
  int i;                                                          \
  for(i = 0; i < 4; ++i) {                                        \
    B[i] = XOR256(B[i], LOADU256(&L[(4*i + (16 - 4*R/256))%16])); \
  }                                                               \
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

#define V1_BLOCKCIPHER_ROTATED_F(B, L, R) do { \
  V1_XOR_ROTATED_MASK(B, L, R);                \
  V1_PERMUTE_F(B);                             \
  V1_XOR_ROTATED_MASK(B, L, R);                \
} while(0)

#define V1_BLOCKCIPHER_ROTATED_B(B, L, R) do { \
  V1_XOR_ROTATED_MASK(B, L, R);                \
  V1_PERMUTE_B(B);                             \
  V1_XOR_ROTATED_MASK(B, L, R);                \
} while(0)

#endif
