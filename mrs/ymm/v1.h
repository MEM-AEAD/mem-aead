/*
    MRS - MEM AEAD source code package

    :copyright: (c) 2015 by Philipp Jovanovic and Samuel Neves
    :license: Creative Commons CC0 1.0
*/
#ifndef MRS_YMM_V1_H
#define MRS_YMM_V1_H

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
  for(i = 0; i < MRS_L; ++i) {         \
    V1_G_F(B[0], B[1], B[2], B[3]);    \
    V1_DIAG_F(B[0], B[1], B[2], B[3]); \
    V1_G_F(B[0], B[1], B[2], B[3]);    \
    V1_DIAG_B(B[0], B[1], B[2], B[3]); \
  }                                    \
} while(0)

#define V1_PERMUTE_B(B) do {           \
  int i;                               \
  for(i = 0; i < MRS_L; ++i) {         \
    V1_DIAG_F(B[0], B[1], B[2], B[3]); \
    V1_G_B(B[0], B[1], B[2], B[3]);    \
    V1_DIAG_B(B[0], B[1], B[2], B[3]); \
    V1_G_B(B[0], B[1], B[2], B[3]);    \
  }                                    \
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

#define V1_LOAD_BLOCK_N(B, m, n) do {  \
  int i;                            \
  for(i = 0; i < n; ++i) {          \
    B[i] = LOADU256(&m[32 * i]);    \
  }                                 \
} while(0)

#define V1_STORE_BLOCK(c, B) do { \
  int i;                          \
  for(i = 0; i < 4; ++i) {        \
    STOREU256(&c[32 * i], B[i]);  \
  }                               \
} while(0)

#define V1_STORE_BLOCK_N(c, B, n) do { \
  int i;                               \
  for(i = 0; i < n; ++i) {             \
    STOREU256(&c[32 * i], B[i]);       \
  }                                    \
} while(0)

#define V1_COPY_BLOCK_N(C, B, n) do { \
  int i;                              \
  for(i = 0; i < n; ++i) {            \
    C[i] = B[i];                      \
  }                                   \
} while(0)

#define V1_COPY_BLOCK(C, B) V1_COPY_BLOCK_N(C, B, 4)

#define V1_ACCUMULATE(T, B) do { \
  int i;                         \
  for(i = 0; i < 4; ++i) {       \
    T[i] = XOR256(B[i], T[i]);   \
  }                              \
} while(0)

#define V1_ACCUMULATE_N(T, B, n) do { \
  int i;                              \
  for(i = 0; i < n; ++i) {            \
    T[i] = XOR256(B[i], T[i]);        \
  }                                   \
} while(0)

#endif
