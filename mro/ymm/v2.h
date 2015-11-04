/*
    MRO - MEM AEAD source code package

    :copyright: (c) 2015 by Philipp Jovanovic and Samuel Neves
    :license: Creative Commons CC0 1.0
*/
#ifndef MRO_YMM_V2_H
#define MRO_YMM_V2_H

#include "v0.h"

#define V2_G_F(A, B, C, D) do {                          \
  A = ADD128(A, B); D = XOR128(D, A); D = ROT128(D, 32); \
  C = ADD128(C, D); B = XOR128(B, C); B = ROT128(B, 24); \
  A = ADD128(A, B); D = XOR128(D, A); D = ROT128(D, 16); \
  C = ADD128(C, D); B = XOR128(B, C); B = ROT128(B, 63); \
} while(0)

#define V2_G_B(A, B, C, D) do {                          \
  B = ROT128(B,  1); B = XOR128(B, C); C = SUB128(C, D); \
  D = ROT128(D, 48); D = XOR128(D, A); A = SUB128(A, B); \
  B = ROT128(B, 40); B = XOR128(B, C); C = SUB128(C, D); \
  D = ROT128(D, 32); D = XOR128(D, A); A = SUB128(A, B); \
} while(0)

#define V2_PERMUTE_F(B) do {            \
  int i;                                \
  for(i = 0; i < MRO_L; ++i) {          \
    /* Column step */                   \
    V2_G_F(B[ 0], B[ 4], B[ 8], B[12]); \
    V2_G_F(B[ 1], B[ 5], B[ 9], B[13]); \
    V2_G_F(B[ 2], B[ 6], B[10], B[14]); \
    V2_G_F(B[ 3], B[ 7], B[11], B[15]); \
    /* Diagonal step */                 \
    V2_G_F(B[ 0], B[ 5], B[10], B[15]); \
    V2_G_F(B[ 1], B[ 6], B[11], B[12]); \
    V2_G_F(B[ 2], B[ 7], B[ 8], B[13]); \
    V2_G_F(B[ 3], B[ 4], B[ 9], B[14]); \
  }                                     \
} while(0)

#define V2_PERMUTE_B(B) do {            \
  int i;                                \
  for(i = 0; i < MRO_L; ++i) {          \
    /* Diagonal step */                 \
    V2_G_B(B[ 0], B[ 5], B[10], B[15]); \
    V2_G_B(B[ 1], B[ 6], B[11], B[12]); \
    V2_G_B(B[ 2], B[ 7], B[ 8], B[13]); \
    V2_G_B(B[ 3], B[ 4], B[ 9], B[14]); \
      /* Column step */                 \
    V2_G_B(B[ 0], B[ 4], B[ 8], B[12]); \
    V2_G_B(B[ 1], B[ 5], B[ 9], B[13]); \
    V2_G_B(B[ 2], B[ 6], B[10], B[14]); \
    V2_G_B(B[ 3], B[ 7], B[11], B[15]); \
  }                                     \
} while(0)

#endif
