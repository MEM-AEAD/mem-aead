/*
   STORM reference source code package - reference C implementations

   Written in 2015 by Philipp Jovanovic <philipp@jovanovic.io>

   To the extent possible under law, the author(s) have dedicated all copyright
   and related and neighboring rights to this software to the public domain
   worldwide. This software is distributed without any warranty.

   You should have received a copy of the CC0 Public Domain Dedication along with
   this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/
#ifndef STORM_H
#define STORM_H

#include <stddef.h>
#include <stdint.h>
#include "storm_config.h"

#if STORM_W == 64
	typedef uint64_t storm_word_t;
#else
	#error "Invalid word size!"
#endif

typedef struct state__
{
    storm_word_t S[16];
} storm_state_t[1];


/* high-level operations */
void storm_aead_encrypt(
        unsigned char *c, size_t *clen,
        const unsigned char *h, size_t hlen,
        const unsigned char *m, size_t mlen,
        const unsigned char *t, size_t tlen,
        const unsigned char *nonce,
        const unsigned char *key);

int storm_aead_decrypt(
        unsigned char *m, size_t *mlen,
        const unsigned char *h, size_t hlen,
        const unsigned char *c, size_t clen,
        const unsigned char *t, size_t tlen,
        const unsigned char *nonce,
        const unsigned char *key);


#endif
