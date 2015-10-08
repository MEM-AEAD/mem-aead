#ifndef MRO_H
#define MRO_H

#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef uint64_t mro_word_t;

typedef struct state__
{
    mro_word_t S[16];
} mro_state_t[1];

typedef enum tag__
{
    ABS_AD  = 0x00,
    ABS_MSG = 0x01
} tag_t;

/* high-level operations */
void mro_aead_encrypt(
        unsigned char *c, size_t *clen,
        const unsigned char *h, size_t hlen,
        const unsigned char *m, size_t mlen,
        const unsigned char *nonce,
        const unsigned char *key);

int mro_aead_decrypt(
        unsigned char *m, size_t *mlen,
        const unsigned char *h, size_t hlen,
        const unsigned char *c, size_t clen,
        const unsigned char *nonce,
        const unsigned char *key);

/*foo*/

#endif
