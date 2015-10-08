#ifndef OPP_XMM_H
#define OPP_XMM_H

#include <stddef.h>
#include <stdint.h>

/* high-level operations */
void opp_aead_encrypt(
        unsigned char *c, size_t *clen,
        const unsigned char *h, size_t hlen,
        const unsigned char *m, size_t mlen,
        const unsigned char *nonce,
        const unsigned char *key);

int opp_aead_decrypt(
        unsigned char *m, size_t *mlen,
        const unsigned char *h, size_t hlen,
        const unsigned char *c, size_t clen,
        const unsigned char *nonce,
        const unsigned char *key);

#endif
