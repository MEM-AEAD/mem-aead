/*
    MEM AEAD source code package

    :copyright: (c) 2015 by Philipp Jovanovic and Samuel Neves
    :license: Creative Commons CC0 1.0
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "kat.h"

void crypto_aead_encrypt(
        unsigned char *c, size_t *clen,
        const unsigned char *h, size_t hlen,
        const unsigned char *p, size_t plen,
        const unsigned char *nonce,
        const unsigned char *key);

int crypto_aead_decrypt(
        unsigned char *p, size_t *plen,
        const unsigned char *h, size_t hlen,
        const unsigned char *c, size_t clen,
        const unsigned char *nonce,
        const unsigned char *key);

int check(const unsigned char *kat)
{
#define MAX_SIZE 768
    unsigned char w[MAX_SIZE];
    unsigned char h[MAX_SIZE];
    unsigned char k[32];
    unsigned char n[16];

    unsigned i;
    int place = 0;

    for(i = 0; i < MAX_SIZE; ++i)
        w[i] = 255 & (i*197 + 123);

    for(i = 0; i < sizeof h; ++i)
        h[i] = 255 & (i*193 + 123);

    for(i = 0; i < sizeof k; ++i)
        k[i] = 255 & (i*191 + 123);

    for(i = 0; i < sizeof n; ++i)
        n[i] = 255 & (i*181 + 123);

    for(i = 0; i < MAX_SIZE; ++i)
    {
        unsigned char m[MAX_SIZE];
        unsigned char c[MAX_SIZE + 32];
        size_t mlen;
        size_t clen;
        size_t hlen;

        memset(m, 0, sizeof m);
        memcpy(m, w, i);

        clen = 0;
        mlen = hlen = i;

        crypto_aead_encrypt(c, &clen, h, hlen, m, mlen, n, k);
        if( 0 != memcmp(kat, c, clen) ) {place = 1; goto fail;}

        memset(m, 0, sizeof m);
        mlen = 0;

        if( 0 != crypto_aead_decrypt(m, &mlen, h, hlen, c, clen, n, k) )
            {place = 2; goto fail;}

        if( 0 != memcmp(m, w, mlen) ) {place = 3; goto fail;}

        kat += clen;
    }
    printf("ok\n");
    return 0;
fail:
    printf("fail at %u:%d\n", i, place);
    return -1;
#undef MAX_SIZE
}

int main()
{
    return check(kat);
}

