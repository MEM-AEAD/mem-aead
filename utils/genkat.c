/*
    MEM AEAD source code package

    :copyright: (c) 2015 by Philipp Jovanovic and Samuel Neves
    :license: Creative Commons CC0 1.0
*/
#include <string.h>
#include <stdio.h>

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

static void genkat(void)
{
#define MAX_SIZE 768
	unsigned char w[MAX_SIZE];
	unsigned char h[MAX_SIZE];
	unsigned char k[32];
	unsigned char n[16];

	unsigned int i, j;

	for(i = 0; i < MAX_SIZE; ++i)
		w[i] = 255 & (i*197 + 123);

	for(i = 0; i < sizeof h; ++i)
		h[i] = 255 & (i*193 + 123);

	for(i = 0; i < sizeof k; ++i)
		k[i] = 255 & (i*191 + 123);

	for(i = 0; i < sizeof n; ++i)
		n[i] = 255 & (i*181 + 123);

	printf("#ifndef STORM_KAT_H\n");
	printf("#define STORM_KAT_H\n");
	printf("static const unsigned char kat[] = \n{\n");
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

		for(j = 0; j < clen; ++j)
			printf("0x%02X%s", c[j], (j + 1 == clen) ? "" : (7 == j%8) ? ",\n" : ", ");

		printf("%s", (i + 1 == sizeof w) ? "\n" : ",\n\n");
	}
	printf("};\n\n");
	printf("#endif\n\n");
#undef MAX_SIZE
}

int main()
{
	genkat();
	return 0;
}

