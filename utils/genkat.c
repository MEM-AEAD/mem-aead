/*
   STORM reference source code package - reference C implementations

   Written in 2014 by Samuel Neves <sneves@dei.uc.pt>
   Modified in 2015 by Philipp Jovanovic <philipp@jovanovic.io>

   To the extent possible under law, the author(s) have dedicated all copyright
   and related and neighboring rights to this software to the public domain
   worldwide. This software is distributed without any warranty.

   You should have received a copy of the CC0 Public Domain Dedication along with
   this software. If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
*/
#include <string.h>
#include <stdio.h>

void storm_aead_encrypt(
        unsigned char *c, size_t *clen,
        const unsigned char *h, size_t hlen,
        const unsigned char *p, size_t plen,
        const unsigned char *nonce,
        const unsigned char *key);

int storm_aead_decrypt(
        unsigned char *p, size_t *plen,
        const unsigned char *h, size_t hlen,
        const unsigned char *c, size_t clen,
        const unsigned char *nonce,
        const unsigned char *key);

static void genkat(void)
{
	unsigned char w[256];
	unsigned char h[256];
	unsigned char k[32];
	unsigned char n[16];

	unsigned int i, j;

	for(i = 0; i < sizeof w; ++i)
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
	for(i = 0; i < sizeof w; ++i)
	{
		unsigned char m[256];
		unsigned char c[256 + 32];
		size_t mlen;
		size_t clen;
		size_t hlen;

		memset(m, 0, sizeof m);
		memcpy(m, w, i);

		clen = 0;
		mlen = hlen = i;

		storm_aead_encrypt(c, &clen, h, hlen, m, mlen, n, k);

		for(j = 0; j < clen; ++j)
			printf("0x%02X%s", c[j], (j + 1 == clen) ? "" : (7 == j%8) ? ",\n" : ", ");

		printf("%s", (i + 1 == sizeof w) ? "\n" : ",\n\n");
	}
	printf("};\n\n");
	printf("#endif\n\n");
}

int main()
{
	genkat();
	return 0;
}

