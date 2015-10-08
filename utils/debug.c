#include <stddef.h>
#include <stdint.h>
#include <inttypes.h>
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

/*
#define FMT "016" PRIX64

typedef uint64_t word_t;

typedef struct state__
{
    word_t S[16];
} state_t[1];

static void print_state(state_t state)
{
    static const char fmt[] = "%" FMT " "
                              "%" FMT " "
                              "%" FMT " "
                              "%" FMT "\n";
    const word_t * S = state->S;
    printf(fmt, S[ 0],S[ 1],S[ 2],S[ 3]);
    printf(fmt, S[ 4],S[ 5],S[ 6],S[ 7]);
    printf(fmt, S[ 8],S[ 9],S[10],S[11]);
    printf(fmt, S[12],S[13],S[14],S[15]);
    printf("\n");
}*/

static void print_bytes(const unsigned char *in, size_t inlen)
{
    size_t i;
    for (i = 0; i < inlen; ++i) {
        printf("%02X ", in[i]);
        if (i%16 == 15) {
            printf("\n");
        }
    }
    printf("\n");
}


#define HSIZE 1024
#define MSIZE 1024

int main() {
    unsigned char k[32] = {0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,0x08,0x09,0x0A,0x0B,0x0C,0x0D,0x0E,0x0F,0xFF,0xEE,0xDD,0xCC,0xBB,0xAA,0x99,0x88,0x77,0x66,0x55,0x44,0x33,0x22,0x11,0x00};
    unsigned char n[16] = {0xF0,0xE0,0xD0,0xC0,0xB0,0xA0,0x90,0x80,0x70,0x60,0x50,0x40,0x30,0x20,0x10,0x00};
    unsigned char h[HSIZE] = {0};
    unsigned char m[MSIZE] = {0};
    unsigned char c[MSIZE + 32] = {0};
    size_t hlen = 129;
    size_t mlen = 256;
    size_t clen = 0;
    size_t i = 0;
    int result = -1;

    for (i = 0; i < hlen; ++i) { h[i] = i & 255; }
    for (i = 0; i < mlen; ++i) { m[i] = i & 255; }

    printf("========== SETUP ==========\n");
    printf("KEY:\n");
    print_bytes(k, sizeof k);
    printf("NONCE:\n");
    print_bytes(n, sizeof n);
    printf("HEADER:\n");
    print_bytes(h, hlen);
    printf("PAYLOAD:\n");
    print_bytes(m, mlen);
    printf("\n");

    printf("========== ENCRYPTION ==========\n");
    crypto_aead_encrypt(c, &clen, h, hlen, m, mlen, n, k);
    printf("ENCRYPTED PAYLOAD + TAG:\n");
    print_bytes(c, clen);
    printf("\n");

    printf("========== DECRYPTION ==========\n");
    result = crypto_aead_decrypt(m, &mlen, h, hlen, c, clen, n, k);
    printf("DECRYPTED PAYLOAD:\n");
    print_bytes(m, mlen);

    printf("verify: %d\n", result);

    return 0;
}
