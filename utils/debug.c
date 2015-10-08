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

void print_bytes(const uint8_t * in, size_t inlen);

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
