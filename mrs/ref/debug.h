#ifndef MRS_DEBUG_H
#define MRS_DEBUG_H

#include <inttypes.h>
#include <stdio.h>

#define FMT "016" PRIX64

static void print_state(mrs_state_t state)
{
    static const char fmt[] = "%" FMT " "
                              "%" FMT " "
                              "%" FMT " "
                              "%" FMT "\n";
    const mrs_word_t * S = state->S;
    printf(fmt, S[ 0],S[ 1],S[ 2],S[ 3]);
    printf(fmt, S[ 4],S[ 5],S[ 6],S[ 7]);
    printf(fmt, S[ 8],S[ 9],S[10],S[11]);
    printf(fmt, S[12],S[13],S[14],S[15]);
    printf("\n");
}

static void print_bytes(const uint8_t * in, size_t inlen)
{
    size_t i;
    for (i = 0; i < inlen; ++i)
    {
        printf("%02X ", in[i]);
        if (i%16 == 15)
        {
            printf("\n");
        }
    }
    printf("\n");
}

#endif
