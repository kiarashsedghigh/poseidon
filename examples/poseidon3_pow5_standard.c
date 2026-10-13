#include <inttypes.h>
#include <stdio.h>
#include "poseidon_pow5_standard.h"

/* Python-style preimage [1,2], initialized as [3,1,2]. */
int main(void)
{
    felt_t state[3] = {{3,0,0,0}, {1,0,0,0}, {2,0,0,0}};

    permutation_3_pow5_standard(state);
    printf("Non-SD fifth-power baseline (C parameters):\n");
    for (int i = 0; i < 3; i++) {
        printf("state[%d] = 0x", i);
        for (int limb = 3; limb >= 0; limb--)
            printf("%016" PRIx64, state[i][limb]);
        printf("\n");
    }
    printf("Python-style hash is state[1].\n");
    return 0;
}
