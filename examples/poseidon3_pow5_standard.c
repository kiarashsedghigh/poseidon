#include <inttypes.h>
#include <stdio.h>
#include "poseidon3_p5_plain.h"

/* Input state supplied directly: all three elements are user inputs. */
int main(void)
{
    felt_t state[3] = {{3,0,0,0}, {1,0,0,0}, {2,0,0,0}};

    poseidon3_ctx_standard ctx;
    poseidon3_init_standard(&ctx);
    permutation_3_pow5_standard(state, &ctx);
    printf("Non-SD fifth-power baseline (C parameters):\n");
    for (int i = 0; i < 3; i++) {
        printf("state[%d] = 0x", i);
        for (int limb = 3; limb >= 0; limb--)
            printf("%016" PRIx64, state[i][limb]);
        printf("\n");
    }
    return 0;
}
