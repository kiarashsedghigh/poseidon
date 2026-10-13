#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include "poseidon_pow5_standard.h"
#include "poseidon_sd.h"

/* Python-style preimage [1,2], initialized as [3,1,2]. */
int main(void)
{
    const felt_t input[3] = {{3,0,0,0}, {1,0,0,0}, {2,0,0,0}};
    felt_t baseline[3], candidate[3], hash;
    const felt_t identity[3] = {{1,0,0,0}, {1,0,0,0}, {1,0,0,0}};
    const felt_t multipliers[3] = {{3,0,0,0}, {5,0,0,0}, {7,0,0,0}};
    poseidon3_ctx_sd ctx;

    memcpy(baseline, input, sizeof(baseline));
    permutation_3_pow5_standard(baseline);
    printf("Non-SD fifth-power baseline (C parameters):\n");
    for (int i = 0; i < 3; i++) {
        printf("state[%d] = 0x", i);
        for (int limb = 3; limb >= 0; limb--)
            printf("%016" PRIx64, baseline[i][limb]);
        printf("\n");
    }
    poseidon_hash_3_pow5_standard(hash, input+1);
    if (memcmp(hash, baseline[1], sizeof(hash))) return 1;
    printf("Python-style hash is state[1].\n");

    if (poseidon3_init_sd(&ctx, identity, false)) return 1;
    memcpy(candidate, input, sizeof(candidate));
    permutation_3_sd(candidate, &ctx);
    if (memcmp(candidate, baseline, sizeof(candidate))) {
        fprintf(stderr, "FAIL: SD identity configuration differs.\n");
        return 1;
    }
    if (poseidon3_init_sd(&ctx, multipliers, true)) return 1;
    memcpy(candidate, input, sizeof(candidate));
    permutation_3_restored_sd(candidate, &ctx);
    if (memcmp(candidate, baseline, sizeof(candidate))) {
        fprintf(stderr, "FAIL: restored SD output differs.\n");
        return 1;
    }
    printf("PASS: baseline matches identity SD and restored SD with A=(3,5,7).\n");
    return 0;
}
