#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "poseidon_sd.h"

/* Demo uses 64-bit multipliers; the library accepts full 256-bit felt_t values.
 * Usage: demo_sd [0|1 [A0 A1 A2]], defaults: 1 1 1 2.
 */
int main(int argc, char **argv)
{
    bool permute = true;
    felt_t multipliers[3] = {{1, 0, 0, 0}, {1, 0, 0, 0}, {2, 0, 0, 0}};
    felt_t state[3] = {{3, 0, 0, 0}, {1, 0, 0, 0}, {2, 0, 0, 0}};
    poseidon3_ctx_sd ctx;
    if (argc != 1 && argc != 2 && argc != 5) {
        fprintf(stderr, "Usage: %s [0|1 [A0 A1 A2]]\n", argv[0]);
        return 1;
    }
    if (argc >= 2) {
        if (strcmp(argv[1], "0") && strcmp(argv[1], "1")) {
            fprintf(stderr, "Permutation flag must be 0 or 1.\n");
            return 1;
        }
        permute = !strcmp(argv[1], "1");
    }
    if (argc == 5) {
        for (int i = 0; i < 3; i++) {
            char *end;
            errno = 0;
            unsigned long long value = strtoull(argv[i+2], &end, 0);
            if (errno || end == argv[i+2] || *end || argv[i+2][0] == '-' ||
                value > UINT64_MAX) {
                fprintf(stderr, "Multipliers must be unsigned 64-bit integers.\n");
                return 1;
            }
            multipliers[i][0] = (uint64_t)value;
        }
    }
    if (poseidon3_init_sd(&ctx, multipliers, permute)) {
        fprintf(stderr, "Multipliers must be nonzero modulo the field prime.\n");
        return 1;
    }
    permutation_3_sd(state, &ctx);
    printf("SD fifth power, first-element S-box; permute=%d, A=(%" PRIu64
           ",%" PRIu64 ",%" PRIu64 ")\n", (int)permute,
           multipliers[0][0], multipliers[1][0], multipliers[2][0]);
    for (int i = 0; i < 3; i++) {
        printf("state[%d] = 0x", i);
        for (int limb = 3; limb >= 0; limb--) printf("%016" PRIx64, state[i][limb]);
        printf("\n");
    }
    return 0;
}
