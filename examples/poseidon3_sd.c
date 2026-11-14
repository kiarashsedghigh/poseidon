#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include "posedion3_p5_sd.h"

/* Demo uses 64-bit multipliers; the library accepts full 256-bit felt_t values.
 * Usage: poseidon3-sd [A0 A1 A2], defaults: 1 1 2.
 * Rotation and scaling are compile-time options.
 */
int main(int argc, char **argv)
{
    felt_t multipliers[3] = {{1, 0, 0, 0}, {1, 0, 0, 0}, {2, 0, 0, 0}};
    felt_t state[3] = {{3, 0, 0, 0}, {1, 0, 0, 0}, {2, 0, 0, 0}};
    poseidon3_ctx_sd ctx;
    if (argc != 1 && argc != 4) {
        fprintf(stderr, "Usage: %s [A0 A1 A2]\n", argv[0]);
        return 1;
    }
    if (argc == 4) {
        for (int i = 0; i < 3; i++) {
            char *end;
            errno = 0;
            unsigned long long value = strtoull(argv[i+1], &end, 0);
            if (errno || end == argv[i+1] || *end || argv[i+1][0] == '-' ||
                value > UINT64_MAX) {
                fprintf(stderr, "Multipliers must be unsigned 64-bit integers.\n");
                return 1;
            }
            multipliers[i][0] = (uint64_t)value;
        }
    }
    if (poseidon3_init_sd(&ctx, multipliers)) {
        fprintf(stderr, "Multipliers must be nonzero modulo the field prime.\n");
        return 1;
    }
    permutation_3_sd(state, &ctx);
    printf("Compiled SD options: rotation=%d, scaling=%d\n", SD_PERMUTE, SD_SCALE);
#if SD_SCALE
    printf("SD fifth power, first-element S-box; A=(%" PRIu64
           ",%" PRIu64 ",%" PRIu64 ")\n",
           multipliers[0][0], multipliers[1][0], multipliers[2][0]);
#else
    printf("Effective multipliers: (1, 1, 1); arguments ignored.\n");
#endif
    for (int i = 0; i < 3; i++) {
        printf("state[%d] = 0x", i);
        for (int limb = 3; limb >= 0; limb--) printf("%016" PRIx64, state[i][limb]);
        printf("\n");
    }
    return 0;
}
