#define _POSIX_C_SOURCE 200809L

#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "config.h"
#include "poseidon3_p5_plain.h"
#include "posedion3_p5_sd.h"

/* Consume every output outside the measured interval. */
static volatile uint64_t result_sink;

static double elapsed_ns(const struct timespec *start, const struct timespec *end)
{
    return (double)(end->tv_sec - start->tv_sec) * 1e9
         + (double)(end->tv_nsec - start->tv_nsec);
}

int main(void)
{
    poseidon3_ctx_sd ctx;
    felt_t state[3];
    struct timespec start, end;
    double plain_ns = 0, sd_ns = 0;
    const uint64_t iterations = (uint64_t)ITER;
    const uint64_t warmup = (uint64_t)WARMUP;

    // Precompute SD constants/matrices once, outside all measured intervals.
    if (poseidon3_init_sd(&ctx, BENCH_MULTIPLIERS, SD_PERMUTE)) {
        fprintf(stderr, "SD multipliers must be nonzero modulo the field prime.\n");
        return 1;
    }

    for (uint64_t i = 0; i < warmup; i++) {
        memcpy(state, BENCH_INPUT, sizeof(state));
        permutation_3_pow5_standard(state);
        result_sink += state[0][0] ^ state[1][0] ^ state[2][0];
        memcpy(state, BENCH_INPUT, sizeof(state));
        permutation_3_sd(state, &ctx);
        result_sink += state[0][0] ^ state[1][0] ^ state[2][0];
    }

    // Reset the input before each timing interval; do not chain outputs.
    for (uint64_t i = 0; i < iterations; i++) {
        memcpy(state, BENCH_INPUT, sizeof(state));
        if (clock_gettime(CLOCK_MONOTONIC, &start)) goto clock_error;
        permutation_3_pow5_standard(state);
        if (clock_gettime(CLOCK_MONOTONIC, &end)) goto clock_error;
        plain_ns += elapsed_ns(&start, &end);
        result_sink += state[0][0] ^ state[1][0] ^ state[2][0];
    }
    for (uint64_t i = 0; i < iterations; i++) {
        memcpy(state, BENCH_INPUT, sizeof(state));
        if (clock_gettime(CLOCK_MONOTONIC, &start)) goto clock_error;
        permutation_3_sd(state, &ctx);
        if (clock_gettime(CLOCK_MONOTONIC, &end)) goto clock_error;
        sd_ns += elapsed_ns(&start, &end);
        result_sink += state[0][0] ^ state[1][0] ^ state[2][0];
    }

    printf("Iterations per version: %" PRIu64 " (warmup: %" PRIu64 ")\n",
           iterations, warmup);
    printf("SD rotation: %s; multipliers (hex):\n", SD_PERMUTE ? "on" : "off");
    for (int i = 0; i < 3; i++) {
        printf("  A%d = 0x", i);
        for (int j = 3; j >= 0; j--) printf("%016" PRIx64, BENCH_MULTIPLIERS[i][j]);
        putchar('\n');
    }
    printf("Plain: total %.3f ms, average %.3f us/permutation\n",
           plain_ns / 1e6, plain_ns / (double)iterations / 1e3);
    printf("SD:    total %.3f ms, average %.3f us/permutation\n",
           sd_ns / 1e6, sd_ns / (double)iterations / 1e3);
    if (plain_ns > 0) printf("SD/plain time ratio: %.3fx\n", sd_ns / plain_ns);
    printf("Result checksum: %016" PRIx64 "\n", result_sink);
    return 0;

clock_error:
    perror("clock_gettime(CLOCK_MONOTONIC)");
    return 1;
}
