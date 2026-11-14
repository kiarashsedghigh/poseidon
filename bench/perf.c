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

/* Benchmark-only helpers: state is already prepared in Montgomery form. */
static void rounds_3_standard(felt_t state[3], const poseidon3_ctx_standard *ctx)
{
    for (int r = 0; r < POSEIDON3_ROUNDS_STANDARD; r++)
        round_3_pow5_standard(state, r,
            r < 4 || r >= 87 ? FULL_ROUND_STANDARD : PARTIAL_ROUND_STANDARD, ctx);
}

static void rounds_3_sd(felt_t state[3], const poseidon3_ctx_sd *ctx)
{
    for (int r = 0; r < POSEIDON3_ROUNDS_SD; r++)
        round_3_sd(state, r,
            r < 4 || r >= 87 ? FULL_ROUND_SD : PARTIAL_ROUND_SD, ctx);
}

int main(void)
{
    poseidon3_ctx_sd ctx;
    poseidon3_ctx_standard plain_ctx;
    felt_t state[3], plain_input[3], sd_input[3], scaled[3];
    struct timespec start, end;
    double plain_ns = 0, sd_ns = 0;
    const uint64_t iterations = (uint64_t)ITER;
    const uint64_t warmup = (uint64_t)WARMUP;

    // Prepare both contexts outside all measured intervals.
    poseidon3_init_standard(&plain_ctx);
    if (poseidon3_init_sd(&ctx, BENCH_MULTIPLIERS)) {
        fprintf(stderr, "SD multipliers must be nonzero modulo the field prime.\n");
        return 1;
    }

    // Prepare each initial state once; only the round loops are timed.
    for (int i = 0; i < 3; i++) {
        f251_to_montgomery(plain_input[i], BENCH_INPUT[i]);
#if SD_SCALE
        f251_montgomery_mult(scaled[i], plain_input[i], ctx.multipliers[i]);
#else
        memcpy(scaled[i], plain_input[i], sizeof(felt_t));
#endif
    }
#if SD_PERMUTE
    const int order[3] = {2, 0, 1};
#else
    const int order[3] = {0, 1, 2};
#endif
    for (int i = 0; i < 3; i++)
        memcpy(sd_input[i], scaled[order[i]], sizeof(felt_t));

    for (uint64_t i = 0; i < warmup; i++) {
        memcpy(state, plain_input, sizeof(state));
        rounds_3_standard(state, &plain_ctx);
        result_sink += state[0][0] ^ state[1][0] ^ state[2][0];
        memcpy(state, sd_input, sizeof(state));
        rounds_3_sd(state, &ctx);
        result_sink += state[0][0] ^ state[1][0] ^ state[2][0];
    }

    // Reset the input before each timing interval; do not chain outputs.
    for (uint64_t i = 0; i < iterations; i++) {
        memcpy(state, plain_input, sizeof(state));
        if (clock_gettime(CLOCK_MONOTONIC, &start)) goto clock_error;
        rounds_3_standard(state, &plain_ctx);
        if (clock_gettime(CLOCK_MONOTONIC, &end)) goto clock_error;
        plain_ns += elapsed_ns(&start, &end);
        result_sink += state[0][0] ^ state[1][0] ^ state[2][0];
    }
    for (uint64_t i = 0; i < iterations; i++) {
        memcpy(state, sd_input, sizeof(state));
        if (clock_gettime(CLOCK_MONOTONIC, &start)) goto clock_error;
        rounds_3_sd(state, &ctx);
        if (clock_gettime(CLOCK_MONOTONIC, &end)) goto clock_error;
        sd_ns += elapsed_ns(&start, &end);
        result_sink += state[0][0] ^ state[1][0] ^ state[2][0];
    }

    printf("Iterations per version: %" PRIu64 " (warmup: %" PRIu64 ")\n",
           iterations, warmup);
    printf("Compiled SD options: rotation=%d, scaling=%d\n", SD_PERMUTE, SD_SCALE);
#if SD_SCALE
    printf("Multipliers (hex):\n");
    for (int i = 0; i < 3; i++) {
        printf("  A%d = 0x", i);
        for (int j = 3; j >= 0; j--) printf("%016" PRIx64, BENCH_MULTIPLIERS[i][j]);
        putchar('\n');
    }
#else
    printf("Effective multipliers: (1, 1, 1); configured multipliers ignored.\n");
#endif
    printf("Plain: total %.3f ms, average %.3f us/91-round block\n",
           plain_ns / 1e6, plain_ns / (double)iterations / 1e3);
    printf("SD:    total %.3f ms, average %.3f us/91-round block\n",
           sd_ns / 1e6, sd_ns / (double)iterations / 1e3);
    if (plain_ns > 0) printf("SD/plain time ratio: %.3fx\n", sd_ns / plain_ns);
    printf("Result checksum: %016" PRIx64 "\n", result_sink);
    return 0;

clock_error:
    perror("clock_gettime(CLOCK_MONOTONIC)");
    return 1;
}
