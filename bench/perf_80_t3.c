#define _POSIX_C_SOURCE 200809L
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "config.h"
#include "poseidon_profiles.h"

static volatile uint64_t result_sink;
static double elapsed_ns(const struct timespec *start, const struct timespec *end)
{
    return (double)(end->tv_sec-start->tv_sec)*1e9 + (double)(end->tv_nsec-start->tv_nsec);
}

int main(void)
{
    const poseidon_profile *profile = NULL;
    for (int i = 0; i < 7; i++)
        if (!strcmp(poseidon_profiles[i].name, "poseidon-80-t3")) profile = &poseidon_profiles[i];
    if (!profile) { fprintf(stderr, "Profile poseidon-80-t3 not found.\n"); return 1; }
    poseidon_profile_ctx plain, sd;
    if (poseidon_profile_init(&plain, profile, 0, NULL) ||
        poseidon_profile_init(&sd, profile, 1, BENCH_MULTIPLIERS)) {
        fprintf(stderr, "Profile initialization failed; multipliers must be nonzero modulo p.\n");
        return 1;
    }
    felt_t plain_input[3], scaled[3], sd_input[3], state[3];
    for (int i = 0; i < 3; i++) {
        f251_to_montgomery(plain_input[i], BENCH_INPUT[i]);
#if SD_SCALE
        f251_montgomery_mult(scaled[i], plain_input[i], sd.multipliers[i]);
#else
        memcpy(scaled[i], plain_input[i], sizeof(felt_t));
#endif
    }
#if SD_PERMUTE
    const int order[3] = {2, 0, 1};
#else
    const int order[3] = {0, 1, 2};
#endif
    for (int i = 0; i < 3; i++) memcpy(sd_input[i], scaled[order[i]], sizeof(felt_t));
    const uint64_t iterations = ITER, warmup = WARMUP;
    for (uint64_t i = 0; i < warmup; i++) {
        memcpy(state, plain_input, sizeof(state));
        poseidon_profile_rounds(&plain, state);
        result_sink += state[0][0] ^ state[1][0] ^ state[2][0];
        memcpy(state, sd_input, sizeof(state));
        poseidon_profile_rounds(&sd, state);
        result_sink += state[0][0] ^ state[1][0] ^ state[2][0];
    }
    double plain_ns = 0, sd_ns = 0;
    struct timespec start, end;
    for (uint64_t i = 0; i < iterations; i++) {
        memcpy(state, plain_input, sizeof(state));
        if (clock_gettime(CLOCK_MONOTONIC, &start)) goto clock_error;
        poseidon_profile_rounds(&plain, state);
        if (clock_gettime(CLOCK_MONOTONIC, &end)) goto clock_error;
        plain_ns += elapsed_ns(&start, &end);
        result_sink += state[0][0] ^ state[1][0] ^ state[2][0];
    }
    for (uint64_t i = 0; i < iterations; i++) {
        memcpy(state, sd_input, sizeof(state));
        if (clock_gettime(CLOCK_MONOTONIC, &start)) goto clock_error;
        poseidon_profile_rounds(&sd, state);
        if (clock_gettime(CLOCK_MONOTONIC, &end)) goto clock_error;
        sd_ns += elapsed_ns(&start, &end);
        result_sink += state[0][0] ^ state[1][0] ^ state[2][0];
    }
    printf("Profile: %s, t=%d, full=%d, partial=%d\n",profile->name,profile->t,profile->rf,profile->rp);
    printf("Iterations per version: %" PRIu64 " (warmup: %" PRIu64 ")\n",iterations,warmup);
    printf("Compiled SD options: rotation=%d, scaling=%d\n",SD_PERMUTE,SD_SCALE);
    printf("Input (hex), and effective SD multipliers (ordinary form):\n");
    for (int i = 0; i < 3; i++) {
        printf("  state[%d] = 0x", i);
        for (int j = 3; j >= 0; j--) printf("%016" PRIx64, BENCH_INPUT[i][j]);
        printf("; A%d = 0x", i);
#if SD_SCALE
        for (int j = 3; j >= 0; j--) printf("%016" PRIx64, BENCH_MULTIPLIERS[i][j]);
#else
        printf("%064x", 1);
#endif
        putchar('\n');
    }
    printf("Plain: total %.3f ms, average %.3f us/41-round block\n",plain_ns/1e6,plain_ns/iterations/1e3);
    printf("SD:    total %.3f ms, average %.3f us/41-round block\n",sd_ns/1e6,sd_ns/iterations/1e3);
    if (plain_ns > 0) {
        printf("SD/plain time ratio: %.3fx\n",sd_ns/plain_ns);
        printf("Overhead: %+.2f%%\n",100*(sd_ns/plain_ns-1));
    }
    printf("Result checksum: %016" PRIx64 "\n",result_sink);
    return 0;
clock_error:
    perror("clock_gettime(CLOCK_MONOTONIC)");
    return 1;
}
