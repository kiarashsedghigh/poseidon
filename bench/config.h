#ifndef POSEIDON_BENCH_CONFIG_H
#define POSEIDON_BENCH_CONFIG_H

#include "f251.h"

/* Measured 91-round blocks per version. Also overridable with -DITER=... . */
#ifndef ITER
#define ITER 10000
#endif

/* Untimed calls per version before measuring. */
#ifndef WARMUP
#define WARMUP 100
#endif

#ifndef A0
#define A0 1ULL
#endif
#ifndef A1
#define A1 1ULL
#endif
#ifndef A2
#define A2 2ULL
#endif

#if ITER <= 0
#error "ITER must be positive"
#endif
#if WARMUP < 0
#error "WARMUP must not be negative"
#endif

/* Same raw input for every call of both versions; no domain tag is inserted. */
static const felt_t BENCH_INPUT[3] = {
    {3, 0, 0, 0}, {1, 0, 0, 0}, {2, 0, 0, 0}
};

/* Ordinary values. Set the upper limbs here for full 256-bit multipliers. */
static const felt_t BENCH_MULTIPLIERS[3] = {
    {A0, 0, 0, 0}, {A1, 0, 0, 0}, {A2, 0, 0, 0}
};

#endif
