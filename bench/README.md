# Plain and SD permutation timing

The two implementations are `sources/poseidon3_p5_plain.c` and
`sources/posedion3_p5_sd.c`, with matching header names. Public function names
remain `permutation_3_pow5_standard` and `permutation_3_sd`.

Both versions use the fixed prime defined in `sources/f251_params.h`:
`p = 2^251 + 17*2^192 + 1`. The field backends expose it as `F251_MODULUS`;
`poseidon3_init_sd` copies that value into `ctx.modulus` for inspection.
This context field is ordinary (non-Montgomery), read-only metadata, not a
runtime modulus setting. Rebuild callers after the SD context layout change.

Edit **config.h**:

- `ITER`: measured calls **per version**, default 10000; must be positive.
- `WARMUP`: untimed calls per version, default 100; may be zero.
- `A0`, `A1`, `A2`: SD multipliers, default `(1,1,2)`.
- `BENCH_INPUT`: the three raw input values, default `[3,1,2]`.
- `BENCH_MULTIPLIERS`: edit upper limbs here for full 256-bit multiplier values.

Each array element is stored as four little-endian 64-bit limbs. The macros
can also be overridden by compiler definitions, e.g. `-DITER=100`.

## Compile-time SD options

`SD_PERMUTE` and `SD_SCALE` are CMake options (both default `ON`). They apply
to the SD implementation, example, and benchmark together. Plain is unchanged.
For direct compilation, pass `-DSD_PERMUTE=0|1 -DSD_SCALE=0|1` to every source;
defaults are in `sources/poseidon3_sd_config.h`.

| Configuration | SD_PERMUTE | SD_SCALE |
|---|---|---|
| Neither (SD implementation overhead) | OFF | OFF |
| Rotation only | ON | OFF |
| Scaling only | OFF | ON |
| Both | ON | ON |

When scaling is disabled, the configured multipliers are ignored and effective
multipliers are `(1,1,1)`. Input multiplication, constant scaling, and matrix
scaling/inversion setup are compiled out. All configurations use the same
general matrix loop as plain; initialization chooses the matrix entries. When rotation is disabled, indices use the identity
order; there is no runtime rotation flag or conditional. Multiplier values remain configurable
through `A0`, `A1`, `A2` when scaling is enabled.

Build and run all four comparisons in separate directories:

```sh
for rotate in OFF ON; do
    for scale in OFF ON; do
        build_dir="build/sd-${rotate}-${scale}"
        cmake -S . -B "$build_dir" -DCMAKE_BUILD_TYPE=Release \
            -DSD_PERMUTE="$rotate" -DSD_SCALE="$scale"
        cmake --build "$build_dir" --target perf
        "$build_dir/perf"
    done
done
```

Each run reports the compiled options and compares its SD round loop against
plain. Both versions use identical general matrix loops: nine Montgomery
multiplications, six additions, and six canonical reductions. The former
optimized standard mix body is retained as a comment in the source.

## Build and run

From the repository root on a POSIX system (e.g. Linux):

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target perf
./build/perf
```

In CLion, select a Release CMake profile, reload CMake, and select target `perf`.
The existing `poseidon3-target` and `poseidon3-sd` examples remain available.

The Makefile also builds the benchmark:

```sh
make -C sources perf CC=gcc
./bench/perf
```

Rebuild after editing `config.h`. Clean before switching the Makefile's compiler
flags, SD options, or backend. For example, after cleaning,
`make -C sources perf CC=gcc SD_PERMUTE=0 SD_SCALE=1` builds scaling only. For ISO C use `-DPOSEIDON_ISO_C=ON` with CMake, or `ISO_C=1`
with Make. Both benchmark implementations always use the same selected backend.

## What is measured

`perf.c` runs each 91-round loop `ITER` times and reports summed elapsed
time in milliseconds, arithmetic mean in microseconds per 91-round block, and
the SD/plain time ratio. Benchmark-only helpers `rounds_3_standard` and
`rounds_3_sd` call the existing round functions with the same 4/83/4 schedule.
The public permutation functions are unchanged.

Only round constant addition, fifth-power S-boxes, and mixing are measured.
Conversion to/from Montgomery form, SD input scaling and optional rotation,
and both versions’ matrix/constant precomputation are excluded. Both initial states are
prepared once before warmup. Each invocation resets from its prepared state
before starting the timer, so outputs are not fed back into later invocations.
Outputs stay in Montgomery form.

Timing uses `clock_gettime(CLOCK_MONOTONIC)` around each round-loop call.
Warmup, input resets, output consumption, and printing are excluded. The timer
boundaries have a small overhead and wall time can include scheduling delays;
no overhead subtraction is applied.

A volatile checksum consumes every output outside the timed intervals. Equivalent
noncanonical Montgomery representations can produce different checksums; use
the correctness tests to compare field values. Results
are timing observations for the current compiler/backend and machine, not a
side-channel-resistance measurement. The benchmark itself does not compare raw
outputs, which generally differ because SD retains its scaling/rotation.
Correctness comparisons remain in the separate Python test suites.

Validate every option combination on default and portable C backends with:

```sh
python3 tests/test_poseidon3_build_options.py
```
