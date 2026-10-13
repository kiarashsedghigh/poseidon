# Three-element SD permutation

`poseidon3-sd` runs `sources/posedion3_p5_sd.c`. The caller provides three state
values; all three output values are returned in place. There is no hash wrapper,
implicit domain tag, selected output element, or automatic inverse SD transform.

## Parameters and controls

The path uses the C prime `p=2^251+17*2^192+1`, original 91x3 constants,
4 full / 83 partial / 4 full rounds, and base matrix
`[[3,1,1],[1,-1,1],[1,1,-2]]`. Both it and the standard baseline use fifth power
and the first-element partial S-box (index 1 after rotation).

```c
#include "posedion3_p5_sd.h"

poseidon3_ctx_sd ctx;
felt_t multipliers[3] = {{1,0,0,0}, {1,0,0,0}, {2,0,0,0}};
if (poseidon3_init_sd(&ctx, multipliers) != 0) return 1;
felt_t state[3] = {{3,0,0,0}, {1,0,0,0}, {2,0,0,0}};
permutation_3_sd(state, &ctx);
```

Rotation and scaling are selected at compile time: CMake options `SD_PERMUTE`
and `SD_SCALE` default to `ON`. Direct compiler definitions use `0` or `1`.
All SD callers must use the same definitions as the implementation. The runtime
rotation argument and `ctx.permute` have been removed; rebuild all callers.
The exported constants `poseidon3_permute_sd` and `poseidon3_scale_sd` report
the options compiled into a linked library for inspection/testing.

When scaling is enabled, multipliers are ordinary 256-bit integers in four
little-endian limbs. Zero modulo p is rejected without changing the context.
When scaling is disabled, multipliers are ignored (the pointer may be NULL),
effective multipliers are `(1,1,1)`, and the scaling/inversion setup is compiled
out. `SD_PERMUTE=1` selects `[s0,s1,s2] -> [s2,s0,s1]`, with no runtime flag.
Initialize once per multiplier configuration and reuse the read-only context.

The public API contains `poseidon3_init_sd`, `permutation_3_sd`, `round_3_sd`,
and `mix_layer_3_sd`. The runtime structure follows the original C code:

```text
permutation_3_sd
  -> round_3_sd
       -> f251_montgomery_pow5
       -> mix_layer_3_sd
```

The permutation converts/scales/rotates the input once, executes the rounds
directly, then converts all three outputs from Montgomery form. Round and mix
functions take Montgomery states. Round indices are 0..90; modes are
`FULL_ROUND_SD` and `PARTIAL_ROUND_SD`. Invalid index/mode returns -1 without
mutation. The round adds transformed constants to every element, applies the
S-box, and calls the full or partial mix layer.

When `SD_SCALE=0`, `mix_layer_3_sd` applies the base matrix using only
additions/subtractions; rotation selects fixed lane indices at compile time.
Both round modes use this path. With scaling enabled, the function uses the
precomputed full/partial matrices and nine Montgomery multiplications. Matrix
coefficients are canonical after initialization. Each input is normalized once;
each completed row is reduced once, for six canonical reductions per mix.
With both operands below p, each Montgomery product is below 2p, and a sum of
three is below 6p < 2^256. This permits accumulating a row without intermediate
canonical reductions. The input copy preserves in-place operation.

## Matrix transformations

Let `D=diag(A0,A1,A2)`, `M` be the base matrix, and `Pi` the optional rotation.
Initialization precomputes:

- Constants: `Pi * D * original_constants`.
- Full matrix: `Pi * D * M * diag(A0^-5,A1^-5,A2^-5) * Pi^-1`.
- Partial matrix: `Pi * D * M * diag(A0^-5,A1^-1,A2^-1) * Pi^-1`.

For input `x`, output equals `Pi * D * standard(x)`. The C API does not undo
this transform. Tests decode it in Python for comparison. Hash, restoration,
custom-exponent, and extra round-runner wrappers added earlier are removed.
Only the field helpers needed for configurable matrix setup and safe arithmetic
remain private.

This retains the uploaded Python SD transformations but uses the C parameters.
Fifth power is not invertible in this field. Functional tests do not establish
physical side-channel resistance, and the code generates no fresh random masks.
`SD_CONST` is unused in the uploaded Python and has no effect here.

## Run and test

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target poseidon3-sd
./build/poseidon3-sd             # rotation on, A=(1,1,2)
./build/poseidon3-sd 3 5 7       # compiled options, A=(3,5,7)
```

For reference tests, build the shared libraries with the existing Makefile:

```sh
make -C sources clean
make -C sources CC=gcc EXTRA_CFLAGS=-DPERF_MEASUREMENT=100
python3 tests/test_poseidon3_sd.py
# On x86_64 with ADX:
python3 tests/test_poseidon3_sd.py sources/lib_pos_sd_asm.so
```

Use `ISO_C=1` after cleaning for the ISO C backend. Tests compare constants,
matrices, every intermediate round, and all final outputs with the uploaded
Python arithmetic/matrix reference adapted to the C parameters. They read the compiled options from the linked library and cover boundary and
random inputs, multiple multiplier sets, invalid arguments, and context reuse.
Run `python3 tests/test_poseidon3_build_options.py` to build and test all four
combinations on default and portable C backends. See `../bench/README.md` for
commands to benchmark each combination against plain.
