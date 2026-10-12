# Three-element SD experiment

The new implementation is in `sources/poseidon_sd.c` and `poseidon_sd.h`.
Original functions remain available unchanged. Every new C function has an
`_sd` suffix. The new libraries are `lib_pos_sd.so` and, on x86_64,
`lib_pos_sd_asm.so` (requires ADX at runtime).

## Parameters

As requested, the default SD path combines:

- **C parameters:** p = 2^251 + 17*2^192 + 1, the original 91x3 round constants,
  4 full / 83 partial / 4 full rounds, and base matrix
  `[[3,1,1],[1,-1,1],[1,1,-2]]`.
- **Uploaded Python SD behavior:** fifth power, the first element active in
  partial rounds, per-element multipliers, optional `[2,0,1]` rotation, and
  separate full/partial matrix adjustments.

This is not the uploaded Python's original BLS12-381 parameter set, and is not
the original C cube permutation. In this C field, `gcd(5,p-1)=5`, so fifth power
is not invertible. This experimental construction must not be described as a
valid Poseidon permutation. Tests establish functional agreement, not resistance
to timing, power, electromagnetic, or other side-channel attacks. Configuration
is public; the code does not generate fresh random masks. `SD_CONST` is unused
in the uploaded Python and has no effect here either.

## Configure and run

```c
#include "poseidon_sd.h"

poseidon3_ctx_sd ctx;
felt_t multipliers[3] = {{1,0,0,0}, {1,0,0,0}, {2,0,0,0}};
bool permute = true;  /* false disables rotation */
if (poseidon3_init_sd(&ctx, multipliers, permute) != 0) {
    /* A multiplier is zero modulo p. Do not use an uninitialized ctx. */
    return 1;
}
felt_t state[3] = {{3,0,0,0}, {1,0,0,0}, {2,0,0,0}};
permutation_3_sd(state, &ctx);
/* state is ordinary field representation, still scaled/permuted. */
```

Initialize once per multiplier/flag selection and reuse the read-only context.
Multipliers can be full 256-bit values, stored in four little-endian limbs. They
are reduced modulo p; zero modulo p is rejected without changing the context.
Use a separate context for each configuration; there is no global mutable flag.

To follow the Python hash interface instead, call
`poseidon_hash_3_sd(output, preimage, &ctx)` with two ordinary field elements.
It initializes `[3, preimage[0], preimage[1]]` and returns final `state[1]`.
There is no final inverse scaling/rotation in this interface, matching Python.
Consequently this output can change when multipliers or `permute` change.

For single-round experiments, call `prepare_state_3_sd` **once**, then
`round_3_sd(state, round_idx, FULL_ROUND_SD or PARTIAL_ROUND_SD, &ctx)`.
The state stays in Montgomery form between rounds. Round indices are 0..90.
Full rounds are 0..3 and 87..90. Invalid indices/modes return -1 without mutation.
To inspect an encoded state, convert each element with `f251_from_montgomery`.
`restore_state_3_sd` instead removes scaling/rotation and converts to ordinary
form. Do not pass already prepared state to `permutation_3_sd`.

## Matrix and constant transformations

Let `D = diag(A0,A1,A2)`, `M` be C's base matrix, and `Pi` rotate `[s0,s1,s2]`
to `[s2,s0,s1]` when enabled (otherwise Pi is the identity).

- Prepared input and each constant vector: `Pi * D * vector`.
- Full-round matrix: `Pi * D * M * diag(A0^-5,A1^-5,A2^-5) * Pi^-1`.
- Partial-round matrix: `Pi * D * M * diag(A0^-5,A1^-1,A2^-1) * Pi^-1`.
- Partial S-box position: index 0 without rotation, index 1 with rotation.

These formulas implement the uploaded per-entry matrix modifications for any
nonzero multipliers. Both matrices and all scaled constants are precomputed
during initialization. The fifth-power helper and matrix products use the
existing C or assembly Montgomery backend, with SD-only reduction helpers.

`poseidon3_init_custom_sd(..., 3, 2)` selects cube and the last-element partial
S-box for comparisons with the old C implementation; its matrices adjust
accordingly. With this configuration, `permutation_3_restored_sd` equals the
old `permutation_3` output. With the default `(5,0)` configuration it returns
the unscaled/unrotated fifth-power result, not the original cube result.

## Build, demo, and tests

From the repository root:

```sh
make -C sources clean
make -C sources CC=gcc EXTRA_CFLAGS=-DPERF_MEASUREMENT=100
cc -Wall -Wextra -Isources examples/poseidon3_sd.c -Lsources \
  -Wl,-rpath,"$PWD/sources" -l_pos_sd -o /tmp/poseidon3_demo_sd
/tmp/poseidon3_demo_sd             # defaults: permute=true, A=(1,1,2)
/tmp/poseidon3_demo_sd 0 3 5 7    # permute=false, A=(3,5,7)
python3 tests/test_poseidon3_sd.py
python3 tests/test_poseidon3_sd.py sources/lib_pos_sd_asm.so
```

Use `ISO_C=1` on the build command after cleaning to exercise the ISO C backend.
Tests compare matrices, all round constants, all intermediate rounds, final
states, and the uploaded Python hash interface after substituting C parameters.
They cover both flags, nontrivial and boundary multipliers, boundary/random
inputs, invalid configurations, context reuse, and cube-mode compatibility.
