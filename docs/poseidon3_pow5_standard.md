# Non-SD fifth-power baseline

`sources/poseidon_pow5_standard.c` is an independent implementation for direct
comparison with the default SD implementation. It does not call SD functions,
scale inputs/constants, rotate the state, or adjust the MDS matrix. All new
library function names end in `_standard`. Existing functions are preserved.

| Implementation | Exponent | Partial S-box position | SD transformations |
|---|---|---|---|
| Existing `round_3_standard` | 3 | Last element (2) | None |
| New `round_3_pow5_standard` | 5 | First element (0) | None |
| Default `round_3_sd` | 5 | 0 without rotation, 1 with rotation | Scaling, rotation, adjusted matrices |

All three use the C field `p=2^251+17*2^192+1`, 4 full / 83 partial / 4 full
rounds and the original 91x3 round constants. The new baseline uses C's base
MDS `[[3,1,1],[1,-1,1],[1,1,-2]]`. It matches the Python round structure and
fifth-power convention with these C parameters, not the uploaded BLS12-381/
Cauchy-matrix/55-partial-round parameter set. As with SD, fifth power is not
invertible over this C field; this remains an experimental construction.

## API

Include `poseidon_pow5_standard.h` and link `lib_pos_pow5_standard.so`.

- `round_3_pow5_standard(state, round_idx, round_mode)`: one round on a Montgomery
  state. Adds all three constants, applies fifth power, then the base matrix.
  Index is 0..90; mode is `FULL_ROUND_STANDARD` or `PARTIAL_ROUND_STANDARD`.
  Invalid index/mode returns -1 without mutation; success returns 0.
- `permutation_3_pow5_standard(state)`: ordinary input/output; handles Montgomery
  conversions and the entire round schedule.
- `poseidon_hash_3_pow5_standard(output, preimage)`: two ordinary input elements;
  starts with `[3,a,b]` and returns the final `state[1]`.
- `state_to_montgomery_3_standard`, `state_from_montgomery_3_standard`, and
  `mix_layer_3_standard`: helpers for single-round experiments. Both conversions
  and all other operations work in place. Do not convert state before calling
  the complete permutation, which already performs conversion internally.

The existing field multiplication/fifth-power helpers are reused unchanged.
The new matrix helper uses basic addition/subtraction: testing found that the
legacy assembly `f251_x_minus_c_times_y` gives an incorrect result at `x=y=0`
for multipliers 2 and 3. The original helper is preserved; this new baseline
avoids it and tests the zero-state case explicitly.

## Comparing with SD

For the same input `x`, let `y = permutation_3_pow5_standard(x)`.

- SD with `A=(1,1,1)` and `permute=false` returns exactly `y`.
- SD with any valid multipliers returns `Pi * D * y`. Raw SD output therefore
  generally differs from the baseline.
- `permutation_3_restored_sd` removes this scaling/rotation and returns exactly
  `y` for every valid configuration using SD's default fifth-power/first-element
  convention. This differs from the optional cube/last-element SD mode.

## Run the C example

From the repository root:

```sh
make -C sources clean
make -C sources CC=gcc EXTRA_CFLAGS=-DPERF_MEASUREMENT=100
gcc -Wall -Wextra -Isources examples/poseidon3_pow5_standard.c -Lsources \
  -Wl,-rpath,"$PWD/sources" -l_pos_pow5_standard -l_pos_sd \
  -o /tmp/poseidon3_demo_standard
/tmp/poseidon3_demo_standard
```

The example prints the baseline state for `[3,1,2]` and checks both identity SD
and restored SD with multipliers `(3,5,7)` and rotation enabled. It exits nonzero
on a mismatch. The final success line starts with `PASS`.

## Tests

```sh
python3 tests/test_poseidon3_pow5_standard.py
# On an x86_64 CPU with ADX:
python3 tests/test_poseidon3_pow5_standard.py \
  sources/lib_pos_pow5_standard_asm.so sources/lib_pos_sd_asm.so
```

For ISO C, clean and rebuild with `ISO_C=1`, then run the first command. Tests
compare all 91 rounds over 104 inputs to an independent integer reference,
exercise the uploaded Python hash with C parameters and SD disabled, check
SD equivalence across ten configurations, and validate invalid round arguments.
