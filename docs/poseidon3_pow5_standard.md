# Standard three-element permutation

`poseidon3-target` runs the standard path in `sources/poseidon_pow5_standard.c`.
It has no scaling, state rotation, SD configuration, or hashing wrapper.
The caller supplies all three state elements, and the function replaces them
with all three output elements. No domain tag is inserted or output selected.

```c
#include "poseidon_pow5_standard.h"

felt_t state[3] = {{3,0,0,0}, {1,0,0,0}, {2,0,0,0}};
permutation_3_pow5_standard(state);
```

The call structure follows the original C implementation:

```text
permutation_3_pow5_standard
  -> round_3_pow5_standard
       -> f251_montgomery_pow5
       -> mix_layer_3_standard
```

The permutation converts to Montgomery form, runs 4 full / 83 partial / 4 full
rounds, and converts back. Each round adds three original constants, applies
fifth power (all elements for full rounds, element 0 for partial rounds), and
multiplies by the base MDS `[[3,1,1],[1,-1,1],[1,1,-2]]`.

These are the same field, round constants, round counts, exponent and partial
S-box position as the SD path: `p=2^251+17*2^192+1`. The older `round_3_standard`
uses cube and the last-element S-box and remains unchanged. The uploaded Python
BLS12-381/Cauchy-matrix/55-partial-round parameters are not used.

The public API contains only:

- `permutation_3_pow5_standard(state)` — complete computation; ordinary input/output.
- `round_3_pow5_standard(state, round_idx, round_mode)` — one round on a Montgomery
  state, index 0..90, mode `FULL_ROUND_STANDARD` or `PARTIAL_ROUND_STANDARD`.
  Returns -1 for invalid arguments without mutation, otherwise 0.
- `mix_layer_3_standard(state)` — base MDS on a Montgomery state.

The matrix helper uses basic addition/subtraction to avoid the legacy assembly
`f251_x_minus_c_times_y` zero-operand edge case. Original functions are unchanged.
As with SD, fifth power is not invertible in this field; this is experimental.

## Run

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target poseidon3-target
./build/poseidon3-target
```

For the Makefile/shared-library workflow:

```sh
make -C sources clean
make -C sources CC=gcc EXTRA_CFLAGS=-DPERF_MEASUREMENT=100
cc -Wall -Wextra -Isources examples/poseidon3_pow5_standard.c -Lsources \
  -Wl,-rpath,"$PWD/sources" -l_pos_pow5_standard -o /tmp/poseidon3_demo_standard
/tmp/poseidon3_demo_standard
python3 tests/test_poseidon3_pow5_standard.py
# On x86_64 with ADX:
python3 tests/test_poseidon3_pow5_standard.py \
  sources/lib_pos_pow5_standard_asm.so sources/lib_pos_sd_asm.so
```

Tests compare all 91 rounds over 104 inputs against an independent integer
reference. Separate SD comparisons verify that output is `Pi * D * baseline`.
With unit multipliers and rotation disabled the raw outputs match; otherwise
only the tests undo scaling/rotation in Python for comparison. There is no
restoration function in the C API. Existing SD/hash/custom-mode convenience
wrappers and redundant conversion helpers added earlier have been removed.
