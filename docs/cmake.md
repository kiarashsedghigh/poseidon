# Build and run the C examples with CMake

From the repository root, with CMake 3.20+ and a C99 compiler installed:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target poseidon3-target poseidon3-sd

./build/poseidon3-target
./build/poseidon3-sd
./build/poseidon3-sd 3 5 7

```

- `poseidon3-target` runs only the non-SD fifth-power example. It calls
  `permutation_3_pow5_standard`, which executes `round_3_pow5_standard` for
  each round. It has no multiplier configuration, state rotation, or SD code.
- `poseidon3-sd` runs the SD example. Optional arguments are `A0 A1 A2`;
  defaults are `1 1 2`. Rotation and scaling are compile-time CMake options
  `SD_PERMUTE` and `SD_SCALE`, both default `ON`. The example input is `[3,1,2]` in both executables.

Both use the C field/constants/base matrix and fifth-power, first-element
partial S-box. The earlier cube-based functions are unchanged. These examples
remain experimental: fifth power is not invertible over the C field.

CMake defines these two example targets plus `perf` for the benchmark and builds
their sources directly; there are no helper library or CTest targets. No prior Makefile build
or shared-library search path is needed. It uses the C backend
with automatic `__int128` selection; enable `-DPOSEIDON_ISO_C=ON` for ISO C.
The existing Makefile remains available for assembly builds.

In CLion, open the root `CMakeLists.txt`, reload CMake, and select either
`poseidon3-target` or `poseidon3-sd` as the run target. For the SD executable,
set program arguments to change multiplier values. Set `-DSD_PERMUTE=OFF`
and/or `-DSD_SCALE=OFF` in the CMake profile to disable those operations, then
reload and rebuild.
With multi-configuration generators, build with `--config Release` and find
the executables in `build/Release`.

SD comparison tests remain separate in `tests/test_poseidon3_pow5_standard.py`;
see `poseidon3_pow5_standard.md` for the Makefile build and Python test commands.
After pulling the simplified configuration, reload CMake in CLion. If old
targets still appear, use **Reset Cache and Reload Project**.

For the timing benchmark, edit `bench/config.h`, build target `perf`, and run
`./build/perf`. See `../bench/README.md` for configuration and measurement details.
