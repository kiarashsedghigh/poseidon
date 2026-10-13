# Build and run the C examples with CMake

From the repository root, with CMake 3.20+ and a C99 compiler installed:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target poseidon3-target poseidon3-sd

./build/poseidon3-target
./build/poseidon3-sd
./build/poseidon3-sd 0 3 5 7

ctest --test-dir build --output-on-failure
```

- `poseidon3-target` runs the non-SD fifth-power example and checks identity SD
  and restored SD outputs against the baseline.
- `poseidon3-sd` runs the SD example. Optional arguments are `PERMUTE A0 A1 A2`;
  defaults are `1 1 1 2`. The example input is `[3,1,2]` in both executables.

Both use the C field/constants/base matrix and fifth-power, first-element
partial S-box. The earlier cube-based functions are unchanged. These examples
remain experimental: fifth power is not invertible over the C field.

CMake builds static implementation libraries directly from source, so no prior
Makefile build or shared-library search path is needed. It uses the C backend
with automatic `__int128` selection; enable `-DPOSEIDON_ISO_C=ON` for ISO C.
The existing Makefile remains available for assembly builds.

In CLion, open the root `CMakeLists.txt`, reload CMake, and select either
`poseidon3-target` or `poseidon3-sd` as the run target. For the SD executable,
set program arguments in its Run/Debug configuration to change the controls.
With multi-configuration generators, build with `--config Release` and find
the executables in `build/Release`; use `ctest -C Release` for that build.
