#!/bin/sh
# Run all profile scripts; optional arguments: iterations, warmup.
set -eu
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
for profile in 80_t3 80_t5 128_t3 128_t5 256_t6 256_t10 x; do
    "$script_dir/run_poseidon_$profile.sh" "$@"
done
