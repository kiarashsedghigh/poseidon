#!/bin/sh
# Build and run the original width-3, 91-round benchmark.
set -eu
repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=${POSEIDON_BUILD_DIR:-"$repo_dir/build/old-bench"}
cmake -S "$repo_dir" -B "$build_dir" -DCMAKE_BUILD_TYPE=Release \
    -DSD_PERMUTE="${SD_PERMUTE:-ON}" -DSD_SCALE="${SD_SCALE:-ON}"
cmake --build "$build_dir" --target perf
exec "$build_dir/perf"
