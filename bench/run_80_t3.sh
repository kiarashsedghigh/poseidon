#!/bin/sh
# Build and run only poseidon-80-t3 in the earlier benchmark's format.
set -eu
repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=${POSEIDON_BUILD_DIR:-"$repo_dir/build/bench-80-t3"}
cmake -S "$repo_dir" -B "$build_dir" -DCMAKE_BUILD_TYPE=Release \
    -DSD_PERMUTE="${SD_PERMUTE:-ON}" -DSD_SCALE="${SD_SCALE:-ON}"
cmake --build "$build_dir" --target perf-80-t3
exec "$build_dir/perf-80-t3"
