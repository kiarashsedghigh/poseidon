#!/bin/sh
# Usage: bench/run_profiles.sh [all|PROFILE [ITERATIONS [WARMUP]]]
set -eu
repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build_dir=${POSEIDON_BUILD_DIR:-"$repo_dir/build/profiles-bench"}
cmake -S "$repo_dir" -B "$build_dir" -DCMAKE_BUILD_TYPE=Release \
    -DSD_PERMUTE="${SD_PERMUTE:-ON}" -DSD_SCALE="${SD_SCALE:-ON}"
cmake --build "$build_dir" --target profiles-perf
exec "$build_dir/profiles-perf" "$@"
