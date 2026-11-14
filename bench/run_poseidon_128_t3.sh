#!/bin/sh
# Run one profile; optional arguments: iterations, warmup.
set -eu
script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec "$script_dir/run_profiles.sh" poseidon-128-t3 "$@"
