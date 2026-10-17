#!/bin/bash
# run_julia_replays.sh: run the four Julia replays serially, single-threaded, with the environments from
# setup_julia.sh. stdout goes to outputs/<name>, stderr (progress and log lines) to outputs/<name>.err.
# ClimaAtmos4637 is run last because it needs the most memory (more than 8 GiB peak).
cd "$(dirname "$0")/.."
export JULIA_DEPOT_PATH=$PWD/env/depot JULIA_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=2 JULIA_PKG_PRECOMPILE_AUTO=0
J=$PWD/env/julia-1.11.7/bin/julia
mkdir -p outputs
run() { local proj=$1 f=$2 out=$3; echo "== $f"; $J --project=env/$proj scripts/$f > outputs/$out 2> outputs/$out.err; echo "exit $? $f"; }
run ca    climaatmos4822_diff.jl     ClimaAtmos4822_diff.out
run ocean Oceananigans6083_replay.jl Oceananigans6083_replay.json
run ocean Oceananigans6024_replay.jl Oceananigans6024_replay.json
CA_COMMIT=ba2bf1225fc48058d709853c065b1db21a412204 run ca ClimaAtmos4637_replay.jl ClimaAtmos4637_replay.json
echo JULIA_REPLAYS_DONE
