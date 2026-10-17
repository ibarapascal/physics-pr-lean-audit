#!/bin/bash
# setup_julia.sh: install Julia 1.11.7 from the official tarball into env/julia-1.11.7 and create two
# project environments with their own depot (env/depot):
#   env/ca     ClimaAtmos.jl at ba2bf122 (merge commit of PR #4822; contains the PR #4637 merge a8efce5a);
#              ClimaComms is added to the project (already in the resolved manifest)
#   env/ocean  Oceananigans.jl at 24de7026 (merge commit of PR #6024; contains the PR #6083 merge 1734ff9e);
#              Adapt is added to the project (already in the resolved manifest)
# JULIA_URL selects the tarball for the platform (default: macOS aarch64). The resolved Manifest.toml files
# of the original runs were not archived; package versions resolve at install time.
set -euo pipefail
cd "$(dirname "$0")"
JV=1.11.7
JULIA_URL=${JULIA_URL:-https://julialang-s3.julialang.org/bin/mac/aarch64/1.11/julia-$JV-macaarch64.tar.gz}
if [ ! -x julia-$JV/bin/julia ]; then
  curl -fsSL -o julia.tar.gz "$JULIA_URL"
  tar -xzf julia.tar.gz
fi
J=$PWD/julia-$JV/bin/julia
export JULIA_DEPOT_PATH=$PWD/depot JULIA_NUM_PRECOMPILE_TASKS=${JULIA_NUM_PRECOMPILE_TASKS:-4} JULIA_PKG_PRECOMPILE_AUTO=0 JULIA_NUM_THREADS=1
$J --version
mkdir -p ca ocean
$J --project=ca -e 'using Pkg; Pkg.add(url="https://github.com/CliMA/ClimaAtmos.jl", rev="ba2bf1225fc48058d709853c065b1db21a412204"); Pkg.add("ClimaComms"; preserve=Pkg.PRESERVE_ALL); Pkg.precompile(); using ClimaAtmos; println("ClimaAtmos ", pkgversion(ClimaAtmos))'
$J --project=ocean -e 'using Pkg; Pkg.add(url="https://github.com/CliMA/Oceananigans.jl", rev="24de7026893af38df065920fbbf80f05024fd1b9"); Pkg.add("Adapt"; preserve=Pkg.PRESERVE_ALL); Pkg.precompile(); using Oceananigans; println("Oceananigans ", pkgversion(Oceananigans))'
echo JULIA_ENV_OK
