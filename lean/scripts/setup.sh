#!/bin/bash
# setup.sh: fetch the pinned Mathlib build cache and build the PhysAudit library.
# Requires elan (https://github.com/leanprover/elan); the toolchain is pinned in lean-toolchain
# (leanprover/lean4:v4.34.1) and Mathlib is pinned in lake-manifest.json (tag v4.34.1).
# Equivalent manual steps, run from the lean/ directory:
#   lake exe cache get     # download prebuilt Mathlib .olean files (avoids compiling Mathlib)
#   lake build             # build the root module PhysAudit (imports all 12 modules)
set -euo pipefail
cd "$(dirname "$0")/.."
lean --version
lake exe cache get
lake build 2>&1 | tee logs/root_build.log
