#!/bin/bash
# build.sh: build PhysAudit module by module (arguments: module short names; default: every .lean file
# under PhysAudit/). Writes one log per module, logs/<module>.log, and a summary logs/summary.txt with
# the exit code, wall time and the number of source lines containing sorry/admit/native_decide.
set -uo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs; : > logs/summary.txt
mods=("$@"); [ ${#mods[@]} -eq 0 ] && mods=($(ls PhysAudit/*.lean | xargs -n1 basename | sed 's/\.lean$//' | sort))
for m in "${mods[@]}"; do
  t0=$(date +%s); lake build "PhysAudit.$m" > "logs/$m.log" 2>&1; rc=$?
  ns=$(grep -cE "sorry|admit|native_decide" "PhysAudit/$m.lean" || true)
  echo "PhysAudit.$m rc=$rc $(( $(date +%s)-t0 ))s sorry_like_lines=$ns" | tee -a logs/summary.txt
done
echo BUILD_DONE | tee -a logs/summary.txt
