#!/bin/sh
# run_full.sh: full judge run inside the OpenModelica container (controlled/ mounted at /work): 5 tasks x
# (pre + controls + 10 samples), stages (i) regression on all regression examples and (ii) N = 10^4 spot checks.
# Output: results/full.jsonl (records) and results/full.stdout (one summary line per variant).
set -u
cd /work
rm -rf work
mkdir -p results
python3 judge/judge_b.py results/full.jsonl --variants pre,ctrl,gen --n 10000 --par 2 --ppar 4 > results/full.stdout 2>&1
echo "rc=$?"; tail -70 results/full.stdout
