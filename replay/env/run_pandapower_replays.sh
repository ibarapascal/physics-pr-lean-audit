#!/bin/bash
# run_pandapower_replays.sh: run pandapower3016_replay.py on the merge-commit and the parent-commit
# environments (env/pp3016, env/pp_parent) with the same arguments, and yzn_base_invariance.py on both.
# The parent-commit run is used only to classify the origin of each inconsistency (introduced vs. pre-existing).
set -uo pipefail
cd "$(dirname "$0")/.."
mkdir -p outputs
env/pp3016/bin/python scripts/pandapower3016_replay.py --n 40 --seed 0 --out outputs/pandapower3016_replay.json > outputs/pandapower3016_replay.log 2>&1
echo "merge replay rc=$?"
env/pp_parent/bin/python scripts/pandapower3016_replay.py --n 40 --seed 0 --out outputs/pandapower3016_replay_parent.json > outputs/pandapower3016_replay_parent.log 2>&1
echo "parent replay rc=$?"
for e in pp3016 pp_parent; do
  env/$e/bin/python scripts/yzn_base_invariance.py $e > outputs/yzn_base_invariance_$e.json 2> outputs/yzn_base_invariance_$e.err
  echo "yzn $e rc=$?"
done
echo PANDAPOWER_REPLAYS_DONE
