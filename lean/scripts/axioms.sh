#!/bin/bash
# axioms.sh: run `#print axioms` on every theorem/lemma of the given modules (default: every module under
# PhysAudit/), write the raw output to logs/axioms_<module>.txt and a summary to logs/axioms_summary.txt.
# propext, Classical.choice and Quot.sound are Lean's standard axioms; any other axiom (including sorryAx)
# is reported as nonstandard. The summary is produced by scripts/summarize_logs.py, which parses the raw
# output (Lean wraps long axiom lists over several lines).
set -uo pipefail
cd "$(dirname "$0")/.."
mkdir -p logs/tmp
mods=("$@"); [ ${#mods[@]} -eq 0 ] && mods=($(ls PhysAudit/*.lean | xargs -n1 basename | sed 's/\.lean$//' | sort))
for m in "${mods[@]}"; do
  f=PhysAudit/$m.lean; t=logs/tmp/ax_$m.lean
  echo "import PhysAudit.$m" > "$t"
  python3 - "$f" >> "$t" <<'PY'
import re, sys
ns, out = [], []
for line in open(sys.argv[1], encoding='utf-8'):
    m = re.match(r'\s*namespace\s+(\S+)', line)
    if m: ns.append(m.group(1)); continue
    m = re.match(r'\s*end\s+(\S+)\s*$', line)
    if m and ns and ns[-1] == m.group(1): ns.pop(); continue
    m = re.match(r'\s*(?:private\s+)?(?:theorem|lemma)\s+([^\s(:{\[]+)', line)
    if m: out.append('.'.join(ns + [m.group(1)]))
for n in out: print(f'#print axioms {n}')
PY
  lake env lean "$t" > "logs/axioms_$m.txt" 2>&1
done
rm -rf logs/tmp
python3 scripts/summarize_logs.py
