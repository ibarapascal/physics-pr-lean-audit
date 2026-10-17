#!/bin/bash
# run_pp_tests.sh: check out the pandapower PR #3016 merge commit 0c6a1f0a and run, with the env/pp3016 Python
# (same commit), the test files changed by the PR (with branch coverage of the two changed source files), the
# short-circuit test directory and the three-phase power-flow tests. Single process. Set UV if uv is not on PATH.
set -uo pipefail
cd "$(dirname "$0")"
M=0c6a1f0a65273af91dd1d5b854498cd89d457939
UV=${UV:-uv}
if [ ! -d pp-src/.git ]; then
  git init -q pp-src && git -C pp-src remote add origin https://github.com/e2nIEE/pandapower.git
  git -C pp-src fetch -q --depth 2 origin "$M" && git -C pp-src checkout -q FETCH_HEAD
fi
P=$PWD/pp3016/bin/python
$P -c "import pytest, pytest_timeout, coverage" 2>/dev/null || "$UV" pip install -q --python "$P" "pytest~=9.0" pytest-timeout coverage
cd pp-src
echo "== files changed by the PR"; git diff --name-only HEAD~1 HEAD | tee ../pp_pr_files.txt
T=$(grep -E '^pandapower/test/.*\.py$' ../pp_pr_files.txt | tr '\n' ' ')
echo "== test files changed by the PR: $T"
[ -n "$T" ] && $P -m coverage run --branch --include="pandapower/pd2ppc_zero.py,pandapower/create/trafo_create.py" -m pytest -q -p no:cacheprovider -o addopts="--timeout=300" $T 2>&1 | tail -15
$P -m coverage report -m 2>&1 | tail -8
echo "== shortcircuit tests"; $P -m pytest -q -p no:cacheprovider -o addopts="--timeout=300" pandapower/test/shortcircuit 2>&1 | tail -8
echo "== three-phase power flow tests"; $P -m pytest -q -p no:cacheprovider -o addopts="--timeout=300" pandapower/test/loadflow/test_runpp_3ph.py 2>&1 | tail -8
echo PP_TESTS_DONE
