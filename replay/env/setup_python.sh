#!/bin/bash
# setup_python.sh: create the Python environments used by the native replays (download and install only).
#   env/pp3016     pandapower at the PR #3016 merge commit 0c6a1f0a (pandapower 3.4.0)
#   env/pp_parent  pandapower at the parent commit be6b67a8 (same version string 3.4.0)
#   env/py2        CoolProp 8.0.0 (contains PR #2851; see CoolProp3339_replay.py for #3339) and
#                  PyBaMM 26.10.1.0 (contains PRs #5888 and #5889), plus numpy and mpmath
# The Cantera environment is built from source by build_cantera.sh; MOOSE is not installed (the two MOOSE
# replays are standard-library Python plus optional mpmath and run in any Python >= 3.8).
# Requires uv (https://github.com/astral-sh/uv); set UV to its path if it is not on PATH.
# Only the listed top-level packages are pinned; transitive dependencies resolve at install time.
set -euo pipefail
cd "$(dirname "$0")"
UV=${UV:-uv}
PP_MERGE=0c6a1f0a65273af91dd1d5b854498cd89d457939
PP_PARENT=be6b67a8562a448d538f385a43c46718956cbc31
[ -x pp3016/bin/python ] || { "$UV" venv --python 3.12 pp3016 && \
  "$UV" pip install --python pp3016/bin/python "pandapower @ git+https://github.com/e2nIEE/pandapower@${PP_MERGE}" numpy; }
pp3016/bin/python -c "import pandapower; print('pandapower (merge)', pandapower.__version__)"
[ -x pp_parent/bin/python ] || { "$UV" venv --python 3.12 pp_parent && \
  "$UV" pip install --python pp_parent/bin/python "pandapower @ git+https://github.com/e2nIEE/pandapower@${PP_PARENT}" numpy; }
pp_parent/bin/python -c "import pandapower; print('pandapower (parent)', pandapower.__version__)"
[ -x py2/bin/python ] || { "$UV" venv --python 3.12 py2 && \
  "$UV" pip install --python py2/bin/python "CoolProp==8.0.0" "pybamm==26.10.1.0" numpy mpmath; }
py2/bin/python -c "import CoolProp, pybamm; print('CoolProp', CoolProp.__version__); print('pybamm', pybamm.__version__)"
echo ENV_OK
