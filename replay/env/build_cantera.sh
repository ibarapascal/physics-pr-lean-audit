#!/bin/bash
# build_cantera.sh: build the Cantera Python package from the PR #2107 merge commit d58a2326 (reports
# version 4.0.0a1). Build dependencies come from conda-forge via micromamba; set MICROMAMBA to its path if
# it is not on PATH and JOBS to the build parallelism (default 4). Creates env/ctb (conda prefix) and
# env/cantera-src (source checkout).
set -euo pipefail
cd "$(dirname "$0")"
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$PWD/mamba}
MM=${MICROMAMBA:-micromamba}
JOBS=${JOBS:-4}
SHA=d58a23265414e1789536cf5332d9fb29de095a43
[ -x ctb/bin/python ] || "$MM" create -y -p ./ctb -c conda-forge python=3.12 scons cython numpy ruamel.yaml pip \
  boost-cpp eigen fmt yaml-cpp sundials highfive hdf5 pytest packaging setuptools wheel compilers doxygen jinja2
# The Python interface build needs doxygen XML (Cantera.tag) and jinja2 for the generated clib sources.
ctb/bin/python -c "import jinja2" 2>/dev/null && [ -x ctb/bin/doxygen ] || "$MM" install -y -p ./ctb -c conda-forge doxygen jinja2
# Source archives do not include submodules: fetch the merge commit with git and only the
# data/example_data submodule (other dependencies come from the conda prefix).
if [ ! -d cantera-src/.git ]; then
  rm -rf cantera-src
  git init -q cantera-src && git -C cantera-src remote add origin https://github.com/Cantera/cantera.git
  git -C cantera-src fetch -q --depth 1 origin "$SHA" && git -C cantera-src checkout -q FETCH_HEAD
  git -C cantera-src submodule update --init --depth 1 data/example_data
fi
cd cantera-src
export PATH=$PWD/../ctb/bin:$PATH
P=$PWD/../ctb
scons build -j"$JOBS" python_package=y f90_interface=n doxygen_docs=n sphinx_docs=n googletest=none \
  system_eigen=y system_fmt=y system_yamlcpp=y system_sundials=y system_highfive=n hdf_support=n \
  extra_inc_dirs=$P/include extra_lib_dirs=$P/lib boost_inc_dir=$P/include prefix=$P python_cmd=$P/bin/python 2>&1 | tail -40
$P/bin/python -m pip install --no-deps --no-build-isolation build/python 2>&1 | tail -5 || true
cd ..
ctb/bin/python -c "import cantera as ct; print('cantera', ct.__version__, getattr(ct,'__git_commit__',''))"
echo BUILD_CT_OK
