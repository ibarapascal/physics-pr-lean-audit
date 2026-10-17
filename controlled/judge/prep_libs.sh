#!/bin/sh
# prep_libs.sh: download the full library trees at the pre-change commit of every task and the MSL versions they
# depend on (GitHub codeload tarballs) into controlled/libs/, print the tarball SHA-256, and unpack. Download and
# unpack only; nothing is compiled.
set -eu
cd "$(dirname "$0")/.."            # controlled/
mkdir -p libs/_dl && cd libs
get() { # name repo ref
  [ -f "_dl/$1.tar.gz" ] || curl -fsSL --retry 3 -o "_dl/$1.tar.gz" "https://codeload.github.com/$2/tar.gz/$3"
  shasum -a 256 "_dl/$1.tar.gz"
}
unpack() { # name keep...
  n=$1; shift; rm -rf "$n"; mkdir -p "$n"
  top=$(tar -tzf "_dl/$n.tar.gz" | head -1 | cut -d/ -f1)
  for k in "$@"; do tar -xzf "_dl/$n.tar.gz" -C "$n" --strip-components 1 \
      --exclude '*/Resources/Images/*' --exclude '*/Resources/www/*' --exclude '*/Resources/Library/*' \
      --exclude '*/Resources/Data/ThermalZones/EnergyPlus*' "$top/$k" || echo "MISSING $n/$k"; done
  du -sh "$n"
}
get mb_309c814f lbl-srg/modelica-buildings 309c814f501b1e7328ce086d29be1f72e5019ed1
get mb_dff251c2 lbl-srg/modelica-buildings dff251c2a154412386ce4dab02a1631766b4ad23
get mb_b83b2278 lbl-srg/modelica-buildings b83b2278b31f61dac4ae8e55b7e9b48cd2b3e049
get msl_v323    modelica/ModelicaStandardLibrary v3.2.3
get msl_a3041258 modelica/ModelicaStandardLibrary a3041258138620322ff5719d82c53f6b3d3f78d3
get msl_7c8ccafa modelica/ModelicaStandardLibrary 7c8ccafa0ba615ab06bc4432a3a80b2aed0ecbfa
for n in mb_309c814f mb_dff251c2 mb_b83b2278; do unpack $n Buildings; done
for n in msl_v323 msl_a3041258 msl_7c8ccafa; do unpack $n Modelica ModelicaServices Complex.mo ModelicaTest; done
# Environment patch: MSL 3.2.3 on GitHub ships the reference ModelicaServices, whose loadResource does not resolve
# modelica:// URIs (Buildings then fails to read weather files); use the OpenModelica implementation instead.
# By MSL convention ModelicaServices is provided by the tool; the patch is independent of the libraries under test
# and identical for all variants.
sed -i.orig 's|fileReference := Modelica.Utilities.Files.fullPathName(uri);|fileReference := OpenModelica.Scripting.uriToFilename(uri);|' msl_v323/ModelicaServices/package.mo
grep -c "OpenModelica.Scripting.uriToFilename" msl_v323/ModelicaServices/package.mo
ls -la */ | head -60
rm -rf _dl/*.tar.gz.keep
echo PREP_DONE
