# Inventory pass 1, stage 3: mark source files (non-doc/test/CI/build) per candidate; writes triaged.json.
import json, re, os
from collections import Counter, defaultdict
D='.'
c=json.load(open('cands_full.json'))
NONSRC = re.compile(r'(^|/)(docs?|documentation|test|tests|testing|unit_test|unit_tests|\.github|\.buildkite|\.ci|ci|benchmarks?|bench|examples?|tutorials?|dev|scripts?|tools|cmake|CMake|\.vscode|assets|figures|www|web|wrappers/(Java|Javascript|Csharp|Octave|MATLAB|Excel|R|Delphi|EES|LabVIEW|Scilab|Maple|Mathematica|PCtools|SMath|Labview|Julia|Python/CoolProp/(tests)?))/|(^|/)test_[^/]*\.(py|jl)$|_test\.(py|jl|C|cpp|h)$|\.(md|rst|txt|yml|yaml|toml|cfg|lock|png|jpg|svg|json|cff|ipynb|gitignore|in|cmake)$|(^|/)(Project\.toml|Manifest\.toml|CHANGELOG[^/]*|LICENSE[^/]*|Makefile|CMakeLists\.txt|requirements[^/]*|setup\.py|pyproject\.toml|CLAUDE\.md|\.gitignore|\.gitmodules)$|/(runtests|tests?)\.jl$', re.I)
res=[]
for x in c:
    names=[f['filename'] for f in x['files']]
    src=[n for n in names if (not NONSRC.search(n)) or re.match(r'dev/(fluids|mixtures|cubics|pcsaft|incompressible)',n)]
    x['src']=src
    res.append(x)
by=defaultdict(lambda:[0,0])
for x in res:
    by[x['repo']][0]+=1
    by[x['repo']][1]+= 1 if x['src'] else 0
for k,v in sorted(by.items(), key=lambda kv:-kv[1][1]): print(k,v)
print(sum(v[1] for v in by.values()))
json.dump(res, open('triaged.json','w'))
