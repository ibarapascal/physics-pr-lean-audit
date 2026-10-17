#!/usr/bin/env python3
"""summarize_logs.py: rebuild logs/summary.txt and logs/axioms_summary.txt from the per-module logs.

Inputs (relative to lean/):
  PhysAudit/<module>.lean     source; lines containing sorry/admit/native_decide are counted
  logs/<module>.log           output of `lake build PhysAudit.<module>`
  logs/axioms_<module>.txt    output of `#print axioms` for every theorem/lemma of the module
Outputs:
  logs/summary.txt            one line per module: build status, module compile time reported by lake,
                              warning count, sorry-like source lines
  logs/axioms_summary.txt     one line per module: theorems checked, axioms used, nonstandard axioms
                              (anything other than propext, Classical.choice, Quot.sound), error count
"""
import os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STANDARD = {'propext', 'Classical.choice', 'Quot.sound'}


def modules():
    src = os.path.join(ROOT, 'PhysAudit')
    return sorted(f[:-5] for f in os.listdir(src) if f.endswith('.lean'))


def build_line(m):
    log = open(os.path.join(ROOT, 'logs', f'{m}.log'), encoding='utf-8').read()
    ok = 'Build completed successfully' in log
    t = re.search(r'Built PhysAudit\.%s \(([0-9.]+m?s)\)' % re.escape(m), log)
    replayed = re.search(r'Replayed PhysAudit\.%s\b' % re.escape(m), log)
    when = t.group(1) if t else ('replayed from cache' if replayed else 'n/a')
    warns = len(re.findall(r'^warning:', log, re.M))
    errs = len(re.findall(r'^error:', log, re.M))
    src = open(os.path.join(ROOT, 'PhysAudit', f'{m}.lean'), encoding='utf-8').read().splitlines()
    sorry = sum(1 for l in src if re.search(r'\b(sorry|admit|native_decide)\b', l))
    status = 'ok' if ok and errs == 0 else 'FAILED'
    return f'PhysAudit.{m} build={status} module_time={when} warnings={warns} errors={errs} sorry_like_lines={sorry}'


def axioms_line(m):
    txt = open(os.path.join(ROOT, 'logs', f'axioms_{m}.txt'), encoding='utf-8').read()
    deps = re.findall(r"'([^']+(?:'[^'\s]*)?)' depends on axioms: \[([^\]]*)\]", txt, re.S)
    none = re.findall(r"'([^']+)' does not depend on any axioms", txt)
    used, nonstd = set(), set()
    for _, lst in deps:
        for a in (x.strip() for x in lst.split(',')):
            if a:
                used.add(a)
                if a not in STANDARD:
                    nonstd.add(a)
    errs = len(re.findall(r'\berror\b', txt))
    n = len(deps) + len(none)
    return (f'{m} theorems={n} axioms_used=[{" ".join(sorted(used))}] '
            f'nonstandard_axioms=[{" ".join(sorted(nonstd))}] errors={errs}'), nonstd, errs


def main():
    mods = modules()
    with open(os.path.join(ROOT, 'logs', 'summary.txt'), 'w', encoding='utf-8') as f:
        for m in mods:
            f.write(build_line(m) + '\n')
    bad = 0
    with open(os.path.join(ROOT, 'logs', 'axioms_summary.txt'), 'w', encoding='utf-8') as f:
        for m in mods:
            line, nonstd, errs = axioms_line(m)
            bad += bool(nonstd) or errs > 0
            f.write(line + '\n')
        f.write(f'modules={len(mods)} modules_with_nonstandard_axioms_or_errors={bad}\n')
    print(open(os.path.join(ROOT, 'logs', 'summary.txt')).read())
    print(open(os.path.join(ROOT, 'logs', 'axioms_summary.txt')).read())
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
