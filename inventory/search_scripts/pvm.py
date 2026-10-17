# Inspection helper: print changed lines of several candidate PRs of one repository.
import json,sys,re
# usage: pvm.py repo cap n1 n2 ...   (shows changed lines of src files, skipping bib/pyi/beads/headers-only noise)
c={ (x['repo'],x['number']):x for x in json.load(open('triaged.json'))}
repo=sys.argv[1]; cap=int(sys.argv[2])
SKIP=re.compile(r'\.(bib|pyi|jsonl)$|\.beads|\.gitattributes')
for n in sys.argv[3:]:
    x=c.get((repo,int(n)))
    if not x: print('MISSING',n); continue
    print('##',n,x['title'][:100])
    fl=[f for f in x['files'] if f['filename'] in x['src'] and not SKIP.search(f['filename'])]
    names=[f['filename'].split('/')[-1] for f in fl]
    print('  files:',len(fl),names[:6])
    budget=cap
    for f in sorted(fl,key=lambda f:-(f['additions']+f['deletions']))[:3]:
        p=f.get('patch','') or ''
        lines=[l[:170] for l in p.split('\n') if l[:1] in '+-@' and not re.match(r'^[+-]\s*(//|#|/\*|\*)',l) and l.strip() not in ('+','-')]
        s='\n'.join(lines)
        print(' --',f['filename'].split('/')[-1],f['additions'],f['deletions'])
        print(s[:budget//2]); budget-=min(len(s),budget//2)
