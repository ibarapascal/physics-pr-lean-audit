# Inspection helper: skim small unread patches of one repository (diff quick-look round).
import json,re,sys
c=json.load(open('triaged.json'))
d=set()
for l in open('decisions.txt'):
    p=l.split('|'); d.add((p[0],int(p[1])))
repo=sys.argv[1]; maxl=int(sys.argv[2])
SKIP=re.compile(r'\.(bib|pyi|jsonl|md)$|\.beads|\.gitattributes|\.clang|ref_counter|Project\.toml')
n=0
for x in c:
    if x['repo']!=repo or (repo,x['number']) in d or not x['src']: continue
    fl=[f for f in x['files'] if f['filename'] in x['src'] and not SKIP.search(f['filename'])]
    ch=sum(f['additions']+f['deletions'] for f in fl)
    if not fl or ch>maxl: continue
    lines=[]
    for f in fl:
        for l in (f.get('patch') or '').split('\n'):
            if l[:1] in '+-' and not re.match(r'^[+-]\s*(//|#|/\*|\*|"""|\'\'\')',l) and l.strip() not in('+','-') and not re.match(r'^[+-]\s*(using|import|#include|export)',l): lines.append(l[:140])
    if not lines: continue
    n+=1
    print('##',x['number'],x['title'][:70],'|',fl[0]['filename'].split('/')[-1]); print('\n'.join(lines[:6]))
print('shown',n)
