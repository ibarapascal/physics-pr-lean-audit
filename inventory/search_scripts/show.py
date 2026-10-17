# Inspection helper: list candidates of one repository with source-file change sizes.
import json,sys
c=json.load(open('triaged.json'))
repo=sys.argv[1]
for x in c:
    if x['repo']!=repo or not x['src']: continue
    fs={f['filename']:f for f in x['files']}
    a=sum(fs[n]['additions'] for n in x['src']); d=sum(fs[n]['deletions'] for n in x['src'])
    print(f"{x['number']} {x['title'][:90]} | {len(x['src'])}f +{a}/-{d} | "+', '.join(n.split('/')[-1] if len(x['src'])>2 else n for n in x['src'][:3]))
