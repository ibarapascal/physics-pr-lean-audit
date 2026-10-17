# Inspection helper: print changed lines of one candidate PR from triaged.json.
import json,sys,re
# usage: pv.py repo num [maxchars] [filefilter-regex]
c=json.load(open('triaged.json'))
repo,num=sys.argv[1],int(sys.argv[2]); mx=int(sys.argv[3]) if len(sys.argv)>3 else 3000
ff=re.compile(sys.argv[4]) if len(sys.argv)>4 else None
for x in c:
    if x['repo']==repo and x['number']==num:
        print('##',x['title'],'|',x['url'])
        tot=0
        for f in x['files']:
            if f['filename'] not in x['src']: continue
            if ff and not ff.search(f['filename']): continue
            p=f.get('patch','') or ''
            # keep only changed lines + hunk headers
            lines=[l for l in p.split('\n') if l[:1] in '+-@' ]
            s='\n'.join(lines)
            print('--',f['filename'],f['status'],f['additions'],f['deletions'])
            print(s[:mx//max(1,min(len(x['src']),4))] )
