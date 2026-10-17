# Similarity of Codex sample patches to the historical correct (and wrong) control diffs
import json,re,difflib,glob,os
B=os.path.dirname(os.path.abspath(__file__))+'/'   # run from anywhere; expects judge/controls/ and samples/ next to this file
def added(diff, only_files=None):
    out=[];cur=None
    for l in diff.split('\n'):
        m=re.match(r'^\+\+\+ (?:b/)?(\S+)',l)
        if m: cur=m.group(1); continue
        if l.startswith('--- '): continue
        if l.startswith('+'):
            s=re.sub(r'\s+','',l[1:])
            s=re.sub(r'//.*','',s)
            if s and not s.startswith('"') and (only_files is None or any(cur.endswith(f) for f in only_files)):
                out.append(s)
    return out
def files_of(diff):
    return [re.sub(r'^(?:b/)','',m) for m in re.findall(r'^\+\+\+ (\S+)',diff,re.M)]
def sim(a,b):
    A='\n'.join(a); Bs='\n'.join(b)
    return difflib.SequenceMatcher(None,A,Bs,autojunk=False).ratio()
for task,ctrls in [('buildings_3027',['correct','wrong']),('msl_1335',['correct','wrong']),('msl_1530',['correct','wrong']),('buildings_2725',['wrong']),('buildings_4775',['wrong'])]:
    for c in ctrls:
        cd=open(B+f'judge/controls/{task}_{c}.diff').read()
        cf=[f.split('/')[-1] for f in files_of(cd)]
        ca=added(cd,cf)
        res=[]
        seen={}
        for p in sorted(glob.glob(B+f'samples/{task}/codex/s*.json')):
            d=json.load(open(p))
            dd=d.get('diff') or ''
            if not dd: continue
            sa=added(dd,cf)
            exact=set(sa)&set(ca)
            res.append((os.path.basename(p)[:3], round(sim(sa,ca),3), len(exact), len(ca)))
        best=max(res,key=lambda x:x[1]) if res else None
        print(task,c,'ctrl_files',cf,'n_ctrl_added',len(ca),'best',best,'all',[r[1] for r in res])
