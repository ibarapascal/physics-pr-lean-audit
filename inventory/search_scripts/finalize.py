# Inventory pass 1, stage 4: merge per-PR decisions (decisions.txt, read-the-patch judgments) with automatic exclusions into pass1_classify.csv.
import json, re, csv, os
from collections import Counter, defaultdict
OUT=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')  # inventory/
c=json.load(open('triaged.json'))
dec={}
for l in open('decisions.txt'):
    l=l.rstrip('\n')
    if not l: continue
    p=l.split('|',4)
    dec[(p[0],int(p[1]))]=(p[2],p[3],p[4])
FRAME={'FEniCS/dolfinx','firedrakeproject/firedrake','SciML/ModelingToolkit.jl','SciML/ModelingToolkitStandardLibrary.jl','devitocodes/devito','OpenMDAO/OpenMDAO'}
# The Chinese labels/patterns below are matched against, or written into, the original Chinese records (decisions.txt and the classification tables); their English counterparts are in the *_en columns / decisions.en.txt.
CATS=[('docs?|documentation|docstring|readme','文档'),('ci\\b|build|cmake|packag|wheel|dep|bump|release|version|compat|dockerfile|spack|install|licen|beads|blame','构建/CI/依赖/发布'),
 ('refactor|style|format|clang|lint|modernize|cleanup|tidy|ambig|piracy|qa\\b|runic|cppcheck|codeql|warning','重构/格式/静态检查'),
 ('perf|speed|faster|optimi|alloc|memory|cache','性能/内存'),
 ('python|nanobind|binding|wrapper|interface|api|stub|serialis|serializ|json|netcdf|zarr|output|writer|reader|mesh|grid|plot|gui|mathcad|bpx|sbtl|svd','接口/IO/网格/绘图/表面模型'),
 ('feat|add|support|new|implement','新功能/基础设施'),('fix|correct|handle|restore|robust|guard|validate','缺陷修补（非物理量路径）')]
def cat(t):
    for rx,n in CATS:
        if re.search(rx,t,re.I): return n
    return '实现层改动'
rows=[]; stat=defaultdict(lambda:[0,0,0,0,0]) # cand, yes, patch_read, auto_nosrc, triage
inv=Counter()
# The Chinese labels/patterns below are matched against, or written into, the original Chinese records (decisions.txt and the classification tables); their English counterparts are in the *_en columns / decisions.en.txt.
for x in c:
    k=(x['repo'],x['number'])
    s=stat[x['repo']]; s[0]+=1
    if k in dec:
        y,ic,r=dec[k]; s[2]+=1
        rows.append([x['repo'],x['number'],x['url'],y,ic,'[读patch] '+r])
        if y=='yes': s[1]+=1; inv[ic]+=1
    elif not x['src']:
        s[3]+=1
        rows.append([x['repo'],x['number'],x['url'],'no','无',f'[自动排除] {len(x["files"])} 个改动文件全部属文档/测试/CI/构建/示例/配置，无物理源码'])
    else:
        s[4]+=1
        names=[f.split('/')[-1] for f in x['src'][:3]]
        extra='；该库为数值/建模框架，库内无物理方程' if x['repo'] in FRAME else ''
        rows.append([x['repo'],x['number'],x['url'],'no','无',f'[路径+标题排除，未逐行读patch] 标题「{x["title"][:70]}」属{cat(x["title"])}，源文件 {", ".join(names)} 等共 {len(x["src"])} 个，未见本构/守恒/单位/符号/热力学函数改动迹象{extra}'])
with open(os.path.join(OUT,'pass1_classify.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['repo','number','url','physics_change','invariant_class','reason']); w.writerows(rows)
tot=[sum(v[i] for v in stat.values()) for i in range(5)]
print('repo cand yes read autoNoSrc triage')
for k,v in sorted(stat.items(), key=lambda kv:-kv[1][0]): print(k,*v)
print('TOTAL',tot); print(inv)
print('yes repos',sum(1 for v in stat.values() if v[1]>0))
json.dump({k:v for k,v in stat.items()}, open('stat.json','w'))
