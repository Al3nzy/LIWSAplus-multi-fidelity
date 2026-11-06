import sys,json,collections
import numpy as np
from scipy import stats
import analyze as A
from paths import *
DEV={'Montage_50','CyberShake_50','Sipht_60','Epigenomics_46','Inspiral_50'}
fam=lambda i:i.split('_')[0]
def go(st):
    runs=[]
    for f in (f'ideal_{st}_small.jsonl',f'ideal_{st}_large.jsonl'):
        runs+=[json.loads(l) for l in open(RESULTS/f)]
    runs=[r for r in runs if r['type']=='run']
    pts=collections.defaultdict(list)
    for r in runs: pts[r['instance']]+=r['oos']
    ref={i:(1.2*max(p[0] for p in v),1.2*max(p[1] for p in v)) for i,v in pts.items()}
    hv=collections.defaultdict(lambda:collections.defaultdict(list)); tm=collections.defaultdict(lambda:collections.defaultdict(list))
    for r in runs:
        hv[r['cond']][r['instance']].append(A.hv2d(r['oos'],ref[r['instance']]))
    return hv
def cmp(hv,a,b,sub=None):
    ia=sorted(set(hv[a])&set(hv[b]))
    if sub: ia=[i for i in ia if sub(i)]
    x=np.array([np.mean(hv[a][i]) for i in ia]);y=np.array([np.mean(hv[b][i]) for i in ia])
    rel=100*(x-y)/y
    return dict(n=len(ia),mean=rel.mean(),med=np.median(rel),W=int((rel>0).sum()),L=int((rel<0).sum()),p=stats.wilcoxon(x,y).pvalue,rel=dict(zip(ia,rel)))
N={'base':'LIWSAP.f.r-plus','faith':'LIWSAP.f.r.san-plus','ideal':'LIWSAP.f.r.sai-plus','rand':'LIWSAP.f.r.sar-plus'}
if __name__=='__main__':
    for st in ('rc','rm'):
        hv=go(st); print('=====',st,{k:len(hv[v]) for k,v in N.items()})
        for a,b in [('faith','base'),('ideal','base'),('rand','base'),('faith','ideal'),('faith','rand'),('ideal','rand')]:
            for lab,sub in (('all',None),('held',lambda i:i not in DEV)):
                r=cmp(hv,N[a],N[b],sub); print(f'{a:5s} vs {b:5s} {lab:4s} n={r["n"]} mean {r["mean"]:+.2f} med {r["med"]:+.2f} W/L {r["W"]}/{r["L"]} p={r["p"]:.3f}')
        r=cmp(hv,N['faith'],N['ideal'])['rel']
        for f in ['Montage','CyberShake','Sipht','Epigenomics','Inspiral']:
            v=[x for i,x in r.items() if fam(i)==f]; print('  faith vs ideal',f,f'{np.mean(v):+.2f}')
