"""Extra analyses requested in review: family-clustered and leave-one-out sensitivity, reference-point rule sensitivity."""
import sys, json, collections
import numpy as np
from scipy import stats
from paths import *
import analyze as A
DEV = {'Montage_50','CyberShake_50','Sipht_60','Epigenomics_46','Inspiral_50'}
FAMS=['Montage','CyberShake','Sipht','Epigenomics','Inspiral']
fam=lambda i:i.split('_')[0]

def load(setting):
    fs=[RESULTS/f'rob_{setting}_small.jsonl',RESULTS/f'rob_{setting}_large.jsonl']
    runs=[]
    for f in fs:
        for l in open(f):
            r=json.loads(l)
            if r['type']=='run': runs.append(r)
    return runs

def gains(runs,a,b,refrule):
    pts=collections.defaultdict(list)
    for r in runs:
        if refrule(r): pts[r['instance']]+=r.get('oos',[])
    ref={i:(1.2*max(p[0] for p in v),1.2*max(p[1] for p in v)) for i,v in pts.items()}
    return ref

def run(setting,a='LIWSAP.f.r.san-plus',b='LIWSAP.f.r-plus'):
    runs=load(setting)
    out={}
    for name,(mult,cond_set) in {'paper (1.2 x worst, all 12+ conditions)':(1.2,None),
                                  '1.2 x worst, compared pair only':(1.2,{a,b}),
                                  '1.05 x worst, all conditions':(1.05,None),
                                  '2.0 x worst, all conditions':(2.0,None)}.items():
        pts=collections.defaultdict(list)
        for r in runs:
            if cond_set is None or r['cond'] in cond_set: pts[r['instance']]+=r.get('oos',[])
        ref={i:(mult*max(p[0] for p in v),mult*max(p[1] for p in v)) for i,v in pts.items()}
        hv=collections.defaultdict(lambda: collections.defaultdict(list))
        for r in runs:
            if r['cond'] in (a,b): hv[r['cond']][r['instance']].append(A.hv2d(r.get('oos',[]),ref[r['instance']]))
        inst=sorted(set(hv[a])&set(hv[b]))
        rel={i:100*(np.mean(hv[a][i])-np.mean(hv[b][i]))/np.mean(hv[b][i]) for i in inst}
        out[name]=rel
    return out

if __name__=='__main__':
    for st,lab in (('rc','RISK-AWARE (CVaR)'),('rm','EXPECTED VALUE')):
        print(f'\n===== {lab}: MF vs LIWSA+ =====')
        res=run(st)
        for name,rel in res.items():
            v=np.array(list(rel.values())); held=np.array([x for i,x in rel.items() if i not in DEV])
            x=np.array([x for i,x in rel.items()])
            print(f'{name:45s} all20 mean {v.mean():+.2f} med {np.median(v):+.2f} W/L {(v>0).sum()}/{(v<0).sum()} | held-out mean {held.mean():+.2f} med {np.median(held):+.2f} W/L {(held>0).sum()}/{(held<0).sum()}')
        rel=res['paper (1.2 x worst, all 12+ conditions)']
        # family clusters
        fm={f:np.mean([x for i,x in rel.items() if fam(i)==f]) for f in FAMS}
        print('family mean gain:', {f:round(v,2) for f,v in fm.items()})
        fv=np.array(list(fm.values())); print(f'five family means: mean {fv.mean():+.2f} median {np.median(fv):+.2f} positive {(fv>0).sum()}/5')
        # held-out by family
        fh={f:np.mean([x for i,x in rel.items() if fam(i)==f and i not in DEV]) for f in FAMS}
        print('held-out family mean gain:', {f:round(v,2) for f,v in fh.items()})
        # leave-one-family-out
        for f in FAMS:
            v=np.array([x for i,x in rel.items() if fam(i)!=f]); print(f'  without {f:12s}: mean {v.mean():+.2f} median {np.median(v):+.2f} W/L {(v>0).sum()}/{(v<0).sum()}')
        # leave-one-workflow-out range
        lo=[(i,np.mean([x for j,x in rel.items() if j!=i])) for i in rel]
        print('  leave-one-workflow-out mean range: %+.2f (omit %s) to %+.2f (omit %s)'%(min(l[1] for l in lo),min(lo,key=lambda t:t[1])[0],max(l[1] for l in lo),max(lo,key=lambda t:t[1])[0]))
        # worst/best workflow
        s=sorted(rel.items(),key=lambda t:t[1]); print('  worst 3:',[(i,round(x,2)) for i,x in s[:3]],' best 3:',[(i,round(x,2)) for i,x in s[-3:]])
        # bootstrap over families (cluster bootstrap)
        rng=np.random.RandomState(1); bs=[]
        byf={f:[x for i,x in rel.items() if fam(i)==f] for f in FAMS}
        for _ in range(5000):
            fs=rng.choice(FAMS,5); vals=[]
            for f in fs:
                vals+=list(rng.choice(byf[f],len(byf[f])))
            bs.append(np.mean(vals))
        print('  cluster-bootstrap 95%% CI of mean gain (resample families, then workflows): [%+.2f, %+.2f]'%(np.percentile(bs,2.5),np.percentile(bs,97.5)))
