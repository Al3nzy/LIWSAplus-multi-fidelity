import sys,json,collections,os
import numpy as np
from scipy import stats
import analyze as A
from paths import *
DEV={'Montage_50','CyberShake_50','Sipht_60','Epigenomics_46','Inspiral_50'}
def hvs(fn):
    p=RESULTS/fn
    if not p.exists(): return None
    runs=[json.loads(l) for l in open(p)]; runs=[r for r in runs if r['type']=='run']
    pts=collections.defaultdict(list)
    for r in runs: pts[r['instance']]+=r['oos']
    ref={i:(1.2*max(q[0] for q in v),1.2*max(q[1] for q in v)) for i,v in pts.items()}
    hv=collections.defaultdict(lambda:collections.defaultdict(list))
    for r in runs: hv[r['cond']][r['instance']].append(A.hv2d(r['oos'],ref[r['instance']]))
    return hv
def cmp(hv,a,b,sub=None):
    ia=sorted(set(hv[a])&set(hv[b]))
    if sub: ia=[i for i in ia if sub(i)]
    x=np.array([np.mean(hv[a][i]) for i in ia]);y=np.array([np.mean(hv[b][i]) for i in ia])
    rel=100*(x-y)/y
    return f'n={len(ia)} mean {rel.mean():+.2f} med {np.median(rel):+.2f} W/L {(rel>0).sum()}/{(rel<0).sum()} p={stats.wilcoxon(x,y).pvalue:.3f}', dict(zip(ia,rel))
B='LIWSAP.f.r-plus';M='LIWSAP.f.r.san-plus';Q='LIWSAP.f.r.saq-plus'
def report():
    for fn in ['rep1_rc_small','rep1_rm_small','rep2_rc_small','rep2_rm_small','quota_rc_small','quota_rm_small']:
        hv=hvs(fn+'.jsonl')
        if hv is None: print(fn,'missing'); continue
        print('==',fn,{k[-12:]:len(v) for k,v in hv.items()})
        if len(hv[B])<15 or len(hv[M])<15: print(' incomplete'); continue
        print(' MF vs base all',cmp(hv,M,B)[0]); print(' MF vs base held',cmp(hv,M,B,lambda i:i not in DEV)[0])
        r=cmp(hv,M,B)[1]; print(' Sipht_30',round(r.get('Sipht_30',float('nan')),2))
        if Q in hv and len(hv[Q])>=15:
            print(' Quota vs base',cmp(hv,Q,B)[0]); print(' Quota vs MF',cmp(hv,Q,M)[0]); print(' Quota vs base held',cmp(hv,Q,B,lambda i:i not in DEV)[0]); print(' Sipht_30 quota vs base',round(cmp(hv,Q,B)[1]['Sipht_30'],2))

if __name__=='__main__':
    report()
