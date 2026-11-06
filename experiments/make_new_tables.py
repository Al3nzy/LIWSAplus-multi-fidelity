import sys,os,json
import numpy as np
from extra_analysis import run, DEV, FAMS, fam
import ideal_analysis as I
import rep_analysis as R
from paths import ANALYSIS
OUT=str(ANALYSIS)+'/'
def main():
    def row(v):
        v=np.array(v); return f'{v.mean():+.2f} ({(v>0).sum()}/{(v<0).sum()})'
    # --- ideal ablation table (main)
    L=[]
    for st,lab in (('rc','Risk-aware'),('rm','Expected value')):
        hv=I.go(st)
        def g(a,b,sub=None):
            r=I.cmp(hv,I.N[a],I.N[b],sub); return r
        for nm,a,b in (('Faithful screening (MF)','faith','base'),('Idealised screening','ideal','base'),('Random pre-selection','rand','base'),('MF over idealised screening','faith','ideal')):
            ra=g(a,b); rh=g(a,b,lambda i:i not in DEV)
            L.append(f'{lab if nm.startswith("Faithful") else ""} & {nm} & {ra["mean"]:+.2f} ({ra["W"]}/{ra["L"]}) & {ra["med"]:+.2f} & {rh["mean"]:+.2f} ({rh["W"]}/{rh["L"]}) & {ra["p"]:.3f} \\\\')
        L.append('\\addlinespace')
    open(OUT+'ideal.tex','w').write('\n'.join(L[:-1])+'\n')
    # --- ideal family breakdown (supp)
    F=[]
    for st,lab in (('rc','Risk-aware'),('rm','Expected value')):
        hv=I.go(st); rel=I.cmp(hv,I.N['faith'],I.N['ideal'])['rel']
        F.append(f'{lab} & '+' & '.join(f'{np.mean([x for i,x in rel.items() if fam(i)==f]):+.1f}' for f in FAMS)+' \\\\')
    open(OUT+'ideal_family.tex','w').write('\n'.join(F)+'\n')
    # --- reference point rule + family clusters (supp)
    RR=[];FC=[]
    for st,lab in (('rc','Risk-aware'),('rm','Expected value')):
        res=run(st)
        for k,(name,rel) in enumerate(res.items()):
            v=np.array(list(rel.values())); h=np.array([x for i,x in rel.items() if i not in DEV])
            nm={'paper (1.2 x worst, all 12+ conditions)':'1.2$\\times$ worst, all conditions (used)','1.2 x worst, compared pair only':'1.2$\\times$ worst, compared pair only','1.05 x worst, all conditions':'1.05$\\times$ worst, all conditions','2.0 x worst, all conditions':'2.0$\\times$ worst, all conditions'}[name]
            RR.append(f'{lab if k==0 else ""} & {nm} & {v.mean():+.2f} ({(v>0).sum()}/{(v<0).sum()}) & {np.median(v):+.2f} & {h.mean():+.2f} ({(h>0).sum()}/{(h<0).sum()}) \\\\')
        rel=res['paper (1.2 x worst, all 12+ conditions)']
        FC.append(f'{lab} & '+' & '.join(f'{np.mean([x for i,x in rel.items() if fam(i)==f]):+.2f}' for f in FAMS)+' \\\\')
    open(OUT+'refpoint.tex','w').write('\n'.join(RR)+'\n'); open(OUT+'family_gain.tex','w').write('\n'.join(FC)+'\n')
    # --- replication table
    RP=[]
    for st,lab in (('rc','Risk-aware'),('rm','Expected value')):
        rel0=run(st)['paper (1.2 x worst, all 12+ conditions)']
        small=[x for i,x in rel0.items() if int(i.split('_')[1])<=100]
        cells=[row(small)]
        for salt in (1,2):
            hv=R.hvs(f'rep{salt}_{st}_small.jsonl'); r=R.cmp(hv,R.M,R.B)[1]; cells.append(row(list(r.values())))
        RP.append(f'{lab} & '+' & '.join(cells)+' \\\\')
    open(OUT+'replication.tex','w').write('\n'.join(RP)+'\n')
    
    
    # --- quota table
    from scipy import stats
    L=[]
    for st,lab in (('rc','Risk-aware'),('rm','Expected value')):
        hv=R.hvs(f'quota_{st}_small.jsonl')
        for nm,a,b in (('MF over LIWSA+',R.M,R.B),('Quota over LIWSA+',R.Q,R.B),('Quota over MF',R.Q,R.M)):
            s,rel=R.cmp(hv,a,b); v=np.array(list(rel.values())); h=np.array([x for i,x in rel.items() if i not in R.DEV])
            x=np.array([np.mean(hv[a][i]) for i in rel]);y=np.array([np.mean(hv[b][i]) for i in rel])
            p=stats.wilcoxon(x,y).pvalue
            L.append(f'{lab if nm.startswith("MF") else ""} & {nm} & {v.mean():+.2f} ({(v>0).sum()}/{(v<0).sum()}) & {np.median(v):+.2f} & {h.mean():+.2f} & {rel["Sipht_30"]:+.2f} & '+('$<$0.001' if p<0.001 else f'{p:.3f}')+' \\\\')
        L.append('\\addlinespace')
    open(OUT+'quota_table.tex','w').write(r"""\begin{table}[h]
    \caption{Quota variant (a fifth of the kept children drawn at random from the screened pool): hypervolume gain (\%), fifteen smaller workflows, five seeds; unadjusted $p$}
    \label{tab:quota}
    \centering\scriptsize
    \setlength{\tabcolsep}{2.5pt}
    \begin{tabular}{@{}llccccc@{}}
    \toprule
    Objective & Comparison & Mean (W/L) & Median & Held-out 10 & Sipht\_30 & $p$ \\
    \midrule
    """+'\n'.join(L[:-1])+r"""
    \bottomrule
    \end{tabular}
    \end{table}
    """)
    # --- workflow-subset table (Table IV of the paper)
    LARGE = lambda i: int(i.split('_')[1]) >= 900
    def srow(rel, sub):
        v = np.array([x for i, x in rel.items() if sub(i)]); return f'{v.mean():+.2f} ({(v>0).sum()}/{(v<0).sum()})'
    SB = []
    for st, lab in (('rc', 'Risk-aware'), ('rm', 'Expected value')):
        for k, (nm, cond) in enumerate((('MF', 'LIWSAP.f.r.san-plus'), ('MF + ridge', 'LIWSAP.f.r.sao-plus'))):
            rel = run(st, a=cond)['paper (1.2 x worst, all 12+ conditions)']
            SB.append(f'{lab if k == 0 else ""} & {nm} & ' + ' & '.join(srow(rel, f) for f in (lambda i: True, lambda i: i in DEV, lambda i: i not in DEV, LARGE)) + ' \\\\')
    open(OUT + 'subsets.tex', 'w').write('\n'.join(SB).rstrip(' \\\\'))

    import re
    for f in ('ideal','ideal_family','family_gain','refpoint','replication'):
        p=OUT+f+'.tex'
        s=open(p).read().replace(' 0.000 \\\\',' $<$0.001 \\\\').replace('MF over idealised screening','MF over idealised').rstrip()
        open(p,'w').write(re.sub(r'\\\\\s*$','',s).rstrip())


if __name__=='__main__':
    main()
