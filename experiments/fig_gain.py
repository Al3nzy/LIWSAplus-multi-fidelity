import sys, json, collections
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from paths import *
import analyze as A
from extra_analysis import run, DEV, FAMS, fam
plt.rcParams.update({"font.family":"serif","font.size":7.5})
col={'Montage':'#1b9e77','CyberShake':'#d95f02','Sipht':'#7570b3','Epigenomics':'#e7298a','Inspiral':'#66a61e'}
mk={'Montage':'o','CyberShake':'s','Sipht':'^','Epigenomics':'D','Inspiral':'v'}
fig,axs=plt.subplots(1,2,figsize=(3.5,2.9),sharey=True)
for ax,(st,lab) in zip(axs,(('rc','Risk-aware (CVaR)'),('rm','Expected value'))):
    rel=run(st)['paper (1.2 x worst, all 12+ conditions)']
    ys=[]
    for yi,f in enumerate(FAMS):
        items=sorted([(i,x) for i,x in rel.items() if fam(i)==f],key=lambda t:int(t[0].split('_')[1]))
        for j,(i,x) in enumerate(items):
            held = i not in DEV
            ax.scatter(x, yi+(j-1.5)*0.16, s=16, marker=mk[f], facecolor=col[f] if held else 'white', edgecolor=col[f], linewidth=0.9, zorder=3)
            if x<-5: ax.annotate(i, (x, yi+(j-1.5)*0.16), xytext=(4,-1), textcoords='offset points', fontsize=6)
    ax.axvline(0,color='k',lw=0.6)
    ax.set_yticks(range(5)); ax.set_yticklabels(FAMS); ax.tick_params(labelleft=(ax is axs[0]))
    ax.set_title(lab,fontsize=7.5); ax.set_xlabel('gain over LIWSA+ (%)'); ax.grid(axis='x',lw=0.3)
    ax.set_xlim(-48,18) if st=='rc' else ax.set_xlim(-2,16)
axs[0].invert_yaxis()
fig.text(0.5,0.005,'filled: held-out (15); open: development (5)',ha='center',fontsize=6.5)
fig.tight_layout(rect=(0,0.03,1,1))
(ANALYSIS/'figures').mkdir(exist_ok=True)
fig.savefig(ANALYSIS/'figures'/'fig_gain.pdf'); fig.savefig(ANALYSIS/'figures'/'fig_gain.png',dpi=160)
