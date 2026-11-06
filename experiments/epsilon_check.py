"""Reference-point-free check of the MF-versus-LIWSA+ comparison: multiplicative epsilon indicator of each run's out-of-sample
front against the per-workflow reference front (non-dominated union of all out-of-sample points of all conditions in the file).
Smaller is better; 1.0 means the run attains the reference front. Reported as a robustness check of the hypervolume result."""
import json, collections
import numpy as np
from scipy import stats
from paths import *
import analyze as A
from extra_analysis import DEV, load

def nd(pts):
    return A.nondom([tuple(p) for p in pts])

def eps(front, ref):
    f = np.array(front); r = np.array(ref)
    # for each reference point the best (smallest) worst-ratio over the run's points
    return float(np.max(np.min(np.max(f[:, None, :] / r[None, :, :], axis=2), axis=0)))

def compute(setting, a='LIWSAP.f.r.san-plus', b='LIWSAP.f.r-plus'):
    runs = load(setting)
    allp = collections.defaultdict(list)
    for r in runs:
        allp[r['instance']] += r.get('oos', [])
    ref = {i: nd(v) for i, v in allp.items()}
    e = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in runs:
        if r['cond'] in (a, b):
            e[r['cond']][r['instance']].append(eps(r['oos'], ref[r['instance']]))
    inst = sorted(set(e[a]) & set(e[b]))
    rel = {i: 100 * (np.mean(e[a][i]) - np.mean(e[b][i])) / np.mean(e[b][i]) for i in inst}   # negative = MF closer to the reference front
    return inst, rel, {i: (np.mean(e[a][i]), np.mean(e[b][i])) for i in inst}

def main():
    L = []
    for st, lab in (('rc', 'Risk-aware'), ('rm', 'Expected value')):
        inst, rel, raw = compute(st)
        for name, sub in (('all 20', inst), ('15 held-out', [i for i in inst if i not in DEV])):
            v = np.array([rel[i] for i in sub]); x = np.array([raw[i][0] for i in sub]); y = np.array([raw[i][1] for i in sub])
            p = stats.wilcoxon(x, y).pvalue
            print(f'{lab:15s} {name:12s} eps change {v.mean():+.2f}% (median {np.median(v):+.2f}%), MF closer on {(v<0).sum()}, further on {(v>0).sum()}, p={p:.3f}; mean eps MF {x.mean():.4f}, LIWSA+ {y.mean():.4f}')
            L.append(f'{lab if name=="all 20" else ""} & {name} & {x.mean():.3f} & {y.mean():.3f} & {v.mean():+.2f} & {(v<0).sum()}/{(v>0).sum()} & ' + ('$<$0.001' if p < 0.001 else f'{p:.3f}'))
    open(str(ANALYSIS) + '/epsilon.tex', 'w').write(' \\\\\n'.join(L))

if __name__ == '__main__':
    main()
