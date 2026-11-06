"""Audit of the multiple-comparison correction: recomputes every Holm-adjusted p value of the scenario-based tables from the raw
per-workflow data with an independent implementation (statsmodels-free step-down Holm) and compares it with the values the tables
were built from. Families: one per objective and per subset table (27 comparisons each)."""
import re, numpy as np
from scipy import stats
from paths import *
import analyze as A
import paperlib as P

def holm_independent(p):
    m = len(p); idx = sorted(range(m), key=lambda i: p[i]); adj = [0.0] * m; prev = 0.0
    for rank, i in enumerate(idx):
        prev = max(prev, (m - rank) * p[i]); adj[i] = min(1.0, prev)
    return adj

def main():
    src = open(EXP / 'paper_statistics.py').read()
    block = src[src.index('SAP = ['):src.index('SALAB')]
    pairs = re.findall(r"\('(LIWSAP[^']*)', '([^']*)'\)", block)
    DEV = {'Montage_50', 'CyberShake_50', 'Sipht_60', 'Epigenomics_46', 'Inspiral_50'}
    out = [f'family size per objective and subset table: {len(pairs)} comparisons']
    worst = 0.0
    for st in ('rc', 'rm'):
        rows = P.sa_rows(st)
        if not rows:
            continue
        for name, sub in (('all 20', None), ('15 held-out', lambda i: i not in DEV), ('five 1000-task', lambda i: i in A.LARGE)):
            raw = []; used = []
            for a, b in pairs:
                c = A.compare(rows, a, b, 'hv_arc', sub)
                raw.append(c['p'] if c else 1.0)
            mine = holm_independent(raw)
            ref = A.holm(np.array(raw))
            d = float(np.max(np.abs(np.array(mine) - ref))); worst = max(worst, d)
            out.append(f'{st} {name:15s}: max |independent Holm - table Holm| = {d:.2e}; smallest raw p {min(raw):.4f}, smallest adjusted {min(mine):.4f}')
    out.append(f'overall maximum difference {worst:.2e}')
    text = '\n'.join(out); print(text)
    open(str(ANALYSIS) + '/report_holm_audit.txt', 'w').write(text + '\n')

if __name__ == '__main__':
    main()
