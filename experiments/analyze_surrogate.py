"""Analysis of the expensive-objective (Monte-Carlo robust) study: hypervolume of the OUT-OF-SAMPLE expected-objective front
(archive genotypes re-evaluated under 100 fresh scenarios), paired statistics, time, and evaluation counts.
Usage: python experiments/analyze_surrogate.py [file ...]  (default: work/results/rob_*.jsonl)"""
import sys, json, glob, collections
import numpy as np
from paths import *
import analyze as A


def load(paths):
    refs, runs = {}, []
    for p in paths:
        for l in open(p):
            r = json.loads(l)
            (refs.__setitem__(r['instance'], r) if r['type'] == 'ref' else runs.append(r))
    return refs, runs


def build(runs):
    pts = collections.defaultdict(list)
    for r in runs:
        pts[r['instance']] += r.get('oos', [])
    ref = {i: (1.2 * max(p[0] for p in v), 1.2 * max(p[1] for p in v)) for i, v in pts.items()}
    rows = []
    for r in runs:
        rows.append(dict(instance=r['instance'], cond=r['cond'], seed=r['seed'], evals=r['evals'], search_ms=r['search_ms'],
                         hv_arc=A.hv2d(r.get('oos', []), ref[r['instance']]), nom_dec=r.get('nom_dec', 0), hv_pop=0.0,
                         n_unique_front=len(r['front'])))
    return rows


def main(paths, pairs=None, subset=None, title=''):
    refs, runs = load(paths)
    rows = build(runs)
    conds = sorted(set(r['cond'] for r in rows))
    print(f'runs {len(rows)} | conditions {conds}')
    for c in conds:
        t = [r['search_ms'] for r in rows if r['cond'] == c]
        e = [r['evals'] for r in rows if r['cond'] == c]
        nd = [r['nom_dec'] for r in rows if r['cond'] == c]
        print(f'{c:>26} search {np.mean(t):9.1f} ms  expensive evals {np.mean(e):6.0f}  nominal decodes {np.mean(nd):8.0f}')
    return rows, conds


if __name__ == '__main__':
    paths = sys.argv[1:] or sorted(glob.glob(str(RESULTS / 'rob_*.jsonl')))
    rows, conds = main(paths)
    base = 'LIWSAP.f.r-plus'
    DEV = {'Montage_50', 'CyberShake_50', 'Sipht_60', 'Epigenomics_46', 'Inspiral_50'}
    pairs = [(c, base) for c in conds if c != base]
    for name, sub in (('ALL', None),):
        print(f'\n== out-of-sample hypervolume gain vs {base} ({name}) ==')
        for a, b in pairs:
            print(f'{a:>26} vs {b}: {A.fmt(A.compare(rows, a, b, "hv_arc", sub))}')
