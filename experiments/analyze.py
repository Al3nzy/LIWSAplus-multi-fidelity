"""Shared analysis functions (loading, hypervolume, paired tests, Holm correction, bootstrap intervals).

Conventions (same as the paper): hypervolume of the (makespan, cost) front, both minimised, against ONE shared
reference point per instance = 1.2 x the worst value seen on any front of any condition and seed for that
instance. Two fronts are scored per run:
  HV_pop : non-dominated set of the FINAL POPULATION (the convention used in the submitted paper)
  HV_arc : non-dominated set of EVERY solution evaluated during the run (anytime archive)
Statistics: per-instance means over seeds, then paired Wilcoxon signed-rank across instances (exact for n<=25),
Holm correction within each family of comparisons, bootstrap 95% CI of the mean relative difference.
"""
import json, glob, sys, os, collections
from paths import *
import numpy as np
from scipy import stats

FAM = ['Montage', 'CyberShake', 'Sipht', 'Epigenomics', 'Inspiral']
LARGE = {'Montage_1000', 'CyberShake_1000', 'Sipht_1000', 'Epigenomics_997', 'Inspiral_1000'}


def hv2d(pts, ref):
    pts = sorted((p for p in pts if p[0] < ref[0] and p[1] < ref[1]))
    hv, best_c = 0.0, ref[1]
    for m, c in pts:
        if c < best_c:
            hv += (ref[0] - m) * (best_c - c)
            best_c = c
    return hv


def nondom(pts):
    pts = sorted(set(map(tuple, pts)))
    out, best_c = [], float('inf')
    for m, c in pts:
        if c < best_c:
            out.append((m, c)); best_c = c
    return out


def load(paths):
    refs, runs = {}, []
    for p in paths:
        for line in open(p):
            r = json.loads(line)
            if r['type'] == 'ref':
                refs[r['instance']] = r
            else:
                runs.append(r)
    return refs, runs


def build(refs, runs):
    allpts = collections.defaultdict(list)
    for inst, r in refs.items():
        allpts[inst] += [r['heft'], r['minmin']]
    for r in runs:
        allpts[r['instance']] += r['front'] + r['archive']
    ref = {i: (1.2 * max(p[0] for p in v), 1.2 * max(p[1] for p in v)) for i, v in allpts.items()}
    rows = []
    for r in runs:
        rp = ref[r['instance']]
        d = dict(instance=r['instance'], cond=r['cond'], seed=r['seed'], evals=r['evals'],
                 search_ms=r['search_ms'], seeding_ms=r.get('seeding_ms', 0.0) or 0.0,
                 hv_pop=hv2d(r['front'], rp),
                 hv_arc=hv2d(r['archive'], rp) if r['archive'] else float('nan'),
                 n_unique_front=len(nondom(r['front'])),
                 rho_m=r.get('rho_m'), rho_c=r.get('rho_c'), trace=r['trace'], ref=rp)
        rows.append(d)
    heur = {}
    for inst, r in refs.items():
        rp = ref[inst]
        heur[inst] = dict(heft=hv2d([r['heft']], rp), minmin=hv2d([r['minmin']], rp))
    return rows, heur, ref


def inst_means(rows, cond, key):
    acc = collections.defaultdict(list)
    for r in rows:
        if r['cond'] == cond and not np.isnan(r[key]):
            acc[r['instance']].append(r[key])
    return {k: float(np.mean(v)) for k, v in acc.items()}


def compare(rows, a, b, key='hv_pop', subset=None):
    """Relative difference of condition a vs b (percent of b), paired over instances."""
    ma, mb = inst_means(rows, a, key), inst_means(rows, b, key)
    inst = sorted(set(ma) & set(mb))
    if subset is not None:
        inst = [i for i in inst if subset(i)]
    if len(inst) < 3:
        return None
    x = np.array([ma[i] for i in inst]); y = np.array([mb[i] for i in inst])
    rel = 100.0 * (x - y) / y
    try:
        p = stats.wilcoxon(x, y).pvalue if np.any(x != y) else 1.0
    except ValueError:
        p = 1.0
    rng = np.random.RandomState(0)
    boots = [rel[rng.randint(len(rel), size=len(rel))].mean() for _ in range(4000)]
    return dict(n=len(inst), mean=rel.mean(), median=float(np.median(rel)), wins=int((rel > 0).sum()),
                losses=int((rel < 0).sum()), p=p, lo=float(np.percentile(boots, 2.5)), hi=float(np.percentile(boots, 97.5)))


def holm(ps):
    order = np.argsort(ps); m = len(ps); out = np.empty(m); run = 0.0
    for k, i in enumerate(order):
        run = max(run, (m - k) * ps[i]); out[i] = min(1.0, run)
    return out


def fmt(c):
    if c is None:
        return 'n/a'
    return (f"{c['mean']:+6.2f}% [{c['lo']:+.2f},{c['hi']:+.2f}] med {c['median']:+.2f}%  "
            f"W/L {c['wins']}/{c['losses']} of {c['n']}  p={c['p']:.4f}")


def table(rows, pairs, key, subset=None, title=''):
    print(f'\n== {title}  [{key}] ==')
    res = [(a, b, compare(rows, a, b, key, subset)) for a, b in pairs]
    ps = np.array([r[2]['p'] if r[2] else 1.0 for r in res])
    hp = holm(ps) if len(ps) else ps
    for (a, b, c), h in zip(res, hp):
        print(f'{a:>24} vs {b:<20} {fmt(c)}  holm={h:.3f}')
