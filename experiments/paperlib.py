"""Shared loading and metric code for paper_statistics.py and make_paper_figures.py (faithful-decoder study)."""
import json, glob, csv, collections
import numpy as np
from paths import *
import analyze as A

NAMES = {'NSGA2.f-plus': 'NSGA-II', 'SPEA2.f-plus': 'SPEA2', 'MOEAD.f-plus': 'MOEA/D', 'LIWSA.f-plus': 'LIWSA (original moves)',
         'LIWSAP.f-plus': 'LIWSA+', 'LIWSAP.nog.f-plus': 'LIWSA+ without guided moves', 'LIWSAP.noe.f-plus': 'LIWSA+ without estimate move',
         'LIWSAP.nol.f-plus': 'LIWSA+ without locust moves', 'LIWSAP.noa.f-plus': 'LIWSA+ without online operator learning',
         'NSGA2.f-warm': 'NSGA-II (HEFT/Min-Min seeds only)', 'LIWSA.f-warm': 'LIWSA (HEFT/Min-Min seeds only)'}
SCEN = {'original': 'faith_', 'linear': 'faithlin_', 'premium': 'faithprem_'}
FAMS = ['Montage', 'CyberShake', 'Sipht', 'Epigenomics', 'Inspiral']
ORDER = ['Montage_25', 'Montage_50', 'Montage_100', 'Montage_1000', 'CyberShake_30', 'CyberShake_50', 'CyberShake_100',
         'CyberShake_1000', 'Sipht_30', 'Sipht_60', 'Sipht_100', 'Sipht_1000', 'Epigenomics_24', 'Epigenomics_46',
         'Epigenomics_100', 'Epigenomics_997', 'Inspiral_30', 'Inspiral_50', 'Inspiral_100', 'Inspiral_1000']
DEV = {'Montage_50', 'CyberShake_50', 'Sipht_60', 'Epigenomics_46', 'Inspiral_50'}


def files(prefix):
    return sorted(glob.glob(str(RESULTS / (prefix + '*.jsonl'))))


def scenario(prefix, extra_files=()):
    """-> rows, heuristics, reference points for a pricing scenario (faithful objective)."""
    refs, runs = A.load(files(prefix) + list(extra_files))
    for r in refs.values():
        r['heft'], r['minmin'] = r['heft_f'], r['minmin_f']
    return A.build(refs, runs)


def raw_runs(prefix):
    out = []
    for p in files(prefix):
        for l in open(p):
            r = json.loads(l)
            if r['type'] == 'run':
                out.append(r)
    return out


def holm_table(rows, pairs, key='hv_arc', subset=None):
    res = []
    for a, b in pairs:
        c = A.compare(rows, a, b, key, subset)
        res.append((a, b, c))
    ps = np.array([r[2]['p'] if r[2] else 1.0 for r in res])
    hp = A.holm(ps) if len(ps) else ps
    return [(a, b, c, h) for (a, b, c), h in zip(res, hp)]


def pfmt(p):
    return '<0.001' if p < 0.001 else f'{p:.3f}'


def sens_data():
    """-> {config: (rows)} with the default configuration (seeds 1-5 of the main run) as 'default'."""
    base_refs, base_runs = A.load(files('faith_'))
    default = [dict(r, cond='default') for r in base_runs if r['cond'] == 'LIWSAP.f-plus' and r['seed'] <= 5]
    out = {}
    for f in sorted(glob.glob(str(RESULTS / 'sens_*.jsonl'))):
        name = Path(f).stem[5:]
        refs, runs = A.load([f])
        out[name] = [dict(r, cond=name) for r in runs]
    return base_refs, default, out


def sens_table():
    base_refs, default, cfgs = sens_data()
    for r in base_refs.values():
        r['heft'], r['minmin'] = r['heft_f'], r['minmin_f']
    runs = list(default)
    for v in cfgs.values():
        runs += v
    rows, heur, ref = A.build(base_refs, runs)
    res = []
    for name in cfgs:
        res.append((name, A.compare(rows, name, 'default', 'hv_arc')))
    return res


def robust_data(path=None):
    rows = []
    for l in open(path or (RESULTS / 'robust.jsonl')):
        r = json.loads(l)
        if r['type'] == 'rob':
            rows.append(r)
    return rows


def robust_table(setting_list=('len10', 'len25', 'mips', 'both')):
    """-> {setting: {cond: {instance: mean HV over scenarios and seeds}}} plus nominal HV."""
    rows = robust_data()
    out = {}
    for setting in ('nominal',) + tuple(setting_list):
        ref = {}
        for r in rows:
            fronts = [r['nominal']] if setting == 'nominal' else r['settings'][setting]
            for fr in fronts:
                for m, c in fr:
                    o = ref.get(r['instance'], (0, 0))
                    ref[r['instance']] = (max(o[0], 1.2 * m), max(o[1], 1.2 * c))
        acc = collections.defaultdict(lambda: collections.defaultdict(list))
        for r in rows:
            fronts = [r['nominal']] if setting == 'nominal' else r['settings'][setting]
            hv = np.mean([A.hv2d(fr, ref[r['instance']]) for fr in fronts])
            acc[r['cond']][r['instance']].append(hv)
        out[setting] = {c: {i: float(np.mean(v)) for i, v in d.items()} for c, d in acc.items()}
    return out


SA_NAMES = {'NSGA2.f.r-plus': 'NSGA-II', 'SPEA2.f.r-plus': 'SPEA2', 'MOEAD.f.r-plus': 'MOEA/D', 'LIWSAP.f.r-plus': 'LIWSA+',
            'LIWSAP.f.r.sar-plus': 'LIWSA+ with random pre-selection', 'LIWSAP.f.r.san-plus': 'MF-LIWSA+ (nominal pre-screening)',
            'LIWSAP.f.r.sao-plus': 'MF-LIWSA+ with online ridge correction', 'LIWSAP.f.r.sao.hvi-plus': 'MF-LIWSA+ with ridge and HVI rule',
            'LIWSAP.f.r.imm-plus': 'LIWSA+ with OLS immigrants',
            'LIWSAP.f.r.sal-plus': 'MF-LIWSA+ with rich-feature ridge', 'LIWSAP.f.r.sag-plus': 'MF-LIWSA+ with boosted trees',
            'LIWSAP.f.r.sal.hvi-plus': 'MF-LIWSA+ with rich ridge and HVI rule', 'LIWSAP.f.r.sag.hvi-plus': 'MF-LIWSA+ with boosted trees and HVI rule'}
SA_SET = {'rc': 'risk-aware (CVaR, stragglers)', 'rm': 'expected value'}


def sa_rows(setting):
    import analyze_surrogate as S
    fs = [str(RESULTS / f'rob_{setting}_small.jsonl'), str(RESULTS / f'rob_{setting}_large.jsonl')]
    fs = [f for f in fs if Path(f).exists()]
    if not fs:
        return None
    refs, runs = S.load(fs)
    return S.build(runs)
