#!/usr/bin/env python3
"""Checks on the learned surrogates (run after `multifidelity`): ranking quality, candidate-pool size, harsher tail risk and hedged selection.
All checks use the five development workflows (and three extra ones for the ranking diagnostics) before any claim about the full study.
    python experiments/run_experiments.py mlchecks [--quick]
Writes experiments/work/analysis/ml_checks.txt."""
import json, sys, collections
import numpy as np
from paths import *
import analyze_surrogate as S

DEV = ['Montage_50', 'CyberShake_50', 'Sipht_60', 'Epigenomics_46', 'Inspiral_50']
DIAG_SET = DEV + ['Sipht_30', 'Montage_25', 'Epigenomics_24']
RC = ['-Dgens=20', '-Drobust=30', '-Dcvar=0.8', '-Dslowp=0.2', '-Dmipslow=0.8']
RM = ['-Dgens=20', '-Drobust=20']
HEAVY = ['-Dgens=20', '-Drobust=30', '-Dcvar=0.9', '-Dslowp=0.4', '-Dmipslow=0.5']
PRE = 'LIWSAP.f.r'


def run_all(runner, quick):
    dev = ['Montage_25'] if quick else DEV
    diag = ['Montage_25'] if quick else DIAG_SET
    ns, nd = (1, 1) if quick else (5, 2)
    out = RESULTS
    for c in ('sao', 'sal', 'sag'):
        runner(out / f'ml_diag_{c}.jsonl', ','.join(diag), nd, [f'{PRE}.{c}-plus'], RC + ['-Ddiag=true'])
        runner(out / f'ml_heavydiag_{c}.jsonl', ','.join(dev), nd, [f'{PRE}.{c}-plus'], HEAVY + ['-Ddiag=true'])
    for K in (10, 20):
        for st, jvm in (('rc', RC), ('rm', RM)):
            runner(out / f'ml_k{K}_{st}.jsonl', ','.join(dev), ns, [f'{PRE}.{c}-plus' for c in ('san', 'sal', 'sag')], jvm + [f'-Dsak={K}'])
    runner(out / 'ml_heavy.jsonl', ','.join(dev), ns, ['NSGA2.f.r-plus', f'{PRE}-plus'] + [f'{PRE}.{c}-plus' for c in ('san', 'sao', 'sal', 'sag')], HEAVY)
    for st, jvm in (('rc', RC), ('rm', RM)):
        runner(out / f'ml_mix_{st}.jsonl', ','.join(dev), ns, [f'{PRE}.{c}.mix-plus' for c in ('sao', 'sal', 'sag')], jvm)


def _recs(path):
    return [r for r in map(json.loads, open(path)) if r['type'] == 'run'] if Path(path).exists() else []


def analyse(quick):
    L = []
    say = lambda *a: L.append(' '.join(str(x) for x in a))
    dev = ['Montage_25'] if quick else DEV
    say('== RANKING QUALITY (mean over workflows; model vs nominal), risk-aware setting, 2 seeds')
    for c, lab in (('sao', 'ridge on 10 descriptors'), ('sal', 'ridge on rich features'), ('sag', 'boosted trees on rich features')):
        d = [r['diag'] for r in _recs(RESULTS / f'ml_diag_{c}.jsonl')]
        if d:
            m = np.mean(d, 0)
            say(f'{lab}: rho makespan nominal {m[0]:.3f} -> model {m[1]:.3f} | rho cost {m[2]:.3f} -> {m[3]:.3f} | overlap {m[4]:.3f} -> {m[5]:.3f} (chance 0.20)')
    base = {st: [r for r in _recs(RESULTS / f'rob_{st}_small.jsonl') if r['instance'] in dev and r['seed'] <= (1 if quick else 5)] for st in ('rc', 'rm')}
    for st in ('rc', 'rm'):
        say(f'\n== CANDIDATE-POOL SIZE K ({st}), development workflows: mean gain over LIWSA+ (%) and mean search time (ms)')
        for K in (5, 10, 20):
            runs = [dict(r) for r in base[st] if r['cond'] in (f'{PRE}-plus', 'NSGA2.f.r-plus') or K == 5]
            if K != 5:
                runs = [r for r in base[st] if r['cond'] in (f'{PRE}-plus', 'NSGA2.f.r-plus')] + _recs(RESULTS / f'ml_k{K}_{st}.jsonl')
            else:
                runs = [r for r in base[st] if r['cond'] in (f'{PRE}-plus', f'{PRE}.san-plus', f'{PRE}.sal-plus', f'{PRE}.sag-plus')]
            if not runs:
                continue
            rows = S.build(runs)
            def m(c):
                d = collections.defaultdict(list)
                for r in rows:
                    if r['cond'] == c: d[r['instance']].append(r['hv_arc'])
                return {i: np.mean(v) for i, v in d.items()}
            Lb = m(f'{PRE}-plus'); parts = []
            for c, lab in (('san', 'nominal'), ('sal', 'rich ridge'), ('sag', 'boosted')):
                mm = m(f'{PRE}.{c}-plus')
                if mm and Lb:
                    g = np.mean([100 * (mm[i] / Lb[i] - 1) for i in Lb if i in mm]); t = np.mean([r['search_ms'] for r in rows if r['cond'] == f'{PRE}.{c}-plus'])
                    parts.append(f'{lab} {g:+.2f}% ({t:.0f} ms)')
            say(f'K={K}: ' + ' | '.join(parts))
    say('\n== HARSHER TAIL RISK (CVaR 0.9, 40% stragglers), development workflows: gain over LIWSA+ | vs nominal pre-screening')
    runs = _recs(RESULTS / 'ml_heavy.jsonl')
    if runs:
        rows = S.build(runs)
        def m(c):
            d = collections.defaultdict(list)
            for r in rows:
                if r['cond'] == c: d[r['instance']].append(r['hv_arc'])
            return {i: np.mean(v) for i, v in d.items()}
        Lb, N = m(f'{PRE}-plus'), m(f'{PRE}.san-plus')
        for c, lab in (('NSGA2.f.r-plus', 'NSGA-II'), (f'{PRE}.san-plus', 'nominal'), (f'{PRE}.sao-plus', 'ridge10'), (f'{PRE}.sal-plus', 'rich ridge'), (f'{PRE}.sag-plus', 'boosted')):
            mm = m(c)
            say(f'  {lab:11s} {np.mean([100 * (mm[i] / Lb[i] - 1) for i in Lb]):+6.2f}% | {np.mean([100 * (mm[i] / N[i] - 1) for i in Lb]):+6.2f}%')
        for c in ('sao', 'sal', 'sag'):
            d = [r['diag'] for r in _recs(RESULTS / f'ml_heavydiag_{c}.jsonl')]
            if d:
                mm = np.mean(d, 0); say(f'  ranking {c}: rho makespan {mm[0]:.3f} -> {mm[1]:.3f}, overlap {mm[4]:.3f} -> {mm[5]:.3f}')
    for st in ('rc', 'rm'):
        say(f'\n== HEDGED SELECTION ({st}), development workflows: gain over LIWSA+ | vs nominal pre-screening | mean search time (ms)')
        runs = [r for r in base[st] if r['cond'] in (f'{PRE}-plus', f'{PRE}.san-plus', f'{PRE}.sao-plus', f'{PRE}.sal-plus', f'{PRE}.sag-plus')] + _recs(RESULTS / f'ml_mix_{st}.jsonl')
        if not runs:
            continue
        rows = S.build(runs)
        def m(c):
            d = collections.defaultdict(list)
            for r in rows:
                if r['cond'] == c: d[r['instance']].append(r['hv_arc'])
            return {i: np.mean(v) for i, v in d.items()}
        Lb, N = m(f'{PRE}-plus'), m(f'{PRE}.san-plus')
        for c, lab in (('san', 'nominal'), ('sao', 'ridge10'), ('sao.mix', 'ridge10 + hedge'), ('sal', 'rich ridge'), ('sal.mix', 'rich ridge + hedge'), ('sag', 'boosted'), ('sag.mix', 'boosted + hedge')):
            mm = m(f'{PRE}.{c}-plus')
            if mm:
                say(f'  {lab:20s} {np.mean([100 * (mm[i] / Lb[i] - 1) for i in Lb]):+6.2f}% | {np.mean([100 * (mm[i] / N[i] - 1) for i in Lb]):+6.2f}% | {np.mean([r["search_ms"] for r in rows if r["cond"] == f"{PRE}.{c}-plus"]):.0f}')
    open(ANALYSIS / 'ml_checks.txt', 'w').write('\n'.join(L) + '\n')
    print('\n'.join(L))
