"""Recomputes every statistic quoted in the paper and writes LaTeX tables (analysis/tables/*.tex) and a text report.
Reads experiments/work/results/*.jsonl, *.csv. Usage: python experiments/paper_statistics.py"""
import sys, json, csv, collections
import numpy as np
import pandas as pd
from paths import *
import analyze as A
from paperlib import *

TAB = ANALYSIS / 'tables'
TAB.mkdir(parents=True, exist_ok=True)
OUT = open(ANALYSIS / 'paper_statistics.txt', 'w')


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    OUT.write(s + '\n')


def cell(c, h):
    if c is None:
        return '--'
    return f"{c['mean']:+.2f} [{c['lo']:+.2f}, {c['hi']:+.2f}] & {c['wins']}/{c['losses']} & {pfmt(h)}"


def write(name, text):
    # the last row carries no row terminator: the .tex file writes "\\input{...}\\\\" so that no empty row is started before \\bottomrule
    t = text.rstrip('\n')
    if t.endswith('\\\\'):
        t = t[:-2].rstrip()
    (TAB / name).write_text(t + '\n')



if __name__ == '__main__' and not globals().get('_CHUNKED'):
    # run each section independently so that a missing input (e.g. in --quick mode) skips that section instead of aborting
    import re as _re
    _src = open(__file__).read()
    _ns = {'_CHUNKED': True, '__file__': __file__, '__name__': '__main__'}
    for _k, _p in enumerate(_re.split(r'\n(?=# -{20,})', _src)):
        try:
            exec(compile(_p, f'{__file__}[section {_k}]', 'exec'), _ns)
        except Exception as _e:
            print(f'WARNING: section {_k} skipped ({type(_e).__name__}: {_e})')
    raise SystemExit(0)

# ------------------------------------------------------------------ 1. fidelity
if (RESULTS / 'fidelity.csv').exists():
    d = pd.read_csv(RESULTS / 'fidelity.csv')
    lines = []
    say('== FIDELITY (vs full WorkflowSim), n =', len(d))
    for k, lab in (('dec', 'Analytic gap-filling decoder'), ('faith', 'Simulator-faithful decoder')):
        em = 100 * (d[k + '_makespan'] - d.sim_makespan).abs() / d.sim_makespan
        ec = 100 * (d[k + '_cost'] - d.sim_cost).abs() / d.sim_cost
        say(f'{lab}: makespan mean {em.mean():.3f}% median {em.median():.3f}% max {em.max():.2f}%; cost mean {ec.mean():.3f}% max {ec.max():.2f}%')
        lines.append(f"{lab} & {em.mean():.2f} & {em.median():.2f} & {em.max():.1f} & {ec.mean():.2f} & {ec.max():.1f} \\\\")
    write('fidelity.tex', '\n'.join(lines) + '\n')
    fam = d.assign(inst=d.instance.str.split('_').str[0])
    rows = []
    for f, g in fam.groupby('inst'):
        e0 = 100 * (g.dec_makespan - g.sim_makespan).abs() / g.sim_makespan
        e1 = 100 * (g.faith_makespan - g.sim_makespan).abs() / g.sim_makespan
        rows.append(f'{f} & {e0.mean():.2f} & {e0.max():.1f} & {e1.mean():.3f} & {e1.max():.2f} \\\\')
    write('fidelity_family.tex', '\n'.join(rows) + '\n')

# ------------------------------------------------------------------ 2. main comparison, three pricing scenarios
PAIRS = [('LIWSAP.f-plus', 'NSGA2.f-plus'), ('LIWSAP.f-plus', 'SPEA2.f-plus'), ('LIWSAP.f-plus', 'MOEAD.f-plus'), ('LIWSAP.f-plus', 'LIWSA.f-plus')]
store = {}
for sname, prefix in SCEN.items():
    if not files(prefix):
        continue
    rows, heur, ref = scenario(prefix)
    store[sname] = (rows, heur, ref)
    say(f'\n== MAIN COMPARISON, {sname} pricing, runs={len(rows)}')
    for key in ('hv_arc', 'hv_pop'):
        for r in holm_table(rows, PAIRS + [('SPEA2.f-plus', 'NSGA2.f-plus')], key):
            say(f'{key} {NAMES[r[0]]:>10} vs {NAMES[r[1]]:<22} {A.fmt(r[2])} holm={pfmt(r[3])}')
    for h in ('heft', 'minmin'):
        m = A.inst_means(rows, 'LIWSAP.f-plus', 'hv_arc')
        rel = [100 * (m[i] - heur[i][h]) / heur[i][h] for i in m if heur[i][h] > 0]
        say(f'LIWSA+ vs {h}: mean {np.mean(rel):+.1f}% median {np.median(rel):+.1f}%')
    m = A.inst_means(rows, 'NSGA2.f-plus', 'hv_arc'); w = A.inst_means(rows, 'NSGA2.f-warm', 'hv_arc')
    if w:
        rel = [100 * (m[i] / w[i] - 1) for i in m if i in w]
        say(f'NSGA-II with sweep seeds vs HEFT/Min-Min seeds only: mean {np.mean(rel):+.2f}% median {np.median(rel):+.2f}%')
if store:
    # main table: rows = baselines, column groups = pricing schemes, archive front
    L = []
    for a, b in PAIRS:
        parts = []
        for sname in ('original', 'linear', 'premium'):
            if sname not in store:
                continue
            tab = {(x[0], x[1]): x for x in holm_table(store[sname][0], PAIRS, 'hv_arc')}
            x = tab[(a, b)]
            parts.append(f"{x[2]['mean']:+.2f} & {x[2]['wins']}/{x[2]['losses']} & {pfmt(x[3])}")
        L.append(f'{NAMES[b]} & ' + ' & '.join(parts) + ' \\\\')
    write('main_comparison.tex', '\n'.join(L) + '\n')
    L = []
    for a, b in PAIRS:
        parts = []
        for sname in ('original', 'linear', 'premium'):
            if sname not in store:
                continue
            tab = {(x[0], x[1]): x for x in holm_table(store[sname][0], PAIRS, 'hv_pop')}
            x = tab[(a, b)]
            parts.append(f"{x[2]['mean']:+.2f} & {x[2]['wins']}/{x[2]['losses']} & {pfmt(x[3])}")
        L.append(f'{NAMES[b]} & ' + ' & '.join(parts) + ' \\\\')
    write('main_comparison_pop.tex', '\n'.join(L) + '\n')

# held-out workflows only (design choices were fixed on the five development workflows)
if 'original' in store:
    sub = lambda i: i not in DEV
    say('\n== MAIN COMPARISON on the 15 workflows not used for any design decision (original pricing)')
    L = []
    for key, tag in (('hv_arc', 'archive'), ('hv_pop', 'population')):
        for a, b, c, h in holm_table(store['original'][0], PAIRS, key, sub):
            say(f'{tag} {NAMES[a]} vs {NAMES[b]}: {A.fmt(c)} holm={pfmt(h)}')
            L.append(f'{tag} & {NAMES[b]} & ' + cell(c, h) + ' \\\\')
    write('main_comparison_test.tex', '\n'.join(L) + '\n')

# ------------------------------------------------------------------ 3. ablation (original pricing)
if 'original' in store:
    rows = store['original'][0]
    abl = [('LIWSAP.f-plus', 'LIWSAP.nog.f-plus'), ('LIWSAP.f-plus', 'LIWSAP.noe.f-plus'), ('LIWSAP.f-plus', 'LIWSAP.nol.f-plus'),
           ('LIWSAP.f-plus', 'LIWSAP.noa.f-plus'), ('NSGA2.f-plus', 'NSGA2.f-warm'), ('LIWSA.f-plus', 'LIWSA.f-warm')]
    say('\n== ABLATION (original pricing), A vs B = effect of component/seeds')
    L = []
    for key in ('hv_arc',):
        for a, b, c, h in holm_table(rows, abl, key):
            say(f'{key} {NAMES[a]} vs {NAMES[b]}: {A.fmt(c)} holm={pfmt(h)}')
            L.append(f'{NAMES[b]} & ' + cell(c, h) + ' \\\\')
    write('ablation.tex', '\n'.join(L) + '\n')
    # per-family
    L = []
    a = A.inst_means(rows, 'LIWSAP.f-plus', 'hv_arc'); b = A.inst_means(rows, 'NSGA2.f-plus', 'hv_arc')
    say('\n== PER FAMILY, LIWSA+ vs NSGA-II (archive front, original pricing)')
    for f in FAMS:
        v = [100 * (a[i] / b[i] - 1) for i in ORDER if i.startswith(f) and i in a and i in b]
        say(f'{f}: mean {np.mean(v):+.2f}% per workflow {[round(x, 2) for x in v]}')
        L.append(f"{f} & " + ' & '.join(f'{x:+.2f}' for x in v) + ' \\\\')
    write('family.tex', '\n'.join(L) + '\n')
    # per-workflow table of mean HV (archive front) for the main algorithms
    cols = ['NSGA2.f-plus', 'SPEA2.f-plus', 'MOEAD.f-plus', 'LIWSA.f-plus', 'LIWSAP.f-plus']
    means = {c: A.inst_means(rows, c, 'hv_arc') for c in cols}
    heur = store['original'][1]
    L = []
    for i in ORDER:
        if i not in means['LIWSAP.f-plus']:
            continue
        best = max(means[c][i] for c in cols)
        cells = []
        for c in cols:
            v = means[c][i] / best * 100
            cells.append(f'\\textbf{{{v:.2f}}}' if abs(v - 100) < 1e-9 else f'{v:.2f}')
        L.append(i.replace('_', '\\_') + ' & ' + ' & '.join(cells) + f' & {heur[i]["heft"] / best * 100:.1f} & {heur[i]["minmin"] / best * 100:.1f} \\\\')
    write('per_workflow.tex', '\n'.join(L) + '\n')
    # time
    L = []
    say('\n== SEARCH TIME (ms), small / 1000-task')
    for c in cols + ['LIWSAP.nog.f-plus']:
        s = [r['search_ms'] for r in rows if r['cond'] == c and r['instance'] not in A.LARGE]
        g = [r['search_ms'] for r in rows if r['cond'] == c and r['instance'] in A.LARGE]
        say(f'{NAMES[c]}: {np.mean(s):.0f} / {np.mean(g):.0f}')
        if c in cols:
            L.append(f'{NAMES[c]} & {np.mean(s):.0f} & {np.mean(g):.0f} \\\\')
    write('time.tex', '\n'.join(L) + '\n')
    u = {c: np.mean([r['n_unique_front'] for r in rows if r['cond'] == c]) for c in cols}
    say('distinct non-dominated points in the final population:', {NAMES[c]: round(v, 1) for c, v in u.items()})

# ------------------------------------------------------------------ 3b. every metric for every workflow (original prices)
if 'original' in store:
    rows = store['original'][0]
    heur = store['original'][1]
    runs = raw_runs('faith_')
    algs = ['HEFT', 'MINMIN', 'NSGA2.f-plus', 'SPEA2.f-plus', 'MOEAD.f-plus', 'LIWSA.f-plus', 'LIWSAP.f-plus']
    refs_, _ = A.load(files('faith_'))
    acc = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in runs:
        a = r['archive']
        if not a:
            continue
        acc[(r['cond'], r['instance'])]['mk'].append(min(p[0] for p in a))
        acc[(r['cond'], r['instance'])]['co'].append(min(p[1] for p in a))
    def metric(alg, inst, what):
        if alg in ('HEFT', 'MINMIN'):
            r = refs_[inst]['heft_f' if alg == 'HEFT' else 'minmin_f']
            return {'mk': r[0], 'co': r[1]}.get(what)
        v = acc.get((alg, inst), {}).get(what)
        return float(np.mean(v)) if v else None
    def meanrow(alg, inst, key):
        v = [r[key] for r in rows if r['cond'] == alg and r['instance'] == inst]
        return float(np.mean(v)) if v else None
    def num(v):
        if v is None:
            return '--'
        if v >= 1000:
            return f'{v:,.0f}'.replace(',', '\\,')
        return f'{v:.1f}' if v >= 100 else f'{v:.2f}'
    def table_rows(getter, best, fmt=num, only=None):
        L = []
        for inst in ORDER:
            vals = [getter(a, inst) for a in algs]
            cand = [v for a, v in zip(algs, vals) if v is not None and (only is None or a in only)]
            if not cand:
                continue
            b = best(cand)
            cells = []
            for a, v in zip(algs, vals):
                t = fmt(v)
                if v is not None and (only is None or a in only) and t == fmt(b):
                    t = '\\textbf{' + t + '}'
                cells.append(t)
            L.append(inst.replace('_', '\\_') + ' & ' + ' & '.join(cells) + ' \\\\')
        return L
    say('\n== PER-WORKFLOW METRIC TABLES written (makespan, cost, hypervolume, distinct points, time)')
    write('metric_makespan.tex', '\n'.join(table_rows(lambda a, i: metric(a, i, 'mk'), min)))
    write('metric_cost.tex', '\n'.join(table_rows(lambda a, i: metric(a, i, 'co'), min)))
    def hv_pct(key):
        def g(a, i):
            if a in ('HEFT', 'MINMIN'):
                v = heur[i]['heft' if a == 'HEFT' else 'minmin']
            else:
                v = meanrow(a, i, key)
            allv = [heur[i]['heft'], heur[i]['minmin']] + [meanrow(x, i, key) for x in algs[2:] if meanrow(x, i, key) is not None]
            return None if v is None else 100 * v / max(allv)
        return g
    write('metric_hv_arc.tex', '\n'.join(table_rows(hv_pct('hv_arc'), max, fmt=lambda v: '--' if v is None else f'{v:.2f}')))
    write('metric_hv_pop.tex', '\n'.join(table_rows(hv_pct('hv_pop'), max, fmt=lambda v: '--' if v is None else f'{v:.2f}')))
    write('metric_points.tex', '\n'.join(table_rows(lambda a, i: (1.0 if a in ('HEFT', 'MINMIN') else meanrow(a, i, 'n_unique_front')), max,
                                                   fmt=lambda v: '--' if v is None else f'{v:.1f}', only=algs[2:])))
    write('metric_time.tex', '\n'.join(table_rows(lambda a, i: (None if a in ('HEFT', 'MINMIN') else meanrow(a, i, 'search_ms')), min,
                                                 fmt=lambda v: '--' if v is None else (f'{v:,.0f}'.replace(',', '\\,') if v >= 1000 else f'{v:.1f}'), only=algs[2:])))

# ------------------------------------------------------------------ 5. sensitivity
if glob.glob(str(RESULTS / 'sens_*.jsonl')):
    res = sens_table()
    say('\n== SENSITIVITY of LIWSA+ (archive HV vs default, all 20 workflows, 5 seeds)')
    names = {'pmin002': r'operator-probability floor 0.02', 'pmin005': r'operator-probability floor 0.05', 'pmin020': r'operator-probability floor 0.20',
             'alpha01': r'credit learning rate 0.1', 'alpha05': r'credit learning rate 0.5', 'mf0005': r'move size 0.5\% of tasks',
             'mf005': r'move size 5\% of tasks', 'mf010': r'move size 10\% of tasks', 'pop20': r'population 20', 'pop50': r'population 50',
             'sweep0': r'no sweep seeds', 'sweep4': r'4 sweep seeds', 'sweep16': r'16 sweep seeds'}
    ps = [c['p'] if c else 1.0 for _, c in res]
    hp = A.holm(np.array(ps))
    L = []
    for (n, c), h in zip(res, hp):
        say(f'{n}: {A.fmt(c)} holm={pfmt(h)}')
        L.append(f'{names.get(n, n)} & ' + cell(c, h) + ' \\\\')
    write('sensitivity.tex', '\n'.join(L) + '\n')

# ------------------------------------------------------------------ 5b. expensive-objective (surrogate-assisted) study
SAP = [('LIWSAP.f.r.san-plus', 'LIWSAP.f.r-plus'), ('LIWSAP.f.r.sar-plus', 'LIWSAP.f.r-plus'), ('LIWSAP.f.r.sao-plus', 'LIWSAP.f.r-plus'),
       ('LIWSAP.f.r.sao.hvi-plus', 'LIWSAP.f.r-plus'), ('LIWSAP.f.r.imm-plus', 'LIWSAP.f.r-plus'),
       ('LIWSAP.f.r.san-plus', 'LIWSAP.f.r.sar-plus'), ('LIWSAP.f.r.san-plus', 'LIWSAP.f.r.sao-plus'), ('LIWSAP.f.r.san-plus', 'LIWSAP.f.r.sao.hvi-plus'),
       ('LIWSAP.f.r.san-plus', 'NSGA2.f.r-plus'), ('LIWSAP.f.r.san-plus', 'SPEA2.f.r-plus'), ('LIWSAP.f.r.san-plus', 'MOEAD.f.r-plus'),
       ('LIWSAP.f.r-plus', 'NSGA2.f.r-plus'), ('LIWSAP.f.r-plus', 'SPEA2.f.r-plus'),
       ('LIWSAP.f.r.sal-plus', 'LIWSAP.f.r-plus'), ('LIWSAP.f.r.sag-plus', 'LIWSAP.f.r-plus'), ('LIWSAP.f.r.sal.hvi-plus', 'LIWSAP.f.r-plus'),
       ('LIWSAP.f.r.sag.hvi-plus', 'LIWSAP.f.r-plus'), ('LIWSAP.f.r.sal-plus', 'LIWSAP.f.r.san-plus'), ('LIWSAP.f.r.sag-plus', 'LIWSAP.f.r.san-plus'),
       ('LIWSAP.f.r.sal.hvi-plus', 'LIWSAP.f.r.san-plus'), ('LIWSAP.f.r.sag.hvi-plus', 'LIWSAP.f.r.san-plus'), ('LIWSAP.f.r.sal-plus', 'LIWSAP.f.r.sao-plus'),
       ('LIWSAP.f.r.sag-plus', 'LIWSAP.f.r.sal-plus'), ('LIWSAP.f.r.sal-plus', 'NSGA2.f.r-plus'), ('LIWSAP.f.r.sal-plus', 'SPEA2.f.r-plus'),
       ('LIWSAP.f.r.sal-plus', 'MOEAD.f.r-plus'), ('LIWSAP.f.r.sal.hvi-plus', 'NSGA2.f.r-plus')]
SALAB = {'LIWSAP.f.r.san-plus': 'MF', 'LIWSAP.f.r-plus': 'LIWSA+', 'LIWSAP.f.r.sar-plus': 'random pre-selection', 'LIWSAP.f.r.sao-plus': 'MF + ridge',
         'LIWSAP.f.r.sao.hvi-plus': 'MF + ridge + HVI', 'LIWSAP.f.r.imm-plus': 'OLS immigrants', 'NSGA2.f.r-plus': 'NSGA-II', 'SPEA2.f.r-plus': 'SPEA2',
         'MOEAD.f.r-plus': 'MOEA/D', 'LIWSAP.f.r.sal-plus': 'MF + rich ridge', 'LIWSAP.f.r.sag-plus': 'MF + boosted trees',
         'LIWSAP.f.r.sal.hvi-plus': 'MF + rich ridge + HVI', 'LIWSAP.f.r.sag.hvi-plus': 'MF + boosted trees + HVI'}
sa = {st: sa_rows(st) for st in ('rc', 'rm')}
if all(v for v in sa.values()):
    for sname, sub, tag in (('ALL 20', None, 'sa_main'), ('15 held-out workflows', lambda i: i not in DEV, 'sa_test'), ('five 1000-task workflows', lambda i: i in A.LARGE, 'sa_large')):
        say(f'\n== EXPENSIVE-OBJECTIVE STUDY ({sname}); out-of-sample hypervolume gain (%)')
        tabs = {st: {(x[0], x[1]): x for x in holm_table(sa[st], SAP, 'hv_arc', sub)} for st in sa}
        L = []
        for a_, b_ in SAP:
            parts = []
            for st in ('rc', 'rm'):
                x = tabs[st][(a_, b_)]
                say(f'{st} {SALAB[a_]} vs {SALAB[b_]}: {A.fmt(x[2])} holm={pfmt(x[3])}')
                parts.append('-- & -- & --' if x[2] is None else f"{x[2]['mean']:+.2f} & {x[2]['wins']}/{x[2]['losses']} & {pfmt(x[3])}")
            L.append(f'{SALAB[a_]} vs {SALAB[b_]} & ' + ' & '.join(parts) + ' \\\\')
        write(tag + '.tex', '\n'.join(L))
    L = []
    say('\n== EXPENSIVE-OBJECTIVE STUDY: mean search time (ms), small / 1000-task; expensive evaluations')
    for c in SA_NAMES:
        cells = []
        for st in ('rc', 'rm'):
            rr = sa[st]
            t = [r['search_ms'] for r in rr if r['cond'] == c and r['instance'] not in A.LARGE]
            g = [r['search_ms'] for r in rr if r['cond'] == c and r['instance'] in A.LARGE]
            cells += [f'{np.mean(t):.0f}', f'{np.mean(g):,.0f}'.replace(',', '\\,')]
            say(f'{st} {SALAB[c]}: {np.mean(t):.0f} / {np.mean(g):.0f}; evaluations {np.mean([r["evals"] for r in rr if r["cond"] == c]):.0f}')
        L.append(SA_NAMES[c] + ' & ' + ' & '.join(cells) + ' \\\\')
    write('sa_time.tex', '\n'.join(L))
    # per-workflow table of out-of-sample hypervolume (% of best), both settings
    cols = list(SA_NAMES)
    for st in ('rc', 'rm'):
        means = {c: A.inst_means(sa[st], c, 'hv_arc') for c in cols}
        L = []
        for i in ORDER:
            if i not in means[cols[0]]:
                continue
            best = max(means[c][i] for c in cols)
            cells = []
            for c in cols:
                v = 100 * means[c][i] / best
                t = f'{v:.2f}'
                cells.append('\\textbf{' + t + '}' if t == f'{100.0:.2f}' else t)
            L.append(i.replace('_', '\\_') + ' & ' + ' & '.join(cells) + ' \\\\')
        write(f'metric_sa_{st}.tex', '\n'.join(L))
    # equal wall-clock control
    ep = RESULTS / 'rob_rc_eqtime.jsonl'
    if ep.exists():
        import analyze_surrogate as S2
        rr = [r for r in S2.load([str(RESULTS / 'rob_rc_small.jsonl')])[1] if r['cond'] in ('LIWSAP.f.r-plus', 'LIWSAP.f.r.san-plus')]
        rr += [dict(r, cond='LIWSAP.f.r-plus@23gens') for r in S2.load([str(ep)])[1]]
        rows_eq = S2.build(rr)
        say('\n== EQUAL WALL-CLOCK CONTROL (risk-aware, 15 small workflows)')
        for c in ('LIWSAP.f.r-plus', 'LIWSAP.f.r-plus@23gens', 'LIWSAP.f.r.san-plus'):
            say(f'{c}: evals {np.mean([r["evals"] for r in rows_eq if r["cond"] == c]):.0f} time {np.mean([r["search_ms"] for r in rows_eq if r["cond"] == c]):.0f} ms')
        say('MF vs LIWSA+ with 23 generations:', A.fmt(A.compare(rows_eq, 'LIWSAP.f.r.san-plus', 'LIWSAP.f.r-plus@23gens', 'hv_arc')))
        say('LIWSA+ 23 vs 20 generations:', A.fmt(A.compare(rows_eq, 'LIWSAP.f.r-plus@23gens', 'LIWSAP.f.r-plus', 'hv_arc')))
        eqp = [A.compare(rows_eq, 'LIWSAP.f.r.san-plus', 'LIWSAP.f.r-plus@23gens', 'hv_arc'), A.compare(rows_eq, 'LIWSAP.f.r-plus@23gens', 'LIWSAP.f.r-plus', 'hv_arc')]
        say('Holm over the two equal wall-clock tests:', [pfmt(h) for h in A.holm(np.array([c['p'] for c in eqp]))])
    # development-subset effect, worst-case behaviour and Holm over the equal wall-clock pair
    say('\n== DEVELOPMENT vs HELD-OUT (MF vs LIWSA+, out-of-sample hypervolume gain %)')
    for st in ('rc', 'rm'):
        for nm, sub in (('five development workflows', lambda i: i in DEV), ('15 held-out workflows', lambda i: i not in DEV)):
            c = A.compare(sa[st], 'LIWSAP.f.r.san-plus', 'LIWSAP.f.r-plus', 'hv_arc', sub)
            say(f'{st} {nm}: {A.fmt(c)}')
    say('\n== WORST CASE ACROSS WORKFLOWS: out-of-sample hypervolume as % of the best algorithm on that workflow')
    for st in ('rc', 'rm'):
        means = {c: A.inst_means(sa[st], c, 'hv_arc') for c in SA_NAMES}
        insts = [i for i in ORDER if i in means[list(SA_NAMES)[0]]]
        best = {i: max(means[c][i] for c in SA_NAMES) for i in insts}
        line = []
        for c in SA_NAMES:
            vals = {i: 100 * means[c][i] / best[i] for i in insts}
            w = min(vals, key=vals.get)
            line.append(f'{SALAB[c]}: min {vals[w]:.1f} ({w}), mean {np.mean(list(vals.values())):.1f}')
        say(st, ' | '.join(line))
        g = {i: 100 * (means['LIWSAP.f.r.san-plus'][i] / means['LIWSAP.f.r-plus'][i] - 1) for i in insts}
        loss = sorted((v, i) for i, v in g.items() if v < 0)
        say(f'{st} workflows where MF is below LIWSA+: ' + ', '.join(f'{i} {v:+.2f}%' for v, i in loss))
        ex = [v for i, v in g.items() if i != 'Sipht_30']
        say(f'{st} mean MF gain over LIWSA+ excluding Sipht_30: {np.mean(ex):+.2f}% (all 20: {np.mean(list(g.values())):+.2f}%)')
    for st in ('rc',):
        means = {c: A.inst_means(sa[st], c, 'hv_arc') for c in SA_NAMES}
        insts = [i for i in ORDER if i in means['LIWSAP.f.r.san-plus']]
        d = [100 * (means['LIWSAP.f.r.sao-plus'][i] / means['LIWSAP.f.r.san-plus'][i] - 1) for i in insts if i != 'Sipht_30']
        say(f'{st} mean MF + ridge vs MF excluding Sipht_30: {np.mean(d):+.2f}%')
    say('\n== SIPHT_30 risk-aware, per seed: out-of-sample hypervolume (archive) as % of the best condition on that seed')
    rr = [r for r in sa['rc'] if r['instance'] == 'Sipht_30']
    for sd in sorted({r['seed'] for r in rr}):
        row = {r['cond']: r['hv_arc'] for r in rr if r['seed'] == sd}
        b = max(row.values())
        say(f'seed {sd}: ' + ', '.join(f"{SALAB[c]} {100 * row[c] / b:.1f}" for c in SA_NAMES if c in row))
    say('\n== SEEDED START: share of the final population hypervolume already held at generation 0 (LIWSA+, nominal objective, original prices)')
    def _hv(P, ref):
        P = sorted(set(map(tuple, P))); out = []; best = 1e300
        for m, c in P:
            if c < best: out.append((m, c)); best = c
        out = [q for q in out if q[0] < ref[0] and q[1] < ref[1]]
        h = 0.0; prev = ref[1]
        for m, c in out:
            h += (ref[0] - m) * (prev - c); prev = c
        return h
    for nm, fn in (('small', 'faith_small'), ('large', 'faith_large')):
        fp = RESULTS / f'{fn}.jsonl'
        if not fp.exists(): continue
        byw = {}
        for l in open(fp):
            r = json.loads(l)
            if r.get('type') == 'run' and r['cond'] == 'LIWSAP.f-plus': byw.setdefault(r['instance'], []).append(r)
        sh = []
        for w, rs in byw.items():
            allp = [q for r in rs for c in r['trace'] for q in c['f']]
            ref = (1.2 * max(q[0] for q in allp), 1.2 * max(q[1] for q in allp))
            sh.append(np.mean([_hv(r['trace'][0]['f'], ref) / _hv(r['trace'][-1]['f'], ref) for r in rs]))
        say(f'{nm}: {100 * np.mean(sh):.1f}% (mean over {len(sh)} workflows)')
    # ranking-quality diagnostic
    L = []
    say('\n== SURROGATE RANKING QUALITY (candidate children, 15 small workflows x 2 seeds)')
    for st in ('rc', 'rm'):
        dp = RESULTS / f'rob_{st}_diag.jsonl'
        if dp.exists():
            d = np.array([json.loads(l)['diag'] for l in open(dp) if '"type":"run"' in l])
            m = d.mean(axis=0)
            say(st, [round(float(x), 3) for x in m])
            L.append(f'{SA_SET[st]} & {m[0]:.3f} & {m[1]:.3f} & {m[2]:.3f} & {m[3]:.3f} & {m[4]:.3f} & {m[5]:.3f} \\\\')
    write('sa_diag.tex', '\n'.join(L))
    wfrows = []
    say('\n== RANKING QUALITY PER WORKFLOW (2 seeds): rho makespan nominal/corrected | rho cost nominal/corrected | overlap nominal/corrected')
    for st in ('rc', 'rm'):
        dp = RESULTS / f'rob_{st}_diag.jsonl'
        if dp.exists():
            per = {}
            for l in open(dp):
                if '"type":"run"' in l:
                    r = json.loads(l)
                    per.setdefault(r['instance'], []).append(r['diag'])
            for i in ORDER:
                if i in per:
                    m = np.mean(per[i], axis=0)
                    say(f'{st} {i}: {m[0]:.3f}/{m[1]:.3f} | {m[2]:.3f}/{m[3]:.3f} | {m[4]:.3f}/{m[5]:.3f}')
                    if st == 'rc': wfrows.append(f"{i.replace('_', chr(92) + '_')} & {m[0]:.2f} & {m[1]:.2f} & {m[2]:.2f} & {m[3]:.2f} & {m[4]:.2f} & {m[5]:.2f} \\\\")
    if wfrows:
        write('sa_diag_wf.tex', '\n'.join(wfrows))

# ------------------------------------------------------------------ 6. robustness
if (RESULTS / 'robust.jsonl').exists():
    rt = robust_table()
    say('\n== ROBUSTNESS: mean hypervolume over scenarios (relative to LIWSA+ in the same setting), and degradation from nominal')
    cols = ['NSGA2.f-plus', 'SPEA2.f-plus', 'MOEAD.f-plus', 'LIWSA.f-plus', 'LIWSAP.f-plus', 'HEFT', 'MINMIN']
    lab = dict(NAMES, HEFT='HEFT', MINMIN='Min-Min')
    L = []
    for setting in ('nominal', 'len10', 'len25', 'mips', 'both'):
        base = rt[setting]['LIWSAP.f-plus']
        line = []
        for c in cols:
            rel = [100 * (rt[setting][c][i] / base[i] - 1) for i in base if i in rt[setting][c] and base[i] > 0]
            line.append(np.mean(rel))
        say(setting, {lab[c]: round(v, 2) for c, v in zip(cols, line)})
        L.append(setting + ' & ' + ' & '.join(f'{v:+.2f}' for v in line) + ' \\\\')
    write('robustness.tex', '\n'.join(L) + '\n')
    # paired test LIWSA+ vs NSGA-II per setting
    for setting in ('nominal', 'len25', 'mips', 'both'):
        a = rt[setting]['LIWSAP.f-plus']; b = rt[setting]['NSGA2.f-plus']
        ii = sorted(set(a) & set(b))
        x = np.array([a[i] for i in ii]); y = np.array([b[i] for i in ii])
        from scipy import stats
        p = stats.wilcoxon(x, y).pvalue
        say(f'robust LIWSA+ vs NSGA-II [{setting}]: mean {100 * np.mean((x - y) / y):+.2f}% W/L {int((x > y).sum())}/{int((x < y).sum())} p={p:.4f}')

# ------------------------------------------------------------------ 7. simulator replay
sp = RESULTS / 'sim_faithful_small.csv'
if sp.exists():
    import sim_report
    say('\n== SIMULATOR REPLAY')
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        sim_report.main(str(sp))
    say(buf.getvalue())
OUT.close()
print('wrote', ANALYSIS / 'paper_statistics.txt', 'and', len(list(TAB.glob('*.tex'))), 'tables in', TAB)
