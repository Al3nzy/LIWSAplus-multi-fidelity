"""Risk-aware validation: schedules from the scenario-based archives replayed in WorkflowSim under R scenarios of the
out-of-sample set (task-length and VM-speed multipliers applied inside the simulator). Reports (i) decoder-versus-simulator
discrepancy per schedule and scenario, (ii) the CVaR / mean estimates of each replayed schedule from the decoder and from the
simulator, and (iii) whether the comparison of algorithms (hypervolume of the replayed sets in the (risk measure, cost) plane)
changes between the two evaluators. Descriptive: T schedules per run and R scenarios are small samples."""
import csv, collections, math, sys
import numpy as np
from scipy import stats
from paths import *
import analyze as A

MAX_SEED = 3      # replays were planned with three seeds per algorithm; rows of further seeds from interrupted runs are ignored
CONDS = {'NSGA2.f.r-plus': 'NSGA-II', 'LIWSAP.f.r-plus': 'LIWSA+', 'LIWSAP.f.r.san-plus': 'MF'}

def load(st):
    rows = []
    for part in ('small', 'large'):
        p = RESULTS / f'scenreplay_{st}_{part}.csv'
        if p.exists():
            rows += [r for r in csv.DictReader(open(p)) if int(r['seed']) <= MAX_SEED]
    return rows

def risk(v, st):
    v = sorted(v)
    if st == 'rc':
        tail = max(1, math.ceil(0.2 * len(v)))
        return float(np.mean(v[-tail:]))
    return float(np.mean(v))

def analyse(st, out):
    rows = load(st)
    if not rows:
        return None
    sched = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in rows:
        k = (r['instance'], r['cond'], int(r['seed']), int(r['sched']))
        for f in ('dec_makespan', 'dec_cost', 'sim_makespan', 'sim_cost'):
            sched[k][f].append(float(r[f]))
    dm = np.array([abs(float(r['sim_makespan']) - float(r['dec_makespan'])) / float(r['dec_makespan']) for r in rows])
    dc = np.array([abs(float(r['sim_cost']) - float(r['dec_cost'])) / float(r['dec_cost']) for r in rows])
    n_runs = len(set(k[:3] for k in sched)); n_inst = len(set(k[0] for k in sched))
    out.append(f'[{st}] {len(rows)} (schedule, scenario) replays; {len(sched)} schedules from {n_runs} runs on {n_inst} workflows; '
               f'schedules per algorithm ' + str({CONDS.get(c, c): sum(1 for k in sched if k[1] == c) for c in CONDS}))
    q = lambda v, p: 100 * np.percentile(v, p)
    out.append(f'  |sim - dec|/dec makespan: mean {100*dm.mean():.3f}%  median {q(dm,50):.3f}%  95th {q(dm,95):.3f}%  max {100*dm.max():.3f}%')
    out.append(f'  |sim - dec|/dec cost:     mean {100*dc.mean():.3f}%  median {q(dc,50):.3f}%  95th {q(dc,95):.3f}%  max {100*dc.max():.3f}%')
    # per-schedule risk estimates
    est = {}
    for k, d in sched.items():
        est[k] = dict(dec=(risk(d['dec_makespan'], st), float(np.mean(d['dec_cost']))), sim=(risk(d['sim_makespan'], st), float(np.mean(d['sim_cost']))))
    de = np.array([abs(e['sim'][0] - e['dec'][0]) / e['dec'][0] for e in est.values()])
    out.append(f'  per-schedule {"CVaR_0.8" if st=="rc" else "mean"} makespan over the replayed scenarios, |sim - dec|/dec: mean {100*de.mean():.3f}%  max {100*de.max():.3f}%')
    xs = np.array([e['dec'][0] for e in est.values()]); ys = np.array([e['sim'][0] for e in est.values()])
    out.append(f'  Spearman correlation of the per-schedule risk estimate, decoder vs simulator: {stats.spearmanr(xs, ys)[0]:.4f}')
    # hypervolume of replayed sets, decoder vs simulator values of the SAME schedules and scenarios
    pts = {'dec': collections.defaultdict(list), 'sim': collections.defaultdict(list)}
    for (inst, c, s, t), e in est.items():
        for ev in ('dec', 'sim'):
            pts[ev][(inst, c, s)].append(e[ev])
    ref = {}
    for ev in pts:
        for (inst, c, s), p in pts[ev].items():
            o = ref.get(inst, (0, 0)); ref[inst] = (max(o[0], 1.2 * max(x[0] for x in p)), max(o[1], 1.2 * max(x[1] for x in p)))
    hv = {ev: collections.defaultdict(lambda: collections.defaultdict(list)) for ev in pts}
    for ev in pts:
        for (inst, c, s), p in pts[ev].items():
            hv[ev][c][inst].append(A.hv2d(A.nondom(p), ref[inst]))
    res = {}
    for a, b in (('LIWSAP.f.r.san-plus', 'LIWSAP.f.r-plus'), ('LIWSAP.f.r.san-plus', 'NSGA2.f.r-plus'), ('LIWSAP.f.r-plus', 'NSGA2.f.r-plus')):
        line = []
        for ev in ('dec', 'sim'):
            ia = sorted(set(hv[ev][a]) & set(hv[ev][b]))
            x = np.array([np.mean(hv[ev][a][i]) for i in ia]); y = np.array([np.mean(hv[ev][b][i]) for i in ia])
            rel = 100 * (x - y) / y
            p = stats.wilcoxon(x, y).pvalue
            res[(a, b, ev)] = (rel.mean(), np.median(rel), int((rel > 0).sum()), int((rel < 0).sum()), p, len(ia))
            line.append(f'{ev}: {rel.mean():+.2f}% (median {np.median(rel):+.2f}%, W/L {(rel>0).sum()}/{(rel<0).sum()} of {len(ia)}, unadjusted p={p:.3f})')
        out.append(f'  HV of replayed sets, {CONDS[a]} vs {CONDS[b]}: ' + ' | '.join(line))
    return res, dict(dm=dm, dc=dc, de=de)

def main():
    out = []; tex = []
    for st, lab in (('rc', 'Risk-aware'), ('rm', 'Expected value')):
        r = analyse(st, out)
        if r is None:
            continue
        res, d = r
        for (a, b, ev), v in res.items():
            if ev == 'sim':
                dd = res[(a, b, 'dec')]
                tex.append(f'{lab if (a,b)==("LIWSAP.f.r.san-plus","LIWSAP.f.r-plus") else ""} & {CONDS[a]} vs {CONDS[b]} & {dd[0]:+.2f} ({dd[2]}/{dd[3]}) & {v[0]:+.2f} ({v[2]}/{v[3]}) & {v[4]:.3f}')
    text = '\n'.join(out)
    print(text)
    open(str(ANALYSIS) + '/report_scenreplay.txt', 'w').write(text + '\n')
    if tex:
        open(str(ANALYSIS) + '/scenreplay.tex', 'w').write(' \\\\\n'.join(tex))

if __name__ == '__main__':
    main()
