"""Simulator-level comparison: hypervolume of the (makespan, cost) values MEASURED in WorkflowSim for representative
solutions replayed from each algorithm's archive front (see SimValidator). Same statistics as analyze.py."""
import sys, csv, glob, collections, numpy as np
from scipy import stats
from paths import *
import analyze as A

def main(path):
    rows = list(csv.DictReader(open(path)))
    pts = collections.defaultdict(list); dec = collections.defaultdict(list)
    for r in rows:
        k = (r['instance'], r['cond'], int(r['seed']))
        pts[k].append((float(r['sim_makespan']), float(r['sim_cost'])))
        dec[k].append((float(r['dec_makespan']), float(r['dec_cost'])))
    ref = {}
    for (inst, c, s), p in pts.items():
        m = max(x[0] for x in p); co = max(x[1] for x in p)
        o = ref.get(inst, (0, 0)); ref[inst] = (max(o[0], 1.2 * m), max(o[1], 1.2 * co))
    out = []
    for (inst, c, s), p in pts.items():
        out.append(dict(instance=inst, cond=c, seed=s, hv=A.hv2d(A.nondom(p), ref[inst]), hv_dec=A.hv2d(A.nondom(dec[(inst, c, s)]), ref[inst])))
    rr = [dict(r, hv_pop=r['hv'], hv_arc=r['hv']) for r in out]
    conds = sorted(set(r['cond'] for r in rr))
    print('simulator-replayed fronts:', len(rr), 'runs,', len(set(r['instance'] for r in rr)), 'instances, conditions', conds)
    base = [c for c in conds if c.startswith('NSGA2')][0] if any(c.startswith('NSGA2') for c in conds) else conds[0]
    print('\nSimulator-measured hypervolume, paired over instances (A vs B)')
    for a in conds:
        for b in conds:
            if a != b and (b.startswith('NSGA2') or b.startswith('SPEA2') or b.startswith('MOEAD') or b.startswith('LIWSA.') or b.startswith('LIWSA-')) and a.startswith('LIWSAP'):
                print(f'{a:>16} vs {b:<14}', A.fmt(A.compare(rr, a, b, 'hv_pop')))
    # decoder vs simulator agreement of the ranking of algorithms
    print('\nPlanning-level (decoder) HV of the SAME replayed points, for reference')
    for a in conds:
        for b in conds:
            if a != b and a.startswith('LIWSAP') and b.startswith('NSGA2'):
                rd = [dict(r, hv_pop=r['hv_dec']) for r in rr]
                print(f'{a:>16} vs {b:<14}', A.fmt(A.compare(rd, a, b, 'hv_pop')))
    dm = [abs(float(r['sim_makespan']) - float(r['dec_makespan'])) / float(r['dec_makespan']) for r in rows]
    dc = [abs(float(r['sim_cost']) - float(r['dec_cost'])) / float(r['dec_cost']) for r in rows]
    print('\nmedian |simulator - decoder| / decoder: makespan %.2f%%, cost %.2f%% (n=%d replayed schedules)' % (100 * np.median(dm), 100 * np.median(dc), len(rows)))

if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else str(RESULTS / 'sim_faithful_small.csv'))
