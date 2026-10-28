"""Decoder-vs-simulator error of the analytic decoder and of the simulator-faithful decoder (output of FidelityTest)."""
import sys, pandas as pd, numpy as np
from paths import *
d = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else str(RESULTS / 'fidelity.csv'))
for k, label in (('dec', 'analytic gap-filling decoder'), ('faith', 'simulator-faithful decoder')):
    em = 100 * (d[k + '_makespan'] - d.sim_makespan).abs() / d.sim_makespan
    ec = 100 * (d[k + '_cost'] - d.sim_cost).abs() / d.sim_cost
    print(f'{label:30s} makespan error: mean {em.mean():6.2f}%  median {em.median():5.2f}%  max {em.max():6.1f}% | cost error: mean {ec.mean():5.2f}%  max {ec.max():5.1f}%  (n={len(d)})')
