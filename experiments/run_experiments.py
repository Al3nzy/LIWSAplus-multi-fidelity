#!/usr/bin/env python3
"""One-command driver for the multi-fidelity risk-aware scheduling study. Run from anywhere:

    python experiments/run_experiments.py paper            # everything the paper reports, resumable (about 2.5 h on one core)
    python experiments/run_experiments.py paper --quick    # 2 workflows, 2 seeds (smoke test, a few minutes)

Steps (each can be run alone and is resumable; finished work is skipped):
    compile        compile sources/ + examples/ into experiments/work/classes (needs a JDK; set JAVAC to override 'javac')
    seeds          run HEFT and Min-Min once per workflow inside WorkflowSim and cache their schedules (shared seeds)
    fidelity       decoder-versus-simulator error of the analytic and the simulator-faithful decoder
    faithful       nominal-objective study: LIWSA+ against NSGA-II, SPEA2, MOEA/D and LIWSA, three price lists, ablation
    sensitivity    parameter sensitivity of LIWSA+ (13 settings)
    robustness     planned schedules re-evaluated under runtime-estimate error and degraded machine speed
    multifidelity  expensive-objective study (expected value and CVaR): pre-screening, learned corrections, OLS immigrants,
                   ranking diagnostics, equal wall-clock control
    simcheck       replay archive-front solutions inside WorkflowSim and compare simulator-measured hypervolume
    mlchecks       checks on the learned surrogates: ranking quality, candidate-pool size, harsher tail risk, hedged selection
    figures        statistics, LaTeX tables and figures (paper_statistics.py, make_paper_figures.py, sim_report.py)
    paper          all of the above in this order
"""
import argparse, os, subprocess, sys, time
from pathlib import Path
if '--quick' in sys.argv:
    os.environ['LIWSA_WORK'] = 'work_quick'      # never mix smoke-test artefacts with the real run
from paths import *

ALL = ['Montage_25', 'Montage_50', 'Montage_100', 'Montage_1000',
       'CyberShake_30', 'CyberShake_50', 'CyberShake_100', 'CyberShake_1000',
       'Sipht_30', 'Sipht_60', 'Sipht_100', 'Sipht_1000',
       'Epigenomics_24', 'Epigenomics_46', 'Epigenomics_100', 'Epigenomics_997',
       'Inspiral_30', 'Inspiral_50', 'Inspiral_100', 'Inspiral_1000']

SIM = ['NSGA2.f-plus', 'SPEA2.f-plus', 'MOEAD.f-plus', 'LIWSAP.f-plus', 'LIWSA.f-plus']
# nominal-objective study ('.f' = simulator-faithful decoder, see FaithfulModel.java)
FAITH = ['NSGA2.f-plus', 'SPEA2.f-plus', 'MOEAD.f-plus', 'LIWSA.f-plus', 'LIWSAP.f-plus', 'LIWSAP.noe.f-plus', 'LIWSAP.nog.f-plus',
         'LIWSAP.nol.f-plus', 'LIWSAP.noa.f-plus', 'NSGA2.f-warm', 'LIWSA.f-warm']
FAITH_PREM = ['NSGA2.f-plus', 'SPEA2.f-plus', 'MOEAD.f-plus', 'LIWSA.f-plus', 'LIWSAP.f-plus', 'LIWSAP.nog.f-plus']
# expensive-objective study: Monte-Carlo objectives, budget 610 expensive evaluations (30 + 20 x 29)
SA_CONDS = ['NSGA2.f.r-plus', 'SPEA2.f.r-plus', 'MOEAD.f.r-plus', 'LIWSAP.f.r-plus', 'LIWSAP.f.r.sar-plus', 'LIWSAP.f.r.san-plus',
            'LIWSAP.f.r.sao-plus', 'LIWSAP.f.r.sao.hvi-plus', 'LIWSAP.f.r.imm-plus',
            'LIWSAP.f.r.sal-plus', 'LIWSAP.f.r.sag-plus', 'LIWSAP.f.r.sal.hvi-plus', 'LIWSAP.f.r.sag.hvi-plus']
SA_SETTINGS = {'rc': ['-Dgens=20', '-Drobust=30', '-Dcvar=0.8', '-Dslowp=0.2', '-Dmipslow=0.8'],   # risk-aware: CVaR_0.8 of makespan, stragglers
               'rm': ['-Dgens=20', '-Drobust=20']}                                                  # expected value over 20 scenarios
ROBUST = ['NSGA2.f-plus', 'SPEA2.f-plus', 'MOEAD.f-plus', 'LIWSA.f-plus', 'LIWSAP.f-plus']
# parameter sensitivity of LIWSA+ (one JVM option set per configuration; population size keeps the 2930-evaluation budget)
SENS = [('pmin002', '-Dpmin=0.02'), ('pmin005', '-Dpmin=0.05'), ('pmin020', '-Dpmin=0.20'), ('alpha01', '-Dalpha=0.1'),
        ('alpha05', '-Dalpha=0.5'), ('mf0005', '-Dmovefrac=0.005'), ('mf005', '-Dmovefrac=0.05'), ('mf010', '-Dmovefrac=0.10'),
        ('pop20', '-Dpop=20'), ('pop50', '-Dpop=50'), ('sweep0', '-Dsweep=0'), ('sweep4', '-Dsweep=4'), ('sweep16', '-Dsweep=16')]


def classpath():
    jars = [str(p) for p in sorted((ROOT / 'lib').glob('*.jar'))]
    return os.pathsep.join([str(CLASSES)] + jars)


def sh(cmd, **kw):
    print('>', ' '.join(map(str, cmd)), flush=True)
    r = subprocess.run(cmd, cwd=str(ROOT), **kw)
    if r.returncode != 0:
        sys.exit(f'command failed with exit code {r.returncode}')


def do_compile(a):
    files = [p.as_posix() for base in ('sources', 'examples') for p in (ROOT / base).rglob('*.java')]
    argfile = WORK / 'javac_args.txt'
    argfile.write_text('\n'.join('"%s"' % f for f in files))
    sh([os.environ.get('JAVAC', 'javac'), '-nowarn', '-encoding', 'UTF-8', '-cp', classpath(), '-d', str(CLASSES), '@' + str(argfile)])


def java_main(main, args, jvm=()):
    sh(['java', '-Xmx2g'] + list(jvm) + ['-cp', classpath(), main] + [str(x) for x in args])


def runner(out, instances, seeds, conds, jvm=()):
    java_main('org.workflowsim.examples.experiments.ExperimentRunner', [out, instances, seeds, ','.join(conds), SEEDS], jvm)


def do_seeds(a):
    todo = [w for w in (a.instances or ALL)
            if not ((SEEDS / f'{w}_HEFT.txt').exists() and (SEEDS / f'{w}_MINMIN.txt').exists())]
    if todo:
        java_main('org.workflowsim.examples.experiments.SeedBuilder', [SEEDS, ','.join(todo)])


def do_fidelity(a):
    inst = 'Montage_25,CyberShake_30' if a.quick else 'Montage_25,Montage_50,CyberShake_30,CyberShake_50,Epigenomics_24,Inspiral_30,Sipht_30,Sipht_60'
    java_main('org.workflowsim.examples.experiments.FidelityTest', [RESULTS / 'fidelity.csv', inst, 4 if a.quick else 10, SEEDS])
    sh([sys.executable, str(EXP / 'fidelity_report.py'), str(RESULTS / 'fidelity.csv')])


def do_faithful(a):
    if a.quick:
        runner(RESULTS / 'faith_quick.jsonl', 'Montage_25,CyberShake_30', 2, ['NSGA2.f-plus', 'LIWSAP.f-plus'])
        return
    runner(RESULTS / 'faith_small.jsonl', 'SMALL', a.small_seeds, FAITH)
    runner(RESULTS / 'faith_large.jsonl', 'LARGE', a.large_seeds, FAITH)
    runner(RESULTS / 'faithprem_small.jsonl', 'SMALL', a.small_seeds, FAITH_PREM, ['-Dpricing=premium'])
    runner(RESULTS / 'faithprem_large.jsonl', 'LARGE', a.large_seeds, FAITH_PREM, ['-Dpricing=premium'])
    runner(RESULTS / 'faithlin_small.jsonl', 'SMALL', a.small_seeds, FAITH_PREM, ['-Dpricing=linear'])
    runner(RESULTS / 'faithlin_large.jsonl', 'LARGE', a.large_seeds, FAITH_PREM, ['-Dpricing=linear'])


def do_sensitivity(a):
    if a.quick:
        runner(RESULTS / 'sens_pop20.jsonl', 'Montage_25,CyberShake_30', 2, ['LIWSAP.f-plus'], ['-Dpop=20'])
        return
    for name, opt in SENS:
        runner(RESULTS / f'sens_{name}.jsonl', 'ALL', 5, ['LIWSAP.f-plus'], [opt])


def do_robustness(a):
    if a.quick:
        java_main('org.workflowsim.examples.experiments.RobustnessTest',
                  [RESULTS / 'robust.jsonl', 'Montage_25,CyberShake_30', 1, 'NSGA2.f-plus,LIWSAP.f-plus', 5, SEEDS])
        return
    java_main('org.workflowsim.examples.experiments.RobustnessTest', [RESULTS / 'robust.jsonl', 'ALL', 3, ','.join(ROBUST), 20, SEEDS])


def do_multifidelity(a):
    if a.quick:
        runner(RESULTS / 'rob_rc_small.jsonl', 'Montage_25,CyberShake_30', 2,
               ['LIWSAP.f.r-plus', 'LIWSAP.f.r.san-plus', 'LIWSAP.f.r.sao-plus'], SA_SETTINGS['rc'])
        return
    for st, jvm in SA_SETTINGS.items():
        runner(RESULTS / f'rob_{st}_small.jsonl', 'SMALL', a.sa_small_seeds, SA_CONDS, jvm)
        runner(RESULTS / f'rob_{st}_large.jsonl', 'LARGE', a.sa_large_seeds, SA_CONDS, jvm)
        runner(RESULTS / f'rob_{st}_diag.jsonl', 'SMALL', 2, ['LIWSAP.f.r.sao-plus'], jvm + ['-Ddiag=true'])
    # equal wall-clock control: LIWSA+ without pre-screening, given three more generations (about the pre-screening overhead)
    runner(RESULTS / 'rob_rc_eqtime.jsonl', 'SMALL', a.sa_small_seeds, ['LIWSAP.f.r-plus'], ['-Dgens=23'] + SA_SETTINGS['rc'][1:])


def do_simcheck(a):
    if a.quick:
        java_main('org.workflowsim.examples.experiments.SimValidator',
                  [RESULTS / 'sim_quick.csv', 'Montage_25', 1, ','.join(['NSGA2.f-plus', 'LIWSAP.f-plus']), 3, SEEDS])
        return
    java_main('org.workflowsim.examples.experiments.SimValidator',
              [RESULTS / 'sim_faithful_small.csv', 'SMALL', 3, ','.join(SIM), 5, SEEDS])


def do_mlchecks(a):
    import ml_checks
    ml_checks.run_all(runner, a.quick)
    ml_checks.analyse(a.quick)


def do_figures(a):
    for script in ('paper_statistics.py', 'make_paper_figures.py'):
        r = subprocess.run([sys.executable, str(EXP / script)], cwd=str(EXP))
        print(f'{script} exit {r.returncode}')
    if (RESULTS / 'sim_faithful_small.csv').exists():
        with open(ANALYSIS / 'report_sim.txt', 'w') as f:
            r = subprocess.run([sys.executable, str(EXP / 'sim_report.py')], cwd=str(EXP), stdout=f, stderr=subprocess.STDOUT)
        print(f'wrote {ANALYSIS / "report_sim.txt"} (exit {r.returncode})')


STEPS = {'compile': do_compile, 'seeds': do_seeds, 'fidelity': do_fidelity, 'faithful': do_faithful, 'sensitivity': do_sensitivity,
         'robustness': do_robustness, 'multifidelity': do_multifidelity, 'simcheck': do_simcheck, 'mlchecks': do_mlchecks, 'figures': do_figures}

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('step', choices=list(STEPS) + ['paper'])
    ap.add_argument('--quick', action='store_true', help='smoke test: 2 workflows, 2 seeds')
    ap.add_argument('--small-seeds', type=int, default=10)
    ap.add_argument('--large-seeds', type=int, default=5)
    ap.add_argument('--sa-small-seeds', type=int, default=5)
    ap.add_argument('--sa-large-seeds', type=int, default=3)
    ap.add_argument('--instances', nargs='*', help='restrict the seeds step to these workflows')
    a = ap.parse_args()
    if a.quick:
        a.instances = ['Montage_25', 'CyberShake_30']
    steps = list(STEPS) if a.step == 'paper' else [a.step]
    t0 = time.time()
    for s in steps:
        print(f'\n=== {s} ===', flush=True)
        STEPS[s](a)
    print(f'\nfinished in {(time.time() - t0) / 60:.1f} min')
