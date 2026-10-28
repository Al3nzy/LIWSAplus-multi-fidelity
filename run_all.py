#!/usr/bin/env python3
"""One command for the whole study: runs every step, then collects the results, tables and figures into RESULTS/ with an index.html.

    python run_all.py            full run (resumable: finished work is skipped, so an interrupted run can simply be started again)
    python run_all.py --quick    smoke test of the whole chain on a few workflows (minutes)
    python run_all.py --collect  only (re)collect the outputs of an earlier run

A step that fails does not stop the others; failures are listed at the end and in RESULTS/run_log.txt.
Windows, Linux and macOS: needs a JDK (javac, java on the PATH, or set JAVAC) and Python 3 with numpy, scipy and matplotlib.
"""
import argparse, html, importlib, os, shutil, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
EXP = ROOT / 'experiments'
RESULTS = ROOT / 'RESULTS'


def collect(work, log):
    ana = work / 'analysis'
    if not ana.exists():
        print('nothing to collect: no analysis folder at', ana)
        return
    for sub in ('figures', 'tables', 'reports'):
        shutil.rmtree(RESULTS / sub, ignore_errors=True)
        (RESULTS / sub).mkdir(parents=True, exist_ok=True)
    for p in ana.rglob('*'):
        if not p.is_file():
            continue
        dest = 'figures' if p.suffix in ('.png', '.pdf', '.svg') else 'tables' if p.suffix == '.tex' else 'reports' if p.suffix in ('.txt', '.csv') else None
        if dest:
            shutil.copy2(p, RESULTS / dest / p.name)
    rows = []
    for name, title in (('figures', 'Figures'), ('tables', 'LaTeX tables'), ('reports', 'Reports')):
        files = sorted((RESULTS / name).glob('*'))
        rows.append(f'<h2>{title} ({len(files)})</h2>')
        if name == 'figures':
            pngs = [f for f in files if f.suffix == '.png']
            rows.append('<div style="display:flex;flex-wrap:wrap;gap:16px">' + ''.join(
                f'<figure style="width:46%;margin:0"><a href="figures/{f.name}"><img src="figures/{f.name}" style="width:100%"></a><figcaption>{html.escape(f.stem)} '
                f'(<a href="figures/{f.stem}.pdf">pdf</a>)</figcaption></figure>' for f in pngs) + '</div>')
        else:
            rows.append('<ul>' + ''.join(f'<li><a href="{name}/{f.name}">{html.escape(f.name)}</a></li>' for f in files) + '</ul>')
    (RESULTS / 'index.html').write_text('<html><head><meta charset="utf-8"><title>Results</title></head><body style="font-family:sans-serif;max-width:1100px;margin:auto">'
                                        f'<h1>Results</h1><p>Collected {time.strftime("%Y-%m-%d %H:%M")} from {html.escape(str(work))}. Main numbers: reports/ (the *_statistics.txt or *_report.txt file).</p>'
                                        + ''.join(rows) + '</body></html>', encoding='utf-8')
    print(f'collected into {RESULTS} (open index.html)')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--collect', action='store_true', help='only collect the outputs of an earlier run')
    ap.add_argument('--steps', nargs='*', help='run only these steps (default: all, in order)')
    a = ap.parse_args()
    if a.quick:
        os.environ['LIWSA_WORK'] = 'work_quick'
    work = EXP / ('work_quick' if a.quick else 'work')
    RESULTS.mkdir(exist_ok=True)
    if a.collect:
        collect(work, None)
        return
    sys.path.insert(0, str(EXP))
    sys.argv = ['run_experiments.py']
    R = importlib.import_module('run_experiments')
    steps = a.steps or list(R.STEPS)
    log = open(RESULTS / 'run_log.txt', 'a', encoding='utf-8')
    log.write(f'\n=== run_all {time.strftime("%Y-%m-%d %H:%M:%S")} quick={a.quick} steps={steps}\n')
    failed, t0 = [], time.time()
    for s in steps:
        t1 = time.time()
        print(f'\n########## step {s}', flush=True)
        cmd = [sys.executable, str(EXP / 'run_experiments.py'), s] + (['--quick'] if a.quick else [])
        proc = subprocess.Popen(cmd, cwd=str(ROOT), stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors='replace')
        for line in proc.stdout:
            if 'Size is negative' in line:
                continue
            sys.stdout.write(line); log.write(line)
        proc.wait()
        msg = f'step {s}: {"ok" if proc.returncode == 0 else "FAILED (exit %d)" % proc.returncode} in {(time.time() - t1) / 60:.1f} min'
        print(msg, flush=True); log.write(msg + '\n'); log.flush()
        if proc.returncode != 0:
            failed.append(s)
    collect(work, log)
    summary = f'finished in {(time.time() - t0) / 60:.1f} min; failed steps: {failed or "none"}'
    print('\n' + summary); log.write(summary + '\n'); log.close()
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
