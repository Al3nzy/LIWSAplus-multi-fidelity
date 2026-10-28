"""Draws every result figure of the paper (vector PDF + PNG) from the experiment outputs.
Usage: python experiments/make_paper_figures.py [--outdir DIR]"""
import sys, json, argparse, collections
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from paths import *
import analyze as A
from paperlib import *

ap = argparse.ArgumentParser()
ap.add_argument('--outdir', default=str(ANALYSIS / 'figures'))
args = ap.parse_args()
OUTD = Path(args.outdir); OUTD.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.size': 7, 'font.family': 'serif', 'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'legend.fontsize': 6.5, 'axes.labelsize': 7.5, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'lines.linewidth': 0.9})
# Okabe-Ito colours + distinct markers (colour AND shape coding)
STY = {'LIWSAP.f-plus': ('#000000', 'o', 'LIWSA+'), 'NSGA2.f-plus': ('#0072B2', 's', 'NSGA-II'), 'SPEA2.f-plus': ('#E69F00', '^', 'SPEA2'),
       'MOEAD.f-plus': ('#009E73', 'D', 'MOEA/D'), 'LIWSA.f-plus': ('#D55E00', 'x', 'LIWSA (original moves)'),
       'HEFT': ('#CC79A7', 'v', 'HEFT'), 'MINMIN': ('#56B4E9', 'P', 'Min-Min')}


def save(fig, name):
    fig.savefig(OUTD / f'{name}.pdf', bbox_inches='tight')
    fig.savefig(OUTD / f'{name}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print('wrote', name)


def short(i):
    f, n = i.split('_')
    return {'Montage': 'Mon', 'CyberShake': 'Cyb', 'Sipht': 'Sip', 'Epigenomics': 'Epi', 'Inspiral': 'Ins'}[f] + n



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

# ---------------------------------------------------------------- Fig: fidelity
if (RESULTS / 'fidelity.csv').exists():
    d = pd.read_csv(RESULTS / 'fidelity.csv')
    famc = {'Montage': '#0072B2', 'CyberShake': '#D55E00', 'Sipht': '#009E73', 'Epigenomics': '#E69F00', 'Inspiral': '#CC79A7'}
    fm = {'Montage': 'o', 'CyberShake': 's', 'Sipht': '^', 'Epigenomics': 'D', 'Inspiral': 'v'}
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.7))
    for k, (key, ttl) in enumerate((('dec', '(a) Analytic gap-filling decoder'), ('faith', '(b) Simulator-faithful decoder'))):
        for f in famc:
            g = d[d.instance.str.startswith(f)]
            ax[k].scatter(g.sim_makespan, g[key + '_makespan'], s=7, c=famc[f], marker=fm[f], label=f, linewidths=0.3, alpha=0.85)
        lim = [d.sim_makespan.min() * 0.8, d.sim_makespan.max() * 1.2]
        ax[k].plot(lim, lim, 'k-', lw=0.5)
        ax[k].set_xscale('log'); ax[k].set_yscale('log'); ax[k].set_xlabel('Makespan measured in WorkflowSim (s)')
        ax[k].set_ylabel('Makespan predicted by the decoder (s)'); ax[k].set_title(ttl, fontsize=7.5)
        e = 100 * (d[key + '_makespan'] - d.sim_makespan).abs() / d.sim_makespan
        ax[k].text(0.04, 0.9, f'mean error {e.mean():.2f}%\nmax error {e.max():.1f}%', transform=ax[k].transAxes, fontsize=6.5, va='top')
    ax[0].legend(loc='lower right', frameon=False, markerscale=1.3)
    save(fig, 'fig_fidelity')

# ---------------------------------------------------------------- Fig: per-workflow differences, three pricing schemes
store = {}
for sname, prefix in SCEN.items():
    if files(prefix):
        store[sname] = scenario(prefix)
if store:
    fig, axs = plt.subplots(1, len(store), figsize=(7.2, 2.5), sharey=False)
    axs = np.atleast_1d(axs)
    for ax, (sname, (rows, heur, ref)) in zip(axs, store.items()):
        a = A.inst_means(rows, 'LIWSAP.f-plus', 'hv_arc')
        for c in ('NSGA2.f-plus', 'SPEA2.f-plus', 'MOEAD.f-plus', 'LIWSA.f-plus'):
            b = A.inst_means(rows, c, 'hv_arc')
            xs = [k for k, i in enumerate(ORDER) if i in a and i in b]
            ys = [100 * (a[ORDER[k]] / b[ORDER[k]] - 1) for k in xs]
            col, mk, lab = STY[c]
            ax.plot(xs, ys, marker=mk, color=col, ls='none', ms=3.2 if mk != 'x' else 3.5, label=lab, mew=0.7, mfc=col if mk != 'x' else 'none')
        ax.axhline(0, color='k', lw=0.5)
        ax.set_yscale('symlog', linthresh=0.2)
        ax.set_xticks(range(len(ORDER))); ax.set_xticklabels([short(i) for i in ORDER], rotation=90, fontsize=5)
        ax.set_title(f'{sname} pricing', fontsize=7.5)
    axs[0].set_ylabel('LIWSA+ hypervolume gain (%)')
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, ['vs ' + x for x in l], loc='upper center', ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.07))
    fig.tight_layout(); save(fig, 'fig_main_diff')

# ---------------------------------------------------------------- Fig: heat map of % of best hypervolume (original pricing)
if 'original' in store:
    rows, heur, ref = store['original']
    cols = ['NSGA2.f-plus', 'SPEA2.f-plus', 'MOEAD.f-plus', 'LIWSA.f-plus', 'LIWSAP.f-plus']
    M = {c: A.inst_means(rows, c, 'hv_arc') for c in cols}
    mat = np.array([[100 * M[c][i] / max(M[x][i] for x in cols) for c in cols] for i in ORDER if i in M[cols[0]]])
    names = [i for i in ORDER if i in M[cols[0]]]
    fig, ax = plt.subplots(figsize=(3.5, 4.4))
    im = ax.imshow(mat, aspect='auto', cmap='Greys_r', vmin=94, vmax=100)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([STY[c][2].replace(' (original moves)', '\n(orig.)') for c in cols], rotation=45, ha='right')
    ax.set_yticks(range(len(names))); ax.set_yticklabels([n.replace('_', ' ') for n in names], fontsize=5.5)
    for r in range(mat.shape[0]):
        for c in range(mat.shape[1]):
            ax.text(c, r, f'{mat[r, c]:.1f}', ha='center', va='center', fontsize=5, color='white' if mat[r, c] < 97 else 'black')
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02, label='% of best hypervolume')
    save(fig, 'fig_heatmap')

# ---------------------------------------------------------------- Fig: convergence
if 'original' in store:
    rows, heur, ref = store['original']
    runs = raw_runs('faith_')
    conds = ['LIWSAP.f-plus', 'NSGA2.f-plus', 'SPEA2.f-plus', 'MOEAD.f-plus', 'LIWSA.f-plus']
    curves = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in runs:
        if r['cond'] in conds and r['trace']:
            seen = {c['g']: (c['e'], A.hv2d(c['f'], ref[r['instance']])) for c in r['trace']}
            gs = sorted(seen)
            curves[(r['cond'], r['instance'])][0].append([seen[g][0] for g in gs])
            curves[(r['cond'], r['instance'])][1].append([seen[g][1] for g in gs])
    base = {}
    for i in ORDER:
        v = [x[-1] for x in curves[('LIWSAP.f-plus', i)][1]] if ('LIWSAP.f-plus', i) in curves else []
        if v:
            base[i] = np.mean(v)
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.6))
    for panel, subset in enumerate((lambda i: i not in A.LARGE, lambda i: i in A.LARGE)):
        for c in conds:
            xs = None; vals = []
            for i in base:
                if not subset(i) or (c, i) not in curves:
                    continue
                e = np.mean(curves[(c, i)][0], axis=0)
                h = np.mean(curves[(c, i)][1], axis=0) / base[i] * 100
                xs = e; vals.append(h)
            if xs is None:
                continue
            col, mk, lab = STY[c]
            ax[panel].plot(xs[1:], np.median(vals, axis=0)[1:], color=col, marker=mk if mk != 'x' else 'x', ms=2.6, label=lab, mew=0.6, mfc=col if mk != 'x' else 'none')
        ax[panel].set_xscale('log'); ax[panel].set_xlabel('Objective evaluations'); ax[panel].set_ylabel('Median hypervolume (% of final LIWSA+)')
        ax[panel].set_title(['(a) 24-100 task workflows', '(b) ~1000 task workflows'][panel], fontsize=7.5)
        ax[panel].set_ylim(92 if panel == 0 else 70, 100.6)
    ax[0].legend(frameon=False, loc='lower right')
    fig.tight_layout(); save(fig, 'fig_convergence')

# ---------------------------------------------------------------- Fig: example fronts
if 'original' in store:
    rows, heur, ref = store['original']
    runs = {(r['instance'], r['cond'], r['seed']): r for r in raw_runs('faith_')}
    refs, _ = A.load(files('faith_'))
    show = ['CyberShake_100', 'Inspiral_100', 'Epigenomics_100', 'Montage_100']
    fig, axs = plt.subplots(2, 2, figsize=(7.1, 4.6))
    for ax, inst in zip(axs.ravel(), show):
        for c in ('NSGA2.f-plus', 'SPEA2.f-plus', 'MOEAD.f-plus', 'LIWSA.f-plus', 'LIWSAP.f-plus'):
            r = runs.get((inst, c, 1))
            if not r:
                continue
            pts = A.nondom(r['archive'])
            col, mk, lab = STY[c]
            ax.plot([p[0] for p in pts], [p[1] for p in pts], marker=mk, color=col, ms=3, lw=0.6, label=lab, mfc=col if mk != 'x' else 'none', mew=0.6)
        ax.plot(*refs[inst]['heft_f'], marker='v', color=STY['HEFT'][0], ls='none', ms=4, label='HEFT')
        ax.plot(*refs[inst]['minmin_f'], marker='P', color=STY['MINMIN'][0], ls='none', ms=4, label='Min-Min')
        ax.set_xlabel('Makespan (s)'); ax.set_ylabel('Cost'); ax.set_title(inst.replace('_', ' '), fontsize=7.5)
    axs[0, 0].legend(frameon=False, fontsize=5.5)
    fig.tight_layout(); save(fig, 'fig_fronts')

# ---------------------------------------------------------------- Fig: ablation + sensitivity
if 'original' in store:
    rows = store['original'][0]
    abl = [('LIWSAP.f-plus', 'LIWSAP.nog.f-plus', 'guided moves'), ('LIWSAP.f-plus', 'LIWSAP.noe.f-plus', 'estimate move'),
           ('LIWSAP.f-plus', 'LIWSAP.nol.f-plus', 'locust moves'), ('LIWSAP.f-plus', 'LIWSAP.noa.f-plus', 'online operator learning'),
           ('NSGA2.f-plus', 'NSGA2.f-warm', 'sweep seeds (NSGA-II)'), ('LIWSA.f-plus', 'LIWSA.f-warm', 'sweep seeds (LIWSA)')]
    res = [(lab, A.compare(rows, a, b, 'hv_arc')) for a, b, lab in abl]
    ps = A.holm(np.array([c['p'] for _, c in res]))
    fig, ax = plt.subplots(figsize=(3.5, 2.3))
    for k, ((lab, c), h) in enumerate(zip(res, ps)):
        ax.errorbar(c['mean'], k, xerr=[[c['mean'] - c['lo']], [c['hi'] - c['mean']]], fmt='o', color='k', ms=3.5, capsize=2, mfc='k' if h < 0.05 else 'white')
    ax.axvline(0, color='k', lw=0.5); ax.set_yticks(range(len(res))); ax.set_yticklabels([l for l, _ in res]); ax.invert_yaxis()
    ax.set_xlabel('Hypervolume gain from the component (%), mean and 95% CI')
    ax.text(0.99, 0.02, 'filled: Holm p < 0.05', transform=ax.transAxes, ha='right', fontsize=5.5)
    save(fig, 'fig_ablation')
if glob.glob(str(RESULTS / 'sens_*.jsonl')):
    res = sens_table()
    lab = {'pmin002': 'floor 0.02', 'pmin005': 'floor 0.05', 'pmin020': 'floor 0.20', 'alpha01': 'learning rate 0.1', 'alpha05': 'learning rate 0.5',
           'mf0005': 'move size 0.5%', 'mf005': 'move size 5%', 'mf010': 'move size 10%', 'pop20': 'population 20', 'pop50': 'population 50',
           'sweep0': 'no sweep seeds', 'sweep4': '4 sweep seeds', 'sweep16': '16 sweep seeds'}
    ps = A.holm(np.array([c['p'] for _, c in res]))
    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    for k, ((n, c), h) in enumerate(zip(res, ps)):
        ax.errorbar(c['mean'], k, xerr=[[c['mean'] - c['lo']], [c['hi'] - c['mean']]], fmt='o', color='k', ms=3.5, capsize=2, mfc='k' if h < 0.05 else 'white')
    ax.axvline(0, color='k', lw=0.5); ax.set_yticks(range(len(res))); ax.set_yticklabels([lab.get(n, n) for n, _ in res]); ax.invert_yaxis()
    ax.set_xlabel('Hypervolume change vs. default (%), mean and 95% CI')
    save(fig, 'fig_sensitivity')

# ---------------------------------------------------------------- Fig: robustness
if (RESULTS / 'robust.jsonl').exists():
    rt = robust_table()
    settings = ['nominal', 'len10', 'len25', 'mips', 'both']
    lab = {'nominal': 'nominal', 'len10': 'length error\n$\\sigma$=0.10', 'len25': 'length error\n$\\sigma$=0.25', 'mips': 'MIPS 60-100%\nof peak', 'both': 'both'}
    conds = ['NSGA2.f-plus', 'SPEA2.f-plus', 'MOEAD.f-plus', 'LIWSA.f-plus', 'LIWSAP.f-plus']
    fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.5), gridspec_kw={'width_ratios': [2.2, 1]})
    for k, c in enumerate(conds):
        ys = []
        for s in settings:
            base = rt[s]['LIWSAP.f-plus']
            ys.append(np.mean([100 * (rt[s][c][i] / base[i] - 1) for i in base if i in rt[s][c]]))
        col, mk, l = STY[c]
        ax[0].plot(range(len(settings)), ys, marker=mk, color=col, label=l, ms=3.5, mfc=col if mk != 'x' else 'none', mew=0.7)
    ax[0].axhline(0, color='k', lw=0.5); ax[0].set_xticks(range(len(settings))); ax[0].set_xticklabels([lab[s] for s in settings])
    ax[0].set_ylabel('Mean HV relative to LIWSA+ (%)'); ax[0].legend(frameon=False, loc='lower left'); ax[0].set_title('(a) Search algorithms', fontsize=7.5)
    for c in ('HEFT', 'MINMIN'):
        ys = []
        for s in settings:
            base = rt[s]['LIWSAP.f-plus']
            ys.append(np.mean([100 * (rt[s][c][i] / base[i] - 1) for i in base if i in rt[s][c]]))
        col, mk, l = STY[c]
        ax[1].plot(range(len(settings)), ys, marker=mk, color=col, label=l, ms=3.5)
    ax[1].set_xticks(range(len(settings))); ax[1].set_xticklabels(['nom.', 'len10', 'len25', 'mips', 'both'], rotation=45)
    ax[1].legend(frameon=False); ax[1].set_title('(b) Heuristics', fontsize=7.5)
    fig.tight_layout(); save(fig, 'fig_robust')

# ---------------------------------------------------------------- Fig: framework (drawn with matplotlib)
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
def _box(ax, x, y, w, h, text, fc='#FFFFFF', ec='#222222', fs=6.6, bold=False, ls='-', lw=0.8):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.25,rounding_size=1.2', fc=fc, ec=ec, lw=lw, ls=ls))
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=fs, fontweight='bold' if bold else 'normal', linespacing=1.25)
def _arrow(ax, p, q, text=None, off=(0, 1.3), rad=0.0, fs=6):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle='-|>', mutation_scale=7, lw=0.9, color='#222222', connectionstyle=f'arc3,rad={rad}'))
    if text:
        ax.text((p[0] + q[0]) / 2 + off[0], (p[1] + q[1]) / 2 + off[1], text, ha='center', va='bottom', fontsize=fs, style='italic')
fig, ax = plt.subplots(figsize=(7.2, 3.7))
ax.set_xlim(0, 100); ax.set_ylim(-6, 51); ax.axis('off')
BLUE, GREY, GREEN, SAND = '#DCE9F7', '#EFEFEF', '#DDF0E4', '#FBEFD5'
_box(ax, 1.5, 32, 15, 9, 'Workflow DAG\ntasks, files', GREY, fs=5.6)
_box(ax, 1.5, 13, 15, 9, 'VM pool, prices\nMIPS, bandwidth', GREY, fs=5.6)
_box(ax, 21, 5, 29, 39, '', BLUE)
ax.text(35.5, 40.8, 'Simulator-faithful\ndecoder', ha='center', va='center', fontsize=6.6, fontweight='bold', linespacing=1.1)
for k, t in enumerate(['FIFO queue per VM', 'stage-in occupies and\nbills the VM', 'per-VM file cache']):
    _box(ax, 23.2, 31.6 - 6.4 * k, 24.6, 4.8, t, '#FFFFFF', fs=5.0, lw=0.5)
_box(ax, 23.2, 11.3, 24.6, 5.0, 'nominal mode (cheap):\none decode', '#E8F3EA', fs=4.8, lw=0.5)
_box(ax, 23.2, 5.6, 24.6, 4.6, 'expensive mode: mean or CVaR\nover S scenarios', '#F7E3E3', fs=4.6, lw=0.5)
_arrow(ax, (16.9, 36.5), (20.6, 36.5)); _arrow(ax, (16.9, 17.5), (20.6, 17.5))
ax.add_patch(FancyBboxPatch((56, 4), 43, 44, boxstyle='round,pad=0.25,rounding_size=1.5', fc='none', ec='#555555', lw=1.0, ls=(0, (4, 2))))
ax.text(98, 5.4, 'LIWSA+ search', ha='right', va='bottom', fontsize=6.6, fontweight='bold')
_box(ax, 59.5, 36, 17, 8, 'Initial population\nHEFT, Min-Min,\n8 list-scheduling seeds', SAND, fs=4.9)
_box(ax, 59.5, 22, 17, 7.2, 'Survival\nrank and crowding', GREEN, fs=5.4)
_box(ax, 59.5, 9.5, 17, 7.2, 'Anytime archive\nnon-dominated set', GREEN, fs=5.4)
_box(ax, 80, 32.5, 16.5, 11.5, 'Variation pool\ncrossover, original\nLIWSA moves, critical-\npath, cost, staging-\naware; chosen by\nonline operator selection', SAND, fs=4.6)
_box(ax, 80, 21, 16.5, 7.5, 'Multi-fidelity pre-screening\nrank K x lambda children\nby nominal decode,\nkeep the best lambda', '#DDF0E4', fs=4.5)
_box(ax, 80, 9.5, 16.5, 7.2, 'Optional online ridge\ncorrection of the\nnominal objective', '#F2F2F2', fs=4.6)
_arrow(ax, (68, 35.6), (68, 29.6)); _arrow(ax, (68, 21.6), (68, 17.1))
_arrow(ax, (76.9, 28.8), (79.6, 33.5)); ax.text(76.2, 31.4, 'parents', fontsize=4.8, style='italic', ha='right')
_arrow(ax, (88.2, 32.1), (88.2, 28.9)); ax.text(89.2, 30.0, 'candidates', fontsize=4.6, style='italic')
_arrow(ax, (88.2, 17.1), (88.2, 20.6)); ax.text(89.2, 18.6, 'correction', fontsize=4.6, style='italic')
ax.annotate('', xy=(59.1, 25.6), xytext=(50.4, 25.6), arrowprops=dict(arrowstyle='-|>', lw=0.9, color='#222222', mutation_scale=7))
ax.text(54.7, 27.2, 'makespan,\ncost', fontsize=4.6, style='italic', ha='center', va='bottom', linespacing=1.1)
ax.plot([96.8, 98.1, 98.1, 35.5], [24.7, 24.7, 46.4, 46.4], color='#222222', lw=0.9)
ax.add_patch(FancyArrowPatch((35.5, 46.4), (35.5, 44.4), arrowstyle='-|>', mutation_scale=7, lw=0.9, color='#222222'))
ax.text(66, 46.8, 'selected children (genotypes)', ha='center', fontsize=5.0, style='italic', va='bottom')
_box(ax, 59.5, -5.2, 17, 5.2, 'Replay in WorkflowSim\nvalidation of the front', GREY, fs=4.9)
_arrow(ax, (68, 9.1), (68, 0.4))
save(fig, 'fig_framework')

# ---------------------------------------------------------------- Fig: makespan and cost as bar charts
if 'original' in store:
    rows, heur, ref = store['original']
    runs = raw_runs('faith_')
    refs_, _ = A.load(files('faith_'))
    ALG = ['HEFT', 'MINMIN', 'NSGA2.f-plus', 'SPEA2.f-plus', 'MOEAD.f-plus', 'LIWSA.f-plus', 'LIWSAP.f-plus']
    LAB = {'HEFT': 'HEFT', 'MINMIN': 'Min-Min', 'NSGA2.f-plus': 'NSGA-II', 'SPEA2.f-plus': 'SPEA2', 'MOEAD.f-plus': 'MOEA/D',
           'LIWSA.f-plus': 'LIWSA (orig.)', 'LIWSAP.f-plus': 'LIWSA+'}
    COL = {'HEFT': '#CC79A7', 'MINMIN': '#56B4E9', 'NSGA2.f-plus': '#0072B2', 'SPEA2.f-plus': '#E69F00', 'MOEAD.f-plus': '#009E73',
           'LIWSA.f-plus': '#D55E00', 'LIWSAP.f-plus': '#000000'}
    HAT = {'HEFT': '//', 'MINMIN': '\\\\', 'NSGA2.f-plus': 'xx', 'SPEA2.f-plus': '..', 'MOEAD.f-plus': '++', 'LIWSA.f-plus': 'oo', 'LIWSAP.f-plus': ''}
    acc = collections.defaultdict(lambda: collections.defaultdict(list))
    for r in runs:
        if r['archive']:
            acc[(r['cond'], r['instance'])]['mk'].append(min(p[0] for p in r['archive']))
            acc[(r['cond'], r['instance'])]['co'].append(min(p[1] for p in r['archive']))
    def val(a, inst, what):
        if a in ('HEFT', 'MINMIN'):
            r = refs_[inst]['heft_f' if a == 'HEFT' else 'minmin_f']
            return r[0] if what == 'mk' else r[1]
        v = acc.get((a, inst), {}).get(what)
        return float(np.mean(v)) if v else np.nan
    def bars(groups, fname, titles):
        fig, axs = plt.subplots(2, len(groups), figsize=(7.2, 4.3))
        for col, (gname, insts) in enumerate(groups):
            for row, what in enumerate(('mk', 'co')):
                ax = axs[row, col]
                n = len(ALG); w = 0.8 / n
                for k, a in enumerate(ALG):
                    ys = [val(a, i, what) for i in insts]
                    ax.bar(np.arange(len(insts)) + (k - (n - 1) / 2) * w, ys, w, color='white' if a != 'LIWSAP.f-plus' else '#000000', edgecolor=COL[a],
                           hatch=HAT[a], lw=0.6, label=LAB[a])
                ax.set_xticks(range(len(insts))); ax.set_xticklabels([i.split('_')[1] for i in insts], fontsize=6)
                if row == 0:
                    ax.set_title(gname, fontsize=7.5)
                if col == 0:
                    ax.set_ylabel('Fastest solution: makespan (s)' if what == 'mk' else 'Cheapest solution: cost')
                ax.ticklabel_format(axis='y', style='sci', scilimits=(-2, 4), useMathText=True)
                ax.tick_params(axis='y', labelsize=5.5)
        h, l = axs[0, 0].get_legend_handles_labels()
        fig.legend(h, l, loc='upper center', ncol=7, frameon=False, bbox_to_anchor=(0.5, 1.03), fontsize=6.3, handlelength=1.4, columnspacing=1.0)
        fig.tight_layout(rect=(0, 0, 1, 0.95)); save(fig, fname)
    small = [(f, [i for i in ORDER if i.startswith(f) and i not in A.LARGE]) for f in FAMS]
    large = [(f, [i for i in ORDER if i.startswith(f) and i in A.LARGE]) for f in FAMS]
    bars(small, 'fig_bars_small', None)
    bars(large, 'fig_bars_large', None)

print('figures in', OUTD)


# ---------------------------------------------------------------- Figs: expensive-objective (surrogate-assisted) study
_sa = {st: sa_rows(st) for st in ('rc', 'rm')}
if all(v for v in _sa.values()):
    order = ['LIWSAP.f.r-plus', 'LIWSAP.f.r.sar-plus', 'LIWSAP.f.r.san-plus', 'LIWSAP.f.r.sao-plus', 'LIWSAP.f.r.sao.hvi-plus', 'LIWSAP.f.r.imm-plus',
             'SPEA2.f.r-plus', 'MOEAD.f.r-plus']
    lab = {'LIWSAP.f.r-plus': 'LIWSA+', 'LIWSAP.f.r.sar-plus': 'random pre-selection', 'LIWSAP.f.r.san-plus': 'MF (nominal pre-screening)',
           'LIWSAP.f.r.sao-plus': 'MF + online ridge', 'LIWSAP.f.r.sao.hvi-plus': 'MF + ridge + HVI rule', 'LIWSAP.f.r.imm-plus': 'LIWSA+ + OLS immigrants',
           'SPEA2.f.r-plus': 'SPEA2', 'MOEAD.f.r-plus': 'MOEA/D'}
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.7), sharey=True)
    for ax, st in zip(axs, ('rc', 'rm')):
        res = [(c, A.compare(_sa[st], c, 'NSGA2.f.r-plus', 'hv_arc')) for c in order]
        hp = A.holm(np.array([r[1]['p'] for r in res]))
        for k, ((c, r), h) in enumerate(zip(res, hp)):
            col = '#000000' if 'san' in c or 'sao' in c else '#555555'
            ax.errorbar(r['mean'], k, xerr=[[max(r['mean'] - r['lo'], 0)], [max(r['hi'] - r['mean'], 0)]], fmt='o' if h < 0.05 else 'o', color=col, ms=3.5,
                        capsize=2, mfc=col if h < 0.05 else 'white')
        ax.axvline(0, color='k', lw=0.5)
        ax.set_yticks(range(len(order))); ax.set_yticklabels([lab[c] for c in order]); ax.invert_yaxis()
        ax.set_xlabel('Out-of-sample hypervolume gain over NSGA-II (%)')
        ax.set_title({'rc': '(a) Risk-aware objective (CVaR, stragglers)', 'rm': '(b) Expected-value objective'}[st], fontsize=7.5)
    axs[0].text(0.01, 0.01, 'filled: Holm p < 0.05', transform=axs[0].transAxes, fontsize=5.5)
    fig.tight_layout(); save(fig, 'fig_sa_forest')

# ---------------------------------------------------------------- Fig: per-workflow gain of MF pre-screening
if all(v for v in _sa.values()):
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.5), sharey=False)
    for ax, st in zip(axs, ('rc', 'rm')):
        a = A.inst_means(_sa[st], 'LIWSAP.f.r.san-plus', 'hv_arc'); b = A.inst_means(_sa[st], 'LIWSAP.f.r-plus', 'hv_arc'); c_ = A.inst_means(_sa[st], 'NSGA2.f.r-plus', 'hv_arc')
        xs = [k for k, i in enumerate(ORDER) if i in a and i in b]
        ax.plot(xs, [100 * (a[ORDER[k]] / b[ORDER[k]] - 1) for k in xs], 'o', color='k', ms=3.4, label='MF vs LIWSA+')
        ax.plot(xs, [100 * (a[ORDER[k]] / c_[ORDER[k]] - 1) for k in xs], 's', color='#0072B2', ms=3.0, mfc='white', label='MF vs NSGA-II')
        ax.axhline(0, color='k', lw=0.5); ax.set_yscale('symlog', linthresh=1.0)
        ax.set_xticks(range(len(ORDER))); ax.set_xticklabels([short(i) for i in ORDER], rotation=90, fontsize=5)
        ax.set_title({'rc': '(a) Risk-aware objective', 'rm': '(b) Expected-value objective'}[st], fontsize=7.5)
    axs[0].set_ylabel('Hypervolume gain (%)'); axs[0].legend(frameon=False, loc='upper left')
    fig.tight_layout(); save(fig, 'fig_sa_perworkflow')

# ---------------------------------------------------------------- Fig: ranking quality of the surrogates
if all((RESULTS / f'rob_{st}_diag.jsonl').exists() for st in ('rc', 'rm')):
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.3))
    names = ['makespan rank corr.', 'cost rank corr.', 'top-set overlap']
    idx = [(0, 1), (2, 3), (4, 5)]
    for ax, nm, (i0, i1) in zip(axs, names, idx):
        for k, st in enumerate(('rc', 'rm')):
            d = np.array([json.loads(l)['diag'] for l in open(RESULTS / f'rob_{st}_diag.jsonl') if '"type":"run"' in l]).mean(axis=0)
            ax.bar(k - 0.18, d[i0], 0.34, color='white', edgecolor='#0072B2', hatch='xx', label='nominal decode' if k == 0 else None)
            ax.bar(k + 0.18, d[i1], 0.34, color='white', edgecolor='#D55E00', hatch='..', label='with ridge correction' if k == 0 else None)
        ax.set_xticks([0, 1]); ax.set_xticklabels(['risk-aware', 'expected value']); ax.set_title(nm, fontsize=7.5); ax.set_ylim(0, 1.05)
    axs[0].legend(frameon=False, fontsize=5.5, loc='lower left')
    fig.tight_layout(); save(fig, 'fig_sa_diag')
