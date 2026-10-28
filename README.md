<div align="center">

# LIWSA-MF2

### Multi-fidelity, risk-aware scheduling for scientific workflows

Reproducible algorithms, WorkflowSim experiments, and analysis tools for studying
the trade-offs between workflow makespan, execution cost, and scheduling risk.

[![Java](https://img.shields.io/badge/Java-8%2B-007396?logo=openjdk&logoColor=white)](https://openjdk.org/)
[![Python](https://img.shields.io/badge/Python-3-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![WorkflowSim](https://img.shields.io/badge/Simulator-WorkflowSim-5865F2)](https://github.com/WorkflowSim/WorkflowSim-1.0)
[![Reproducible experiments](https://img.shields.io/badge/Research-reproducible-2E8B57)](#reproducing-and-interpreting-results)

</div>

---

## At a glance

LIWSA-MF2 is a research codebase for evaluating workflow schedulers under
multi-objective and risk-aware settings. It includes the original LIWSA method,
the guided **LIWSA+** search, simulator-faithful scheduling, competing
multi-objective algorithms, benchmark workflows, and scripts to regenerate
reports, tables, and plots.

The package contains **20 Pegasus workflow definitions** from five workflow
families, along with a complete set of reference experiment outputs.

## Results

The figures below are generated from the included experiment outputs. Click a
figure to open its PDF version.

<div align="center">
  <a href="RESULTS/figures/fig_framework.pdf">
    <img src="RESULTS/figures/fig_framework.png" alt="Overview of the LIWSA-MF2 scheduling framework" width="90%">
  </a>
  <br>
  <em>Overview of the scheduling and evaluation framework.</em>
</div>

<table>
  <tr>
    <td align="center" width="50%">
      <a href="RESULTS/figures/fig_main_diff.pdf">
        <img src="RESULTS/figures/fig_main_diff.png" alt="Main algorithm comparison" width="100%">
      </a>
      <br><em>Main algorithm comparison</em>
    </td>
    <td align="center" width="50%">
      <a href="RESULTS/figures/fig_robust.pdf">
        <img src="RESULTS/figures/fig_robust.png" alt="Risk-aware and robustness results" width="100%">
      </a>
      <br><em>Risk-aware and robustness results</em>
    </td>
  </tr>
  <tr>
    <td align="center" width="50%">
      <a href="RESULTS/figures/fig_fidelity.pdf">
        <img src="RESULTS/figures/fig_fidelity.png" alt="Decoder fidelity results" width="100%">
      </a>
      <br><em>Decoder fidelity assessment</em>
    </td>
    <td align="center" width="50%">
      <a href="RESULTS/figures/fig_convergence.pdf">
        <img src="RESULTS/figures/fig_convergence.png" alt="Search convergence results" width="100%">
      </a>
      <br><em>Search convergence across workflow sizes</em>
    </td>
  </tr>
</table>

In the simulator-replay analysis, LIWSA+ recorded a **+1.84% mean hypervolume
difference** versus original LIWSA across 15 instances (225 runs). This
comparison was **not statistically significant** in the reported paired test
(*p* = 0.1354); see the complete report before interpreting the result.

Explore all available plots, tables, and reports in
[`RESULTS/index.html`](RESULTS/index.html) or browse the
[`RESULTS/` folder](RESULTS/).

## What’s included

| Area | Contents |
|---|---|
| **Scheduling algorithms** | Original LIWSA, LIWSA+, simulator-faithful decoder, NSGA-II, SPEA2, MOEA/D, and HEFT |
| **Risk and search variants** | Monte Carlo risk objectives, CVaR, candidate pre-screening, learned surrogates, and ablations |
| **Benchmarks** | 20 Pegasus DAX workflows from Montage, CyberShake, Sipht, Epigenomics, and Inspiral |
| **Experiment tooling** | Java experiment runners and Python scripts for running studies, analysis, statistics, and figures |
| **Reference outputs** | Figures, LaTeX tables, statistical reports, and simulation checks |

The reduced `sources/` tree contains the CloudSim/WorkflowSim classes used by
the experiments and the LIWSA planning implementations. The repository does
not include the earlier GNN warm-start code or earlier baselines.

## Quick start

### Requirements

- JDK 8 or later (`javac` must be available; set `JAVAC` if it has another name)
- Python 3 with `numpy`, `scipy`, and `matplotlib`
- No GPU or PyTorch required

### Run the quick end-to-end check

From the repository root:

```bash
python run_all.py --quick
```

On Windows, you can also run:

```bat
run_all.bat --quick
```

### Run the complete workflow

```bash
python run_all.py
```

To run the experiment driver directly:

```bash
python experiments/run_experiments.py paper
```

Add `--quick` to the direct experiment command for a smaller run. The complete
paper run takes about 1.5 hours on one core; runtime varies by machine. The
driver is resumable and skips completed work.

The pipeline runs compilation, seed generation, fidelity and robustness
experiments, analysis, and figure generation. Outputs are collected in
`RESULTS/`; intermediate experiment data is written to
`experiments/work/`. Individual stages can also be run directly:
`compile`, `seeds`, `fidelity`, `faithful`, `sensitivity`, `robustness`,
`multifidelity`, `simcheck`, `mlchecks`, and `figures`.

## Reproducing and interpreting results

- The included `experiments/reference_results/` contains a complete reference
  run with reports, tables, and figures.
- `RESULTS/index.html` is a convenient local index for browsing the collected
  output.
- Hypervolume reference points are computed from the algorithms present in
  each result file. Adding or removing conditions can therefore change
  percentage gains; compare results from the same complete run.
- Reported timings were measured on one core of the build machine and will vary
  across systems.
- The full `mlchecks` run was completed once for this package. The candidate
  pool-size sweep, harsher tail-risk setting, and hedged selection are limited
  to five development workflows; consult `experiments/reference_results/ml_checks.txt`
  for details.
- Learned surrogates improved ranking quality in the reported checks, but did
  not improve the final result over nominal pre-screening.

## Repository map

```text
config/dax/                       Pegasus workflow inputs
experiments/                      Experiment driver and analysis scripts
experiments/reference_results/    Complete reference reports, tables, figures
sources/                          Reduced simulator sources and scheduling code
RESULTS/                          Collected figures, reports, and LaTeX tables
run_all.py                        Cross-platform pipeline entry point
run_all.bat                       Windows pipeline shortcut
MANIFEST.txt                      File checksums
```

To verify a copy against the included manifest, run:

```bash
python experiments/check_install.py
```

## License and third-party notices

See [`license.txt`](license.txt) and
[`license-workflowsim.txt`](license-workflowsim.txt) for the licensing
information and third-party notices included with the project.
