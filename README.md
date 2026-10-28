# Multi-fidelity risk-aware workflow scheduling: code and simulations

Everything needed to run the algorithms of the study inside WorkflowSim and to regenerate every number, table and figure.
Nothing else from the research repository is included (no GNN warm-start code, no earlier baselines).

## Contents
| path | what |
|---|---|
| `sources/` | CloudSim 3 and WorkflowSim 1.1.0 sources reduced to the classes the experiments use, plus `org/workflowsim/planning/`: `LIWSAPlanningAlgorithm` (original LIWSA and the shared decoder, objectives and scenario model), `FaithfulModel` (simulator-faithful decoder), `LIWSAPlus` (guided evolutionary search, pre-screening, learned surrogates), `OnlineRidge`, `BoostedTrees`, `MOEAHosts` (NSGA-II, SPEA2, MOEA/D), HEFT |
| `examples/.../experiments/` | `ExperimentRunner` (one JSON line per instance, condition and seed), `SeedBuilder`, `FidelityTest`, `RobustnessTest`, `SimValidator` |
| `config/dax/` | the 20 Pegasus workflows (Montage, CyberShake, Sipht, Epigenomics, Inspiral; four sizes each) |
| `experiments/` | `run_experiments.py` (driver), `paper_statistics.py`, `make_paper_figures.py`, `ml_checks.py`, analysis helpers |
| `experiments/reference_results/` | output of a complete run: `paper_statistics.txt`, `ml_checks.txt`, `report_sim.txt`, 25 LaTeX tables, figures |
| `run` | one-command entry point (Code Ocean or shell); `QUICK=1 ./run` is a one-minute smoke test |
| `MANIFEST.txt`, `experiments/check_install.py` | checksums of every file and a checker: run `python experiments/check_install.py` after extracting to confirm the folder is an unmodified copy |

## One command (recommended)
    python run_all.py            # runs every step (compile seeds fidelity faithful sensitivity robustness multifidelity simcheck mlchecks figures), then collects everything into RESULTS/
    python run_all.py --quick    # smoke test of the whole chain
    run_all.bat --quick          # Windows shortcut
`RESULTS/index.html` shows all figures and links to every table and report; `RESULTS/reports/` holds the numbers, `RESULTS/tables/` the LaTeX tables,
`RESULTS/figures/` the figures (png and pdf), `RESULTS/run_log.txt` the log. A failing step does not stop the others. A full run is resumable.

## Requirements
JDK 8 or later (set `JAVAC` if the compiler is not called `javac`), Python 3 with numpy, scipy and matplotlib. No GPU, no PyTorch.

## Run
    python experiments/run_experiments.py paper            # everything, resumable, about 1.5 h on one core
    python experiments/run_experiments.py paper --quick    # 2 workflows, 2 seeds
Steps (each can be run alone): `compile seeds fidelity faithful sensitivity robustness multifidelity simcheck mlchecks figures`.
Generated files go to `experiments/work/` (`results/*.jsonl`, `analysis/`). Finished work is skipped, so an interrupted run can be restarted.

## Conditions (`HOST[.flag]*-START`)
Hosts: `LIWSA` (original), `LIWSAP` (LIWSA+, the guided search), `NSGA2`, `SPEA2`, `MOEAD`. Start: `plus` (HEFT, Min-Min and 8 weighted list-scheduling seeds), `warm` (HEFT, Min-Min), `cold`.
Flags: `f` simulator-faithful decoder; `r` expensive Monte-Carlo objective (`-Drobust=S`, `-Dcvar=a`, `-Dslowp`, `-Dmipslow`); pre-screening `san` nominal, `sao` online ridge on 10 descriptors, `sal` ridge on rich features, `sag` boosted trees on rich features, `sar` random; `hvi` hypervolume-improvement rule; `mix` hedged (alternating nominal and learned) selection; `imm` OLS immigrants; `nog nol noc noa noe` ablations. Pre-screening factor `-Dsak=K` (default 5).

## Verification performed
* The reduced sources reproduce the full research code bit for bit on 52 run records (nominal and risk-aware, all pre-screening variants); the added surrogates leave all earlier conditions unchanged (28 records).
* Searches and HEFT/Min-Min seeds do not depend on the JVM identity hash (tested with three hash modes).
* A complete run from a fresh unzip is in `experiments/reference_results`. Hypervolume reference points are computed from all algorithms present in a result file, so adding a condition changes percentage gains slightly (for example the pre-screening gain over LIWSA+ moved from +2.36% to +2.25% when four conditions were added); always quote numbers from one complete run.

## Known limits
* `mlchecks` was run in full once from this package; the K sweep, harsher tail risk and hedged selection use the five development workflows only.
* The learned surrogates improve ranking quality but not the end result over nominal pre-screening; see `reference_results/ml_checks.txt`.
* Timings were measured on one core of the build machine; absolute times will differ.
