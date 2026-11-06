# Code changes made for the revised manuscript (copy into the repository, same relative paths)
- sources/org/workflowsim/planning/LIWSAPlanningAlgorithm.java: decodeIdeal() extracted (analytic gap-filling decoder, used as the idealised screening model); CONFIG_SCEN_SALT offsets the scenario seeds of the search and of the out-of-sample set (0 = published runs).
- sources/org/workflowsim/planning/LIWSAPlus.java: CONFIG_SA_IDEAL (screen with the idealised decoder), CONFIG_SA_QUOTA (exploration quota); selectCandidates() wraps selectCandidatesCore().
- examples/org/workflowsim/examples/experiments/ExperimentRunner.java: condition tokens `sai` (idealised screening) and `saq` (quota 0.2); -Dscensalt=k.
- experiments/: run_ideal.py, run_rep.py (new runs); extra_analysis.py, ideal_analysis.py, rep_analysis.py, make_new_tables.py, fig_gain.py (analyses, new tables, Fig. 2); fig_framework.py (Fig. 1).
- Defaults are unchanged: with no new flag every published result is reproduced.
- The Code Ocean capsule (10.24433/CO.3382903.v1) does not contain these changes yet; update it before submission or the new tables and Fig. 2 cannot be regenerated from it.
