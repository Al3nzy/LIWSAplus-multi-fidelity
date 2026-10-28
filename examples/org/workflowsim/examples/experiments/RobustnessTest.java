package org.workflowsim.examples.experiments;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Random;
import org.workflowsim.CondorVM;
import org.workflowsim.Task;
import org.workflowsim.planning.FaithfulModel;
import org.workflowsim.planning.LIWSAPlanningAlgorithm;
import org.workflowsim.planning.LIWSASeededPlanningAlgorithm;

/**
 * Robustness to the idealised inputs of the planning model. Schedules are planned on nominal task lengths and peak VM speeds
 * (the model assumption criticised in the reviews) and then re-evaluated under perturbed inputs:
 *   len10 / len25 : every task length multiplied by an independent log-normal factor (mean 1, sigma 0.10 / 0.25),
 *   mips          : every VM delivers a random fraction U(0.6, 1.0) of its peak MIPS,
 *   both          : len25 and mips together.
 * The same scenarios are used for every algorithm, seed and workflow. Output: one JSON line per (workflow, algorithm, seed,
 * setting) with the non-dominated front of every scenario (see experiments/robustness_report.py).
 * Usage: RobustnessTest <out.jsonl> <instances|SMALL|LARGE|ALL> <numSeeds> <cond1,...> <numScenarios> [seedsDir]
 */
public class RobustnessTest {
    static class P extends LIWSASeededPlanningAlgorithm {
        void init() { setupProblem(); }
        List<Task> order() { return taskOrder; }
        List<CondorVM> vms() { return vmList; }
        List<int[]> seeds() { return generateSeedGenotypes(); }
    }

    static final String[] SETTINGS = {"len10", "len25", "mips", "both"};

    public static void main(String[] a) throws Exception {
        List<String> instances = new ArrayList<>();
        if (a[1].equals("SMALL") || a[1].equals("LARGE") || a[1].equals("ALL")) {
            for (String[] fam : ExperimentRunner.FAMILIES) {
                for (int s = 0; s < 4; s++) {
                    boolean large = s == 3;
                    if (a[1].equals("ALL") || (a[1].equals("SMALL") && !large) || (a[1].equals("LARGE") && large)) {
                        instances.add(fam[s]);
                    }
                }
            }
        } else {
            instances = java.util.Arrays.asList(a[1].split(","));
        }
        int seeds = Integer.parseInt(a[2]);
        String[] conds = a[3].split(",");
        int S = Integer.parseInt(a[4]);
        String seedsDir = a.length > 5 ? a[5] : "experiments/work/seeds";
        java.util.Set<String> done = new java.util.HashSet<>();
        java.io.File f = new java.io.File(a[0]);
        if (f.exists()) {
            try (java.io.BufferedReader r = new java.io.BufferedReader(new java.io.FileReader(f))) {
                String l;
                java.util.regex.Pattern pat = java.util.regex.Pattern.compile("\"instance\":\"([^\"]+)\",\"cond\":\"([^\"]+)\",\"seed\":(\\d+)");
                while ((l = r.readLine()) != null) {
                    java.util.regex.Matcher m = pat.matcher(l);
                    if (m.find()) {
                        done.add(m.group(1) + "|" + m.group(2) + "|" + m.group(3));
                    }
                }
            }
        }
        try (PrintWriter out = new PrintWriter(new java.io.FileWriter(f, true))) {
            for (String wf : instances) {
                List<Map<Integer, Integer>> warm = new ArrayList<>();
                warm.add(ExperimentRunner.readSeed(seedsDir + "/" + wf + "_HEFT.txt"));
                warm.add(ExperimentRunner.readSeed(seedsDir + "/" + wf + "_MINMIN.txt"));
                List<Task> tasks = StandaloneHarness.loadTasks("config/dax/" + wf + ".xml");
                LIWSAPlanningAlgorithm.CONFIG_SEED_ASSIGNMENTS = warm;
                LIWSASeededPlanningAlgorithm.CONFIG_SWEEP_SEEDS = 0;
                LIWSAPlanningAlgorithm.CONFIG_FAITHFUL = false;
                P p = new P();
                p.setTaskList(tasks);
                p.setVmList(StandaloneHarness.makeVms());
                p.init();
                FaithfulModel fm = new FaithfulModel(p.order(), p.vms());
                int n = p.order().size();
                int m = p.vms().size();
                // common scenarios for this workflow
                double[][][] lm = new double[SETTINGS.length][S][n];
                double[][][] mm = new double[SETTINGS.length][S][m];
                Random rs = new Random(97L * wf.hashCode() + 5);
                for (int st = 0; st < SETTINGS.length; st++) {
                    for (int s = 0; s < S; s++) {
                        double sigma = SETTINGS[st].equals("len10") ? 0.10 : (SETTINGS[st].equals("len25") || SETTINGS[st].equals("both")) ? 0.25 : 0.0;
                        for (int k = 0; k < n; k++) {
                            lm[st][s][k] = sigma == 0 ? 1.0 : Math.exp(sigma * rs.nextGaussian() - sigma * sigma / 2);
                        }
                        boolean degrade = SETTINGS[st].equals("mips") || SETTINGS[st].equals("both");
                        for (int v = 0; v < m; v++) {
                            mm[st][s][v] = degrade ? 0.6 + 0.4 * rs.nextDouble() : 1.0;
                        }
                    }
                }
                // heuristics as references
                List<int[]> sd = p.seeds();
                emit(out, wf, "HEFT", 0, java.util.Collections.singletonList(sd.get(0)), fm, lm, mm, S, done);
                emit(out, wf, "MINMIN", 0, java.util.Collections.singletonList(sd.get(1)), fm, lm, mm, S, done);
                for (String cond : conds) {
                    for (long seed = 1; seed <= seeds; seed++) {
                        if (done.contains(wf + "|" + cond + "|" + seed)) {
                            continue;
                        }
                        ExperimentRunner.runOne(wf, cond, seed, warm);
                        List<int[]> gen = new ArrayList<>(LIWSAPlanningAlgorithm.lastRun.archiveGenotypes);
                        emit(out, wf, cond, seed, gen, fm, lm, mm, S, done);
                    }
                    System.out.println("done " + wf + " " + cond);
                }
                LIWSAPlanningAlgorithm.CONFIG_SEED_ASSIGNMENTS = warm;
            }
        }
    }

    static void emit(PrintWriter out, String wf, String cond, long seed, List<int[]> gen, FaithfulModel fm,
                     double[][][] lm, double[][][] mm, int S, java.util.Set<String> done) {
        if (done.contains(wf + "|" + cond + "|" + seed)) {
            return;
        }
        StringBuilder sb = new StringBuilder();
        sb.append("{\"type\":\"rob\",\"instance\":\"").append(wf).append("\",\"cond\":\"").append(cond).append("\",\"seed\":")
          .append(seed).append(",\"nominal\":");
        List<double[]> nom = new ArrayList<>();
        for (int[] g : gen) {
            nom.add(fm.decode(g));
        }
        sb.append(front(nom)).append(",\"settings\":{");
        for (int st = 0; st < SETTINGS.length; st++) {
            if (st > 0) {
                sb.append(',');
            }
            sb.append('"').append(SETTINGS[st]).append("\":[");
            for (int s = 0; s < S; s++) {
                List<double[]> pts = new ArrayList<>();
                for (int[] g : gen) {
                    pts.add(fm.decodePerturbed(g, lm[st][s], mm[st][s]));
                }
                if (s > 0) {
                    sb.append(',');
                }
                sb.append(front(pts));
            }
            sb.append(']');
        }
        sb.append("}}");
        out.println(sb);
        out.flush();
    }

    static String front(List<double[]> pts) {
        List<double[]> nd = new ArrayList<>();
        for (double[] q : pts) {
            boolean dom = false;
            for (double[] r : pts) {
                if (r != q && r[0] <= q[0] && r[1] <= q[1] && (r[0] < q[0] || r[1] < q[1])) {
                    dom = true;
                    break;
                }
            }
            if (!dom) {
                boolean dup = false;
                for (double[] r : nd) {
                    if (r[0] == q[0] && r[1] == q[1]) {
                        dup = true;
                        break;
                    }
                }
                if (!dup) {
                    nd.add(q);
                }
            }
        }
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < nd.size(); i++) {
            if (i > 0) {
                sb.append(',');
            }
            sb.append('[').append(nd.get(i)[0]).append(',').append(nd.get(i)[1]).append(']');
        }
        return sb.append(']').toString();
    }
}
