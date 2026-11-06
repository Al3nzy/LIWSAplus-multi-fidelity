package org.workflowsim.examples.experiments;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.workflowsim.planning.LIWSAPlanningAlgorithm;

/**
 * Scenario-level validation inside WorkflowSim. For each (instance, algorithm, seed) the scenario-based search is run
 * as in the main study; T archive solutions (evenly spaced from fastest to cheapest by the search objective) are then
 * replayed in the simulator under R of the out-of-sample scenarios: task lengths are scaled by the scenario's length
 * multipliers and VM speeds by its speed multipliers inside WorkflowSim, and the simulator-measured (makespan, cost)
 * are written next to the decoder values for the same schedule and scenario.
 * Usage: ScenarioReplay <out.csv> <instances|SMALL|LARGE> <numSeeds> <cond1,cond2,...> <T> <R> [seedsDir]
 * JVM options are those of the corresponding study (e.g. -Dgens=20 -Drobust=30 -Dcvar=0.8 -Dslowp=0.2 -Dmipslow=0.8).
 */
public class ScenarioReplay {
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
        int T = Integer.parseInt(a[4]);
        int R = Integer.parseInt(a[5]);
        String seedsDir = a.length > 6 ? a[6] : "experiments/work/seeds";
        boolean exists = new java.io.File(a[0]).exists();
        java.util.Set<String> done = new java.util.HashSet<>();
        if (exists) {
            try (java.io.BufferedReader r = new java.io.BufferedReader(new java.io.FileReader(a[0]))) {
                String l = r.readLine();
                while ((l = r.readLine()) != null) {
                    String[] p = l.split(",");
                    done.add(p[0] + "|" + p[1] + "|" + p[2]);
                }
            }
        }
        try (PrintWriter out = new PrintWriter(new java.io.FileWriter(a[0], true))) {
            if (!exists) {
                out.println("instance,cond,seed,sched,scen,dec_makespan,dec_cost,sim_makespan,sim_cost,arch_makespan,arch_cost");
            }
            for (String wf : instances) {
                List<Map<Integer, Integer>> warm = new ArrayList<>();
                warm.add(ExperimentRunner.readSeed(seedsDir + "/" + wf + "_HEFT.txt"));
                warm.add(ExperimentRunner.readSeed(seedsDir + "/" + wf + "_MINMIN.txt"));
                for (String cond : conds) {
                    for (long seed = 1; seed <= seeds; seed++) {
                        if (done.contains(wf + "|" + cond + "|" + seed)) {
                            continue;
                        }
                        ExperimentRunner.runOne(wf, cond, seed, warm);
                        LIWSAPlanningAlgorithm.LastRunMetrics m = LIWSAPlanningAlgorithm.lastRun;
                        LIWSAPlanningAlgorithm alg = LIWSAPlanningAlgorithm.lastAlgorithm;
                        double[][][] sc = alg.replayScenarios(R);
                        List<double[]> pts = m.archiveFront;
                        List<int[]> gen = m.archiveGenotypes;
                        Integer[] order = new Integer[pts.size()];
                        for (int i = 0; i < order.length; i++) {
                            order[i] = i;
                        }
                        java.util.Arrays.sort(order, (x, y) -> Double.compare(pts.get(x)[0], pts.get(y)[0]));
                        int take = Math.min(T, order.length);
                        for (int t = 0; t < take; t++) {
                            int idx = order[take == 1 ? 0 : (int) Math.round((double) t * (order.length - 1) / (take - 1))];
                            int[] g = gen.get(idx);
                            for (int s = 0; s < R; s++) {
                                double[] dec = alg.decodeScenario(g, sc[0][s], sc[1][s]);
                                LIWSAPlanningAlgorithm.CONFIG_SIM_LENMULT = sc[0][s];
                                SeedBuilder.SIM_MIPSMULT = sc[1][s];
                                double[] sim;
                                try {
                                    sim = SeedBuilder.simulateForced(wf, g);
                                } finally {
                                    LIWSAPlanningAlgorithm.CONFIG_SIM_LENMULT = null;
                                    SeedBuilder.SIM_MIPSMULT = null;
                                }
                                out.println(wf + "," + cond + "," + seed + "," + t + "," + s + "," + dec[0] + "," + dec[1] + ","
                                        + sim[0] + "," + sim[1] + "," + pts.get(idx)[0] + "," + pts.get(idx)[1]);
                            }
                        }
                        out.flush();
                    }
                    System.out.println("done " + wf + " " + cond);
                }
            }
        }
    }
}
