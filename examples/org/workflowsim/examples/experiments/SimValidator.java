package org.workflowsim.examples.experiments;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.workflowsim.planning.LIWSAPlanningAlgorithm;

/**
 * End-to-end check inside the WorkflowSim simulator. For each (instance, algorithm, seed) the algorithm is run with the
 * planning-level decoder, up to K representative solutions of its archive front (evenly spaced from fastest to cheapest)
 * are replayed in the full simulation, and the simulator-measured (makespan, cost) are written out so that
 * simulator-level hypervolumes can be compared (experiments/sim_report.py).
 * Usage: SimValidator <out.csv> <instances|SMALL|LARGE> <numSeeds> <cond1,cond2,...> <K> [seedsDir]
 */
public class SimValidator {
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
        int K = Integer.parseInt(a[4]);
        String seedsDir = a.length > 5 ? a[5] : "experiments/work/seeds";
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
                out.println("instance,cond,seed,point,dec_makespan,dec_cost,sim_makespan,sim_cost");
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
                        List<double[]> pts = m.archiveFront;
                        List<int[]> gen = m.archiveGenotypes;
                        Integer[] order = new Integer[pts.size()];
                        for (int i = 0; i < order.length; i++) {
                            order[i] = i;
                        }
                        java.util.Arrays.sort(order, (x, y) -> Double.compare(pts.get(x)[0], pts.get(y)[0]));
                        int take = Math.min(K, order.length);
                        for (int t = 0; t < take; t++) {
                            int idx = order[take == 1 ? 0 : (int) Math.round((double) t * (order.length - 1) / (take - 1))];
                            double[] sim = SeedBuilder.simulateForced(wf, gen.get(idx));
                            out.println(wf + "," + cond + "," + seed + "," + t + "," + pts.get(idx)[0] + "," + pts.get(idx)[1]
                                    + "," + sim[0] + "," + sim[1]);
                        }
                        out.flush();
                    }
                    System.out.println("done " + wf + " " + cond);
                }
            }
        }
    }
}
