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
 * Compares three evaluations of the same task-to-VM assignments: the analytic decoder, the simulator-faithful decoder
 * (FaithfulModel) and the full WorkflowSim simulation. Assignments: the HEFT and Min-Min schedules, random assignments and
 * random perturbations of HEFT. Usage: FidelityTest <out.csv> <instances> <numRandom> [seedsDir]
 */
public class FidelityTest {
    static class P extends LIWSASeededPlanningAlgorithm {
        void init() { setupProblem(); }
        List<Task> order() { return taskOrder; }
        double[] dec(int[] g) { return decode(g); }
        List<int[]> seeds() { return generateSeedGenotypes(); }
        List<CondorVM> vms() { return vmList; }
    }

    public static void main(String[] a) throws Exception {
        String seedsDir = a.length > 3 ? a[3] : "experiments/work/seeds";
        int nRand = Integer.parseInt(a[2]);
        try (PrintWriter out = new PrintWriter(a[0])) {
            out.println("instance,kind,sim_makespan,sim_cost,dec_makespan,dec_cost,faith_makespan,faith_cost");
            for (String wf : a[1].split(",")) {
                List<Map<Integer, Integer>> warm = new ArrayList<>();
                warm.add(ExperimentRunner.readSeed(seedsDir + "/" + wf + "_HEFT.txt"));
                warm.add(ExperimentRunner.readSeed(seedsDir + "/" + wf + "_MINMIN.txt"));
                List<Task> tasks = StandaloneHarness.loadTasks("config/dax/" + wf + ".xml");
                List<CondorVM> vms = StandaloneHarness.makeVms();
                LIWSAPlanningAlgorithm.CONFIG_SEED_ASSIGNMENTS = warm;
                LIWSASeededPlanningAlgorithm.CONFIG_SWEEP_SEEDS = 4;
                P p = new P();
                p.setTaskList(tasks);
                p.setVmList(vms);
                p.init();
                FaithfulModel fm = new FaithfulModel(p.order(), p.vms());
                List<int[]> gs = new ArrayList<>();
                List<String> kinds = new ArrayList<>();
                List<int[]> sd = p.seeds();
                String[] sk = {"heft", "minmin", "sweep1", "sweep2", "sweep3", "sweep4"};
                for (int i = 0; i < sd.size() && i < sk.length; i++) {
                    gs.add(sd.get(i));
                    kinds.add(sk[i]);
                }
                Random r = new Random(7);
                int n = p.order().size();
                for (int i = 0; i < nRand; i++) {
                    int[] g = new int[n];
                    if (i % 2 == 0) {
                        for (int k = 0; k < n; k++) {
                            g[k] = r.nextInt(16);
                        }
                        kinds.add("random");
                    } else {
                        g = sd.get(0).clone();
                        for (int k = 0; k < n; k++) {
                            if (r.nextDouble() < 0.1) {
                                g[k] = r.nextInt(16);
                            }
                        }
                        kinds.add("heft+10%");
                    }
                    gs.add(g);
                }
                for (int i = 0; i < gs.size(); i++) {
                    double[] d = p.dec(gs.get(i));
                    double[] f = fm.decode(gs.get(i));
                    double[] s = SeedBuilder.simulateForced(wf, gs.get(i));
                    out.println(wf + "," + kinds.get(i) + "," + s[0] + "," + s[1] + "," + d[0] + "," + d[1] + "," + f[0] + "," + f[1]);
                }
                out.flush();
                System.out.println("done " + wf);
            }
        }
    }
}
