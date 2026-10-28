package org.workflowsim.examples.experiments;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.workflowsim.CondorVM;
import org.workflowsim.Task;
import org.workflowsim.planning.LIWSAPlanningAlgorithm;
import org.workflowsim.planning.LIWSAPlus;
import org.workflowsim.planning.LIWSASeededPlanningAlgorithm;
import org.workflowsim.planning.MOEAHosts;

/**
 * Experiment driver. One JSON line per (instance, condition, seed).
 *
 * Condition grammar:  HOST[.flag]*-START
 *   HOST   : LIWSA | NSGA2 | SPEA2 | MOEAD | LIWSAP (LIWSA+)
 *   START  : plus (HEFT, Min-Min and weighted list-scheduling seeds injected) | warm (HEFT and Min-Min only) | cold (random only)
 *   flags  : f (simulator-faithful decoder), r (expensive Monte-Carlo objective), san | sao | sar (pre-screening: nominal,
 *            online ridge correction, random), hvi (hypervolume-improvement rule), imm (OLS immigrants),
 *            nog | nol | noc | noa | noe (ablations of LIWSA+)
 *
 * Usage: ExperimentRunner <out.jsonl> <instances|ALL|SMALL|LARGE> <numSeeds> <cond1,cond2,...> [seedsDir]
 */
public class ExperimentRunner {

    static final String[][] FAMILIES = {
        {"Montage_25", "Montage_50", "Montage_100", "Montage_1000"},
        {"CyberShake_30", "CyberShake_50", "CyberShake_100", "CyberShake_1000"},
        {"Sipht_30", "Sipht_60", "Sipht_100", "Sipht_1000"},
        {"Epigenomics_24", "Epigenomics_46", "Epigenomics_100", "Epigenomics_997"},
        {"Inspiral_30", "Inspiral_50", "Inspiral_100", "Inspiral_1000"}};

    public static void main(String[] a) throws Exception {
        String outPath = a[0];
        List<String> instances = new ArrayList<>();
        for (String[] fam : FAMILIES) {
            for (int s = 0; s < 4; s++) {
                boolean large = (s == 3);
                if (a[1].equals("ALL") || (a[1].equals("SMALL") && !large) || (a[1].equals("LARGE") && large)) {
                    instances.add(fam[s]);
                }
            }
        }
        if (!a[1].equals("ALL") && !a[1].equals("SMALL") && !a[1].equals("LARGE")) {
            instances = java.util.Arrays.asList(a[1].split(","));
        }
        int numSeeds = Integer.parseInt(a[2]);
        String[] conds = a[3].split(",");
        String seedsDir = a.length > 4 ? a[4] : "experiments/work/seeds";
        java.util.Set<String> done = new java.util.HashSet<>();
        java.util.Set<String> refDone = new java.util.HashSet<>();
        if (new File(outPath).exists()) {
            try (BufferedReader r = new BufferedReader(new FileReader(outPath))) {
                String l;
                java.util.regex.Pattern pat = java.util.regex.Pattern.compile("\"instance\":\"([^\"]+)\",\"cond\":\"([^\"]+)\",\"seed\":(\\d+)");
                while ((l = r.readLine()) != null) {
                    java.util.regex.Matcher mm = pat.matcher(l);
                    if (mm.find()) {
                        done.add(mm.group(1) + "|" + mm.group(2) + "|" + mm.group(3));
                    } else if (l.startsWith("{\"type\":\"ref\"")) {
                        int i0 = l.indexOf("\"instance\":\"") + 12;
                        refDone.add(l.substring(i0, l.indexOf('"', i0)));
                    }
                }
            }
        }
        try (PrintWriter out = new PrintWriter(new java.io.FileWriter(outPath, true))) {
            for (String wf : instances) {
                List<Map<Integer, Integer>> warm = new ArrayList<>();
                warm.add(readSeed(seedsDir + "/" + wf + "_HEFT.txt"));
                warm.add(readSeed(seedsDir + "/" + wf + "_MINMIN.txt"));
                // reference record: decoder-level objective values of the two heuristic schedules
                if (!refDone.contains(wf)) {
                    out.println(refRecord(wf, warm));
                    out.flush();
                }
                for (String cond : conds) {
                    long seed0 = Long.getLong("seed0", 0L);
                    for (long seed = seed0 + 1; seed <= seed0 + numSeeds; seed++) {
                        if (done.contains(wf + "|" + cond + "|" + seed)) {
                            continue;
                        }
                        try {
                            out.println(runOne(wf, cond, seed, warm));
                        } catch (Exception e) {
                            System.err.println("FAILED " + wf + " " + cond + " " + seed + ": " + e);
                            e.printStackTrace();
                        }
                        out.flush();
                    }
                    System.out.println("done " + wf + " " + cond);
                }
            }
        }
    }

    static Map<Integer, Integer> readSeed(String path) throws Exception {
        Map<Integer, Integer> m = new HashMap<>();
        try (BufferedReader r = new BufferedReader(new FileReader(path))) {
            String l;
            while ((l = r.readLine()) != null) {
                String[] p = l.trim().split("\\s+");
                if (p.length == 2) {
                    m.put(Integer.parseInt(p[0]), Integer.parseInt(p[1]));
                }
            }
        }
        return m;
    }

    static class Probe extends LIWSASeededPlanningAlgorithm {
        boolean skipCheck = false;

        List<double[]> seedObjectives() {
            setupProblem();
            List<double[]> out = new ArrayList<>();
            for (int[] g : generateSeedGenotypes()) {
                double[] mc = decode(g);
                double[] ref = skipCheck ? mc : decodeReference(g);
                if (mc[0] != ref[0] || mc[1] != ref[1]) {
                    throw new IllegalStateException("fast decoder mismatch");
                }
                out.add(mc);
            }
            return out;
        }
    }

    static String refRecord(String wf, List<Map<Integer, Integer>> warm) {
        List<Task> tasks = StandaloneHarness.loadTasks("config/dax/" + wf + ".xml");
        List<CondorVM> vms = StandaloneHarness.makeVms();
        LIWSAPlanningAlgorithm.CONFIG_RANDOM_SEED = 1L;
        LIWSAPlanningAlgorithm.CONFIG_SEED_ASSIGNMENTS = warm;
        LIWSAPlanningAlgorithm.CONFIG_FAITHFUL = false;
        Probe p = new Probe();
        p.setTaskList(tasks);
        p.setVmList(vms);
        List<double[]> pts = p.seedObjectives();
        LIWSAPlanningAlgorithm.CONFIG_FAITHFUL = true;
        Probe pf = new Probe();
        pf.skipCheck = true;
        pf.setTaskList(StandaloneHarness.loadTasks("config/dax/" + wf + ".xml"));
        pf.setVmList(StandaloneHarness.makeVms());
        List<double[]> ptsF = pf.seedObjectives();
        LIWSAPlanningAlgorithm.CONFIG_FAITHFUL = false;
        return "{\"type\":\"ref\",\"instance\":\"" + wf + "\",\"n_tasks\":" + tasks.size()
                + ",\"heft\":" + pt(pts.get(0)) + ",\"minmin\":" + pt(pts.get(1))
                + ",\"heft_f\":" + pt(ptsF.get(0)) + ",\"minmin_f\":" + pt(ptsF.get(1)) + "}";
    }

    static String pt(double[] p) {
        return "[" + p[0] + "," + p[1] + "]";
    }

    static String pts(List<double[]> l) {
        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < l.size(); i++) {
            if (i > 0) {
                sb.append(',');
            }
            sb.append(pt(l.get(i)));
        }
        return sb.append(']').toString();
    }

    static String runOne(String wf, String cond, long seed, List<Map<Integer, Integer>> warm) throws Exception {
        String[] hs = cond.split("\\+");
        String[] hostStart = hs[0].split("-");
        String host = hostStart[0];
        boolean isWarm = !hostStart[1].equals("cold");
        boolean plus = hostStart[1].equals("plus");     // warm + weighted list-scheduling sweep seeds
        String[] hostFlags = host.split("\\.");
        host = hostFlags[0];
        LIWSAPlus.CONFIG_POOL = 31;
        LIWSAPlus.CONFIG_AOS = true;
        LIWSAPlanningAlgorithm.CONFIG_FAITHFUL = false;
        LIWSAPlanningAlgorithm.CONFIG_ROBUST_S = 0;
        LIWSAPlanningAlgorithm.CONFIG_ROBUST_CVAR = Double.parseDouble(System.getProperty("cvar", "0"));
        LIWSAPlanningAlgorithm.CONFIG_ROBUST_SLOWP = Double.parseDouble(System.getProperty("slowp", "0"));
        LIWSAPlanningAlgorithm.CONFIG_ROBUST_MIPS_LOW = Double.parseDouble(System.getProperty("mipslow", "0.6"));
        LIWSAPlanningAlgorithm.CONFIG_ROBUST_SIGMA = Double.parseDouble(System.getProperty("sigma", "0.25"));
        LIWSAPlus.CONFIG_SA_MODE = 0;
        LIWSAPlus.CONFIG_IMM = false;
        LIWSAPlus.CONFIG_SA_HVI = false;
        LIWSAPlus.CONFIG_SA_MIX = false;
        LIWSAPlus.CONFIG_DIAG = Boolean.getBoolean("diag");
        LIWSAPlus.lastDiag = null;
        LIWSAPlus.CONFIG_SA_K = Integer.getInteger("sak", 5);
        LIWSAPlus.CONFIG_SA_RIDGE = Double.parseDouble(System.getProperty("ridge", "0.1"));
        LIWSAPlus.CONFIG_SA_RIDGE_RICH = Double.parseDouble(System.getProperty("ridgerich", "1.0"));
        LIWSAPlus.CONFIG_GBM_TREES = Integer.getInteger("gbmtrees", 60);
        LIWSAPlus.CONFIG_GBM_DEPTH = Integer.getInteger("gbmdepth", 3);
        LIWSAPlus.CONFIG_IMM_COUNT = Integer.getInteger("imm", 4);
        LIWSAPlus.CONFIG_IMM_PERIOD = Integer.getInteger("immperiod", 10);
        for (int i = 1; i < hostFlags.length; i++) {
            switch (hostFlags[i]) {
                case "nog": LIWSAPlus.CONFIG_POOL &= ~(LIWSAPlus.OP_CRIT | LIWSAPlus.OP_COST | LIWSAPlus.OP_EST); break;
                case "nol": LIWSAPlus.CONFIG_POOL &= ~LIWSAPlus.OP_LIWSA; break;
                case "noc": LIWSAPlus.CONFIG_POOL &= ~LIWSAPlus.OP_CANON; break;
                case "noa": LIWSAPlus.CONFIG_AOS = false; break;
                case "noe": LIWSAPlus.CONFIG_POOL &= ~LIWSAPlus.OP_EST; break;
                case "f": LIWSAPlanningAlgorithm.CONFIG_FAITHFUL = true; break;
                case "r": LIWSAPlanningAlgorithm.CONFIG_ROBUST_S = Integer.getInteger("robust", 20); break;
                case "san": LIWSAPlus.CONFIG_SA_MODE = 1; break;
                case "sao": LIWSAPlus.CONFIG_SA_MODE = 2; break;
                case "sar": LIWSAPlus.CONFIG_SA_MODE = 3; break;
                case "sal": LIWSAPlus.CONFIG_SA_MODE = 4; break;
                case "sag": LIWSAPlus.CONFIG_SA_MODE = 5; break;
                case "mix": LIWSAPlus.CONFIG_SA_MIX = true; break;
                case "imm": LIWSAPlus.CONFIG_IMM = true; break;
                case "hvi": LIWSAPlus.CONFIG_SA_HVI = true; break;
                default: throw new IllegalArgumentException("unknown flag " + hostFlags[i]);
            }
        }
        LIWSASeededPlanningAlgorithm.CONFIG_SWEEP_SEEDS = plus ? Integer.getInteger("sweep", 8) : 0;
        LIWSAPlus.CONFIG_MOVE_FRAC = Double.parseDouble(System.getProperty("movefrac", "0.02"));
        LIWSAPlus.CONFIG_AOS_PMIN = Double.parseDouble(System.getProperty("pmin", "0.10"));
        LIWSAPlus.CONFIG_AOS_ALPHA = Double.parseDouble(System.getProperty("alpha", "0.30"));

        List<Task> tasks = StandaloneHarness.loadTasks("config/dax/" + wf + ".xml");
        List<CondorVM> vms = StandaloneHarness.makeVms();

        int pop = Integer.getInteger("pop", 30);
        LIWSAPlanningAlgorithm.CONFIG_POPULATION_SIZE = pop;
        // the evaluation budget stays at 2930 whatever the population size: gens = (2930 - pop) / (pop - 1)
        LIWSAPlanningAlgorithm.CONFIG_GENERATION_COUNT = Integer.getInteger("gens", (2930 - pop) / (pop - 1));
        LIWSAPlanningAlgorithm.CONFIG_RANDOM_SEED = seed;
        LIWSAPlanningAlgorithm.CONFIG_SEED_ASSIGNMENTS = isWarm ? warm : null;
        LIWSAPlanningAlgorithm.CONFIG_TRACE = true;

        long searchMs;
        List<double[]> front;
        List<double[]> archive = new ArrayList<>();
        long evals = -1;
        List<double[]> oos = null;
        long nomDec = 0;
        String diag = null;
        List<LIWSAPlanningAlgorithm.Checkpoint> trace = new ArrayList<>();
        long wall0 = System.nanoTime();
        LIWSASeededPlanningAlgorithm alg;
        switch (host) {
            case "LIWSA":
                alg = new LIWSASeededPlanningAlgorithm();
                break;
            case "NSGA2":
                alg = new MOEAHosts.NSGA2();
                break;
            case "SPEA2":
                alg = new MOEAHosts.SPEA2();
                break;
            case "MOEAD":
                alg = new MOEAHosts.MOEAD();
                break;
            case "LIWSAP":
                alg = new LIWSAPlus();
                break;
            default:
                throw new IllegalArgumentException(host);
        }
        alg.setTaskList(tasks);
        alg.setVmList(vms);
        alg.run();
        LIWSAPlanningAlgorithm.LastRunMetrics m = LIWSAPlanningAlgorithm.lastRun;
        front = m.paretoFrontPoints;
        archive = m.archiveFront;
        searchMs = m.searchWallClockMillis;
        if (LIWSAPlanningAlgorithm.CONFIG_ROBUST_S > 0) {
            oos = alg.outOfSampleFront(m.archiveGenotypes, 100);
        }
        nomDec = alg.nominalDecodeCount();
        if (LIWSAPlus.lastDiag != null) {
            diag = java.util.Arrays.toString(LIWSAPlus.lastDiag);
        }
        evals = m.evaluations;
        trace = m.trace;
        long wallNs = System.nanoTime() - wall0;

        StringBuilder sb = new StringBuilder();
        sb.append("{\"type\":\"run\",\"instance\":\"").append(wf).append("\",\"cond\":\"").append(cond)
          .append("\",\"seed\":").append(seed).append(",\"evals\":").append(evals)
          .append(",\"search_ms\":").append(searchMs).append(",\"wall_ms\":").append(wallNs / 1000000.0);
        if (oos != null) {
            sb.append(",\"oos\":").append(pts(oos)).append(",\"nom_dec\":").append(nomDec);
            if (diag != null) {
                sb.append(",\"diag\":").append(diag);
            }
        }
        sb.append(",\"front\":").append(pts(front)).append(",\"archive\":").append(pts(archive)).append(",\"trace\":[");
        for (int i = 0; i < trace.size(); i++) {
            LIWSAPlanningAlgorithm.Checkpoint c = trace.get(i);
            if (i > 0) {
                sb.append(',');
            }
            sb.append("{\"g\":").append(c.generation).append(",\"e\":").append(c.evaluations)
              .append(",\"ms\":").append(c.nanos / 1e6).append(",\"f\":").append(pts(c.front)).append('}');
        }
        sb.append("]}");
        return sb.toString();
    }
}
