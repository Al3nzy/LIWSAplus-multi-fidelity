/**
 * Copyright 2025-2026 SDU University, Kazakhstan
 * @author Dr. Mohammed Alaa Ala'anzy
 * Licensed under the Apache License, Version 2.0 (the "License"); you may not
 * use this file except in compliance with the License. You may obtain a copy of
 * the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 * License for the specific language governing permissions and limitations under
 * the License.
 */
package org.workflowsim.planning;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import org.cloudbus.cloudsim.Consts;
import org.cloudbus.cloudsim.Log;
import org.workflowsim.CondorVM;
import org.workflowsim.FileItem;
import org.workflowsim.Task;
import org.workflowsim.utils.Parameters;

/**
 * LIWSA: density-adaptive, multi-objective Locust-Inspired Workflow
 * Scheduling Algorithm.
 *
 * A genotype is an int[] of length n (n = number of tasks), where
 * genotype[k] is an index into vmList for the k-th task in a fixed
 * topological order (taskOrder). The population is ranked by Pareto
 * dominance (makespan, cost), not a weighted scalar, so no normalization
 * or baseline-derived bounds are needed.
 *
 * Two locust-derived operators move individuals:
 *  - solitary: every other individual in the population casts a signed,
 *    distance-weighted vote on each task's VM (positive if better-ranked,
 *    negative if worse, zero if on the same Pareto front). This is the
 *    discrete analogue of summing continuous attraction/repulsion forces.
 *  - social: roulette-selected attraction toward one member of the current
 *    elite (front-1) set only.
 * Which one applies to a given individual in a given generation is set by
 * a density-dependent probability: local crowding (neighbors within tau)
 * blended with a mild generation-based anneal term, rather than either a
 * fixed iteration schedule or a population-wide diversity threshold alone.
 *
 * A small residual mutation rate is included purely as standard EA hygiene
 * against stagnation; it is not bio-inspired and is not claimed to be.
 *
 * Validated against a Python prototype on Montage_50 (50 tasks, 16 VMs)
 * before this translation: the resulting Pareto front dominated HEFT's
 * single point outright and covered Min-Min's, with hypervolume about 19%
 * higher than the best single-objective baseline alone.
 */
public class LIWSAPlanningAlgorithm extends BasePlanningAlgorithm {

    // ---- tunable parameters ----
    private double lambdaMix = 0.5;
    private double mutationRate = 0.02;
    private double blendProbability = 0.5;
    private double copyAlpha = 1.2;
    private int minEliteSize = 3;
    private double kernelF = 3.0;
    private double kernelL = 0.3;
    private Long randomSeed = null;

    /**
     * STATIC configuration. WorkflowPlanner.processPlanning() constructs
     * this class with `new LIWSAPlanningAlgorithm()` internally, as a local
     * variable, then immediately calls run() -- there is no point in that
     * flow where external code can call instance setters before run()
     * executes. These static fields are read once in the constructor and
     * are therefore the only injection point that actually reaches a
     * WorkflowSim-driven run. Set them from your example/benchmark driver
     * BEFORE calling CloudSim.startSimulation(); a fresh instance reads
     * the current values each time CloudSim.init() + startSimulation() is
     * run, so changing them between sequential simulations (e.g. to sweep
     * seeds) works correctly.
     *
     * The instance setters below (setPopulationSize, etc.) are still
     * provided for direct unit-testing / manual instantiation outside
     * WorkflowSim, but have no effect when this class is run through
     * WorkflowPlanner -- use the static fields for that path.
     */
    public static int CONFIG_POPULATION_SIZE = 30;
    public static int CONFIG_GENERATION_COUNT = 100;
    public static Long CONFIG_RANDOM_SEED = null;

    /**
     * Snapshot of the most recently completed run, published at the end of
     * run(). CloudSim/WorkflowSim is single-threaded and each simulation
     * (CloudSim.init() ... CloudSim.stopSimulation()) constructs and runs
     * exactly one planning algorithm instance, so reading this immediately
     * after CloudSim.stopSimulation() and before starting the next
     * simulation is safe. Not safe for concurrent/parallel simulation runs
     * within a single JVM.
     */
    public static volatile LastRunMetrics lastRun = null;

    public static class LastRunMetrics {
        public double chosenMakespan;
        public double chosenCost;
        public int paretoFrontSize;
        /** Each element is {makespan, cost} for one member of the final Pareto front. */
        public List<double[]> paretoFrontPoints;
        public long searchWallClockMillis;
        public int populationSizeUsed;
        public int generationCountUsed;
        /** Number of objective evaluations (decoder calls made for population members). */
        public long evaluations;
        /** Non-dominated set of EVERY solution evaluated during the run (anytime archive). */
        public List<double[]> archiveFront;
        public List<int[]> archiveGenotypes;
        /** Convergence checkpoints: {generation, evaluations, cumulativeNanos} plus the population front. */
        public List<Checkpoint> trace;
    }

    public static class Checkpoint {
        public int generation;
        public long evaluations;
        public long nanos;
        public List<double[]> front;
    }

    /** Evaluate with the simulator-faithful decoder (FaithfulModel) instead of the analytic gap-filling decoder. */
    public static boolean CONFIG_FAITHFUL = false;
    protected FaithfulModel faithful;
    /** If > 0 (and the faithful decoder is on) every objective evaluation is the MEAN over this many fixed perturbed scenarios
     *  (task-length factors log-normal with the given sigma, machine speeds U(low,1) of peak): an expensive, uncertainty-aware objective. */
    public static int CONFIG_ROBUST_S = 0;
    public static double CONFIG_ROBUST_SIGMA = 0.25;
    public static double CONFIG_ROBUST_MIPS_LOW = 0.6;
    /** Risk-aware variant: > 0 makes the makespan objective the mean of the worst (1 - alpha) share of scenarios (CVaR_alpha). */
    public static double CONFIG_ROBUST_CVAR = 0.0;
    /** Probability that a machine is a straggler in a scenario (speed factor U(0.2, 0.5)); 0 = no stragglers. */
    public static double CONFIG_ROBUST_SLOWP = 0.0;
    protected double[][] robLen, robMips;
    protected long nominalDecodes = 0;

    public long nominalDecodeCount() {
        return nominalDecodes;
    }

    protected double[][][] scenarios(int S, long seed) {
        java.util.Random rs = new java.util.Random(seed);
        double[][] lm = new double[S][nTasks];
        double[][] mm = new double[S][nVms];
        for (int s = 0; s < S; s++) {
            for (int k = 0; k < nTasks; k++) {
                lm[s][k] = Math.exp(CONFIG_ROBUST_SIGMA * rs.nextGaussian() - CONFIG_ROBUST_SIGMA * CONFIG_ROBUST_SIGMA / 2);
            }
            for (int v = 0; v < nVms; v++) {
                if (CONFIG_ROBUST_SLOWP > 0 && rs.nextDouble() < CONFIG_ROBUST_SLOWP) {
                    mm[s][v] = 0.2 + 0.3 * rs.nextDouble();
                } else {
                    mm[s][v] = CONFIG_ROBUST_MIPS_LOW + (1.0 - CONFIG_ROBUST_MIPS_LOW) * rs.nextDouble();
                }
            }
        }
        return new double[][][]{lm, mm};
    }

    protected double[] expectedObjective(int[] g, double[][] lm, double[][] mm) {
        int S = lm.length;
        double[] ms = new double[S];
        double sc = 0, sm = 0;
        for (int s = 0; s < S; s++) {
            double[] r = faithful.decodePerturbed(g, lm[s], mm[s]);
            ms[s] = r[0];
            sm += r[0];
            sc += r[1];
        }
        if (CONFIG_ROBUST_CVAR > 0) {
            java.util.Arrays.sort(ms);
            int tail = Math.max(1, (int) Math.ceil((1.0 - CONFIG_ROBUST_CVAR) * S));
            double t = 0;
            for (int i = S - tail; i < S; i++) {
                t += ms[i];
            }
            return new double[]{t / tail, sc / S};
        }
        return new double[]{sm / S, sc / S};
    }

    /** Non-dominated front of the expected objectives of the given genotypes under FRESH scenarios (out-of-sample). */
    public List<double[]> outOfSampleFront(List<int[]> genos, int S) {
        double[][][] sc = scenarios(S, 424242L + 31L * nTasks);
        List<double[]> pts = new ArrayList<>();
        for (int[] g : genos) {
            if (g != null) {
                pts.add(expectedObjective(g, sc[0], sc[1]));
            }
        }
        List<double[]> out = new ArrayList<>();
        for (double[] q : pts) {
            boolean dom = false;
            for (double[] r : pts) {
                if (r != q && r[0] <= q[0] && r[1] <= q[1] && (r[0] < q[0] || r[1] < q[1])) {
                    dom = true;
                    break;
                }
            }
            if (!dom) {
                out.add(q);
            }
        }
        return out;
    }

    /** If non-null, run() only commits this genotype (used to replay a schedule inside the simulator). */
    public static int[] CONFIG_FORCED_GENOTYPE = null;

    /** Set true by experiment drivers to record convergence checkpoints. */
    public static boolean CONFIG_TRACE = false;
    protected static final int[] TRACE_GENS = {0, 1, 2, 3, 5, 7, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    protected long evalCount = 0;
    protected List<double[]> archiveFront = new ArrayList<>();
    protected List<Checkpoint> trace = new ArrayList<>();
    protected long runStartNanos;

    /** Counted objective evaluation: decode + evaluation counter + anytime non-dominated archive. */
    protected double[] evaluate(int[] genotype) {
        double[] mc = decode(genotype);
        evalCount++;
        lastAdded = observe(mc[0], mc[1], genotype);
        return mc;
    }

    /** Genotypes of the anytime archive, parallel to archiveFront (kept so fronts can be re-evaluated in the simulator). */
    protected List<int[]> archiveGeno = new ArrayList<>();
    /** True if the most recent evaluate() call added a new non-dominated point to the anytime archive. */
    protected boolean lastAdded = false;

    protected boolean observe(double m, double c) {
        return observe(m, c, null);
    }

    protected boolean observe(double m, double c, int[] g) {
        for (int a = 0; a < archiveFront.size(); a++) {
            double[] q = archiveFront.get(a);
            if (q[0] <= m && q[1] <= c) {
                return false;   // dominated or duplicate
            }
        }
        for (int a = archiveFront.size() - 1; a >= 0; a--) {
            double[] q = archiveFront.get(a);
            if (m <= q[0] && c <= q[1]) {
                archiveFront.remove(a);
                archiveGeno.remove(a);
            }
        }
        archiveFront.add(new double[]{m, c});
        archiveGeno.add(g == null ? null : g.clone());
        return true;
    }

    /** Non-dominated subset of the current population (decoder values). */
    protected List<double[]> currentFront() {
        List<double[]> out = new ArrayList<>();
        for (int i = 0; i < populationSize; i++) {
            boolean dom = false;
            for (int j = 0; j < populationSize && !dom; j++) {
                if (j != i && dominates(makespans[j], costs[j], makespans[i], costs[i])) {
                    dom = true;
                }
            }
            if (!dom) {
                out.add(new double[]{makespans[i], costs[i]});
            }
        }
        return out;
    }

    protected void recordCheckpoint(int gen) {
        if (!CONFIG_TRACE) {
            return;
        }
        for (int g : TRACE_GENS) {
            if (g == gen) {
                Checkpoint cp = new Checkpoint();
                cp.generation = gen;
                cp.evaluations = evalCount;
                cp.nanos = System.nanoTime() - runStartNanos;
                cp.front = currentFront();
                trace.add(cp);
            }
        }
    }

    /**
     * Optional warm-start schedules to seed the initial population with,
     * e.g. HEFT's and Min-Min's actual assignments. Keyed by CLOUDLET ID
     * (Task.getCloudletId()) rather than array position, because each
     * schedule is typically produced by a SEPARATE simulation run (its own
     * WorkflowPlanner, its own freshly-constructed planning algorithm
     * instance with its own internal task ordering) -- cloudlet IDs are
     * the one thing guaranteed to line up across separate runs on the
     * same DAX file, since WorkflowParser assigns them deterministically,
     * in document order, starting fresh from a per-instance counter every
     * time the file is parsed.
     *
     * Each entry in the outer list is one complete schedule: a map from
     * cloudlet ID to the VM ID (CondorVM.getId(), not a list index) that
     * task was assigned to. Set this from the benchmark driver right
     * before CloudSim.startSimulation(), after computing the schedule(s)
     * you want to seed with in a prior simulation run.
     *
     * If a seed schedule is missing an entry for some task in the current
     * workflow, or maps a task to a VM ID not present in the current VM
     * pool, that whole seed is skipped (logged, not silently dropped) and
     * search proceeds with whatever seeds remain plus random fill-in.
     */
    public static List<Map<Integer, Integer>> CONFIG_SEED_ASSIGNMENTS = null;

    // ---- problem data, set once at the start of run() ----
    protected List<Task> taskOrder;
    protected List<CondorVM> vmList;
    protected double averageBandwidth;
    protected Map<Task, Map<CondorVM, Double>> computationCosts;
    protected Map<Task, Map<Task, Double>> transferCosts;
    protected Random random;

    // ---- search state ----
    protected List<int[]> population;
    protected double[] makespans;
    protected double[] costs;
    protected int populationSize = CONFIG_POPULATION_SIZE;
    protected int generationCount = CONFIG_GENERATION_COUNT;
    private List<int[]> seedGenotypes;

    private static class Event {

        double start;
        double finish;

        Event(double start, double finish) {
            this.start = start;
            this.finish = finish;
        }
    }

    public LIWSAPlanningAlgorithm() {
        this.populationSize = CONFIG_POPULATION_SIZE;
        this.generationCount = CONFIG_GENERATION_COUNT;
        this.randomSeed = CONFIG_RANDOM_SEED;
    }

    public void setPopulationSize(int populationSize) {
        this.populationSize = populationSize;
    }

    public void setGenerationCount(int generationCount) {
        this.generationCount = generationCount;
    }

    public void setLambdaMix(double lambdaMix) {
        this.lambdaMix = lambdaMix;
    }

    public void setMutationRate(double mutationRate) {
        this.mutationRate = mutationRate;
    }

    public void setRandomSeed(long seed) {
        this.randomSeed = seed;
    }

    /**
     * Optional warm-start seeds (e.g. HEFT's and Min-Min's actual
     * assignments, translated into genotypes in this algorithm's task
     * order). Matches how the validated Python comparison was run, and
     * mirrors MLEAOPlanningAlgorithm's setSeedGenotypes for a fair,
     * identically-seeded comparison between the two.
     */
    public void setSeedGenotypes(List<int[]> seedGenotypes) {
        this.seedGenotypes = seedGenotypes;
    }

    @Override
    public void run() {
        Log.printLine("LIWSA planner running with " + getTaskList().size()
                + " tasks, " + getVmList().size() + " VMs.");
        if (CONFIG_FORCED_GENOTYPE != null) {
            setupProblem();
            commitAssignment(CONFIG_FORCED_GENOTYPE);
            return;
        }
        long searchStartMillis = System.currentTimeMillis();
        runStartNanos = System.nanoTime();

        setupProblem();

        initializePopulation();
        recordCheckpoint(0);
        double tau = calibrateTau();

        for (int gen = 0; gen < generationCount; gen++) {
            int[] frontNumber = new int[populationSize];
            List<List<Integer>> fronts = nonDominatedSort(frontNumber);

            List<Integer> elite = new ArrayList<>(fronts.get(0));
            int idx = 1;
            while (elite.size() < minEliteSize && idx < fronts.size()) {
                elite.addAll(fronts.get(idx));
                idx++;
            }

            int bestIndex = bestOf(fronts.get(0));

            for (int i = 0; i < populationSize; i++) {
                if (i == bestIndex) {
                    continue;
                }
                double density = localDensity(i, tau);
                double pSocial = (1 - lambdaMix) * ((double) gen / Math.max(generationCount, 1))
                        + lambdaMix * density;

                int[] child;
                if (random.nextDouble() > pSocial) {
                    child = solitaryMove(i, frontNumber);
                } else {
                    child = socialMove(i, frontNumber, elite);
                }
                mutate(child);

                double[] mc = evaluate(child);
                if (!dominates(makespans[i], costs[i], mc[0], mc[1])) {
                    population.set(i, child);
                    makespans[i] = mc[0];
                    costs[i] = mc[1];
                }
            }
            recordCheckpoint(gen + 1);
        }

        int[] finalFrontNumber = new int[populationSize];
        List<List<Integer>> finalFronts = nonDominatedSort(finalFrontNumber);
        int chosen = bestOf(finalFronts.get(0));

        // Capture metrics before committing, while makespans/costs/population
        // still hold the full final state.
        LastRunMetrics metrics = new LastRunMetrics();
        metrics.chosenMakespan = makespans[chosen];
        metrics.chosenCost = costs[chosen];
        metrics.paretoFrontSize = finalFronts.get(0).size();
        metrics.paretoFrontPoints = new ArrayList<>();
        for (int i : finalFronts.get(0)) {
            metrics.paretoFrontPoints.add(new double[]{makespans[i], costs[i]});
        }
        metrics.searchWallClockMillis = System.currentTimeMillis() - searchStartMillis;
        metrics.populationSizeUsed = populationSize;
        metrics.generationCountUsed = generationCount;
        metrics.evaluations = evalCount;
        metrics.archiveFront = new ArrayList<>(archiveFront);
        metrics.archiveGenotypes = new ArrayList<>(archiveGeno);
        metrics.trace = trace;
        lastRun = metrics;

        commitAssignment(population.get(chosen));

        Log.printLine("LIWSA finished. Pareto front size: " + finalFronts.get(0).size()
                + ", chosen makespan=" + makespans[chosen] + ", cost=" + costs[chosen]);
    }

    /** Builds all per-instance data (VM list, task order, cost tables, fast decoder tables). */
    protected void setupProblem() {
        vmList = new ArrayList<>();
        for (Object vmObject : getVmList()) {
            vmList.add((CondorVM) vmObject);
        }
        random = (randomSeed != null) ? new Random(randomSeed) : new Random();
        averageBandwidth = calculateAverageBandwidth();
        taskOrder = topologicalOrder(getTaskList());
        calculateComputationCosts();
        calculateTransferCosts();
        buildFastTables();
        faithful = CONFIG_FAITHFUL ? new FaithfulModel(taskOrder, vmList) : null;
        robLen = null;
        robMips = null;
        if (faithful != null && CONFIG_ROBUST_S > 0) {
            double[][][] sc = scenarios(CONFIG_ROBUST_S, 7919L * nTasks + 131L * nVms + 17L);
            robLen = sc[0];
            robMips = sc[1];
        }
    }

    protected int bestOf(List<Integer> indices) {
        int best = indices.get(0);
        for (int i : indices) {
            if (makespans[i] < makespans[best]
                    || (makespans[i] == makespans[best] && costs[i] < costs[best])) {
                best = i;
            }
        }
        return best;
    }

    // ---------------------------------------------------------------
    // Problem data setup (computation/transfer costs mirror
    // HEFTPlanningAlgorithm exactly, so LIWSA and HEFT share identical
    // low-level timing assumptions)
    // ---------------------------------------------------------------
    protected double calculateAverageBandwidth() {
        double avg = 0.0;
        for (CondorVM vm : vmList) {
            avg += vm.getBw();
        }
        return avg / vmList.size();
    }

    protected int computeDepth(Task t, Map<Task, Integer> depth) {
        if (depth.containsKey(t)) {
            return depth.get(t);
        }
        int d = 0;
        for (Object parentObj : t.getParentList()) {
            Task p = (Task) parentObj;
            d = Math.max(d, computeDepth(p, depth) + 1);
        }
        depth.put(t, d);
        return d;
    }

    protected List<Task> topologicalOrder(List<Task> tasks) {
        final Map<Task, Integer> depth = new HashMap<>();
        for (Task t : tasks) {
            computeDepth(t, depth);
        }
        List<Task> order = new ArrayList<>(tasks);
        Collections.sort(order, (a, b) -> {
            int da = depth.get(a);
            int db = depth.get(b);
            if (da != db) {
                return Integer.compare(da, db);
            }
            double la = a.getCloudletTotalLength();
            double lb = b.getCloudletTotalLength();
            if (la != lb) {
                return Double.compare(lb, la);
            }
            return Integer.compare(a.getCloudletId(), b.getCloudletId());
        });
        return order;
    }

    protected void calculateComputationCosts() {
        computationCosts = new HashMap<>();
        for (Task task : getTaskList()) {
            Map<CondorVM, Double> costsForTask = new HashMap<>();
            for (CondorVM vm : vmList) {
                costsForTask.put(vm, task.getCloudletTotalLength() / vm.getMips());
            }
            computationCosts.put(task, costsForTask);
        }
    }

    protected void calculateTransferCosts() {
        transferCosts = new HashMap<>();
        for (Task task : getTaskList()) {
            transferCosts.put(task, new HashMap<Task, Double>());
        }
        for (Task parent : getTaskList()) {
            for (Object childObj : parent.getChildList()) {
                Task child = (Task) childObj;
                transferCosts.get(parent).put(child, calculateTransferCost(parent, child));
            }
        }
    }

    private double calculateTransferCost(Task parent, Task child) {
        List<FileItem> parentFiles = parent.getFileList();
        List<FileItem> childFiles = child.getFileList();
        double acc = 0.0;
        for (FileItem parentFile : parentFiles) {
            if (parentFile.getType() != Parameters.FileType.OUTPUT) {
                continue;
            }
            for (FileItem childFile : childFiles) {
                if (childFile.getType() == Parameters.FileType.INPUT
                        && childFile.getName().equals(parentFile.getName())) {
                    acc += childFile.getSize();
                    break;
                }
            }
        }
        // file size is in bytes, acc converted to MB, then to Mbit -- identical
        // to HEFTPlanningAlgorithm.calculateTransferCost
        acc = acc / Consts.MILLION;
        return acc * 8 / averageBandwidth;
    }

    // ---------------------------------------------------------------
    // Decoder: genotype -> (makespan, cost), feasible by construction
    // (topological order + insertion-based per-VM scheduling, so there is
    // no separate repair step)
    // ---------------------------------------------------------------
    private double findFinishTime(List<Event> sched, double readyTime, double duration, boolean occupySlot) {
        if (sched.isEmpty()) {
            if (occupySlot) {
                sched.add(new Event(readyTime, readyTime + duration));
            }
            return readyTime + duration;
        }

        if (sched.size() == 1) {
            double start;
            int pos;
            if (readyTime >= sched.get(0).finish) {
                pos = 1;
                start = readyTime;
            } else if (readyTime + duration <= sched.get(0).start) {
                pos = 0;
                start = readyTime;
            } else {
                pos = 1;
                start = sched.get(0).finish;
            }
            if (occupySlot) {
                sched.add(pos, new Event(start, start + duration));
            }
            return start + duration;
        }

        double start = Math.max(readyTime, sched.get(sched.size() - 1).finish);
        double finish = start + duration;
        int pos = sched.size();
        int i = sched.size() - 1;
        int j = sched.size() - 2;
        while (j >= 0) {
            Event current = sched.get(i);
            Event previous = sched.get(j);
            if (readyTime > previous.finish) {
                if (readyTime + duration <= current.start) {
                    start = readyTime;
                    finish = readyTime + duration;
                }
                break;
            }
            if (previous.finish + duration <= current.start) {
                start = previous.finish;
                finish = previous.finish + duration;
                pos = i;
            }
            i--;
            j--;
        }

        if (readyTime + duration <= sched.get(0).start) {
            pos = 0;
            start = readyTime;
            if (occupySlot) {
                sched.add(pos, new Event(start, start + duration));
            }
            return start + duration;
        }

        if (occupySlot) {
            sched.add(pos, new Event(start, finish));
        }
        return finish;
    }

    // ---- fast decoder tables (built once per instance; arithmetic identical to the map-based decoder) ----
    protected int nTasks;
    protected int nVms;
    protected int[][] parentIdx;
    protected double[][] parentTc;
    protected double[][] dur;
    protected double[] vmRate;
    private double[][] evStart;
    private double[][] evFinish;
    private int[] evCount;
    private double[] finishBuf;

    protected void buildFastTables() {
        nTasks = taskOrder.size();
        nVms = vmList.size();
        Map<Task, Integer> pos = new HashMap<>();
        for (int k = 0; k < nTasks; k++) {
            pos.put(taskOrder.get(k), k);
        }
        parentIdx = new int[nTasks][];
        parentTc = new double[nTasks][];
        dur = new double[nTasks][nVms];
        for (int k = 0; k < nTasks; k++) {
            Task task = taskOrder.get(k);
            List<?> parents = task.getParentList();
            parentIdx[k] = new int[parents.size()];
            parentTc[k] = new double[parents.size()];
            for (int q = 0; q < parents.size(); q++) {
                Task parent = (Task) parents.get(q);
                parentIdx[k][q] = pos.get(parent);
                Double tc = transferCosts.get(parent).get(task);
                parentTc[k][q] = (tc != null) ? tc : 0.0;
            }
            for (int v = 0; v < nVms; v++) {
                dur[k][v] = computationCosts.get(task).get(vmList.get(v));
            }
        }
        vmRate = new double[nVms];
        for (int v = 0; v < nVms; v++) {
            vmRate[v] = vmList.get(v).getCost();
        }
        evStart = new double[nVms][nTasks + 1];
        evFinish = new double[nVms][nTasks + 1];
        evCount = new int[nVms];
        finishBuf = new double[nTasks];
    }

    private void evInsert(int v, int pos, double start, double finish) {
        int cnt = evCount[v];
        double[] s = evStart[v];
        double[] f = evFinish[v];
        if (pos < cnt) {
            System.arraycopy(s, pos, s, pos + 1, cnt - pos);
            System.arraycopy(f, pos, f, pos + 1, cnt - pos);
        }
        s[pos] = start;
        f[pos] = finish;
        evCount[v] = cnt + 1;
    }

    /** Insertion-based gap-filling placement, identical logic to findFinishTime() on arrays. */
    private double placeFast(int v, double ready, double duration) {
        int cnt = evCount[v];
        double[] s = evStart[v];
        double[] f = evFinish[v];
        if (cnt == 0) {
            evInsert(v, 0, ready, ready + duration);
            return ready + duration;
        }
        if (cnt == 1) {
            double start;
            int pos;
            if (ready >= f[0]) {
                pos = 1;
                start = ready;
            } else if (ready + duration <= s[0]) {
                pos = 0;
                start = ready;
            } else {
                pos = 1;
                start = f[0];
            }
            evInsert(v, pos, start, start + duration);
            return start + duration;
        }
        double start = Math.max(ready, f[cnt - 1]);
        double finish = start + duration;
        int pos = cnt;
        int i = cnt - 1;
        int j = cnt - 2;
        while (j >= 0) {
            double curStart = s[i];
            double prevFinish = f[j];
            if (ready > prevFinish) {
                if (ready + duration <= curStart) {
                    start = ready;
                    finish = ready + duration;
                }
                break;
            }
            if (prevFinish + duration <= curStart) {
                start = prevFinish;
                finish = prevFinish + duration;
                pos = i;
            }
            i--;
            j--;
        }
        if (ready + duration <= s[0]) {
            evInsert(v, 0, ready, ready + duration);
            return ready + duration;
        }
        evInsert(v, pos, start, finish);
        return finish;
    }

    protected double[] decode(int[] genotype) {
        if (faithful != null) {
            if (robLen != null) {
                return expectedObjective(genotype, robLen, robMips);
            }
            return faithful.decode(genotype);
        }
        java.util.Arrays.fill(evCount, 0);
        double cost = 0.0;
        for (int k = 0; k < nTasks; k++) {
            int v = genotype[k];
            double ready = 0.0;
            int[] par = parentIdx[k];
            double[] ptc = parentTc[k];
            for (int q = 0; q < par.length; q++) {
                int pk = par[q];
                double pf = finishBuf[pk];
                if (genotype[pk] != v) {
                    pf += ptc[q];
                }
                if (pf > ready) {
                    ready = pf;
                }
            }
            double d = dur[k][v];
            finishBuf[k] = placeFast(v, ready, d);
            cost += d * vmRate[v];
        }
        double makespan = 0.0;
        for (int k = 0; k < nTasks; k++) {
            if (finishBuf[k] > makespan) {
                makespan = finishBuf[k];
            }
        }
        return new double[]{makespan, cost};
    }

    /** Reference (map-based) decoder kept only for the regression test of the fast one. */
    protected double[] decodeReference(int[] genotype) {
        Map<CondorVM, List<Event>> schedules = new HashMap<>();
        for (CondorVM vm : vmList) {
            schedules.put(vm, new ArrayList<Event>());
        }
        Map<Task, Double> finish = new HashMap<>();
        Map<Task, CondorVM> assignedVm = new HashMap<>();
        double cost = 0.0;
        for (int k = 0; k < taskOrder.size(); k++) {
            Task task = taskOrder.get(k);
            CondorVM vm = vmList.get(genotype[k]);
            double ready = 0.0;
            for (Object parentObj : task.getParentList()) {
                Task parent = (Task) parentObj;
                double pf = finish.get(parent);
                if (assignedVm.get(parent) != vm) {
                    Double tc = transferCosts.get(parent).get(task);
                    pf += (tc != null) ? tc : 0.0;
                }
                ready = Math.max(ready, pf);
            }
            double duration = computationCosts.get(task).get(vm);
            double fin = findFinishTime(schedules.get(vm), ready, duration, true);
            finish.put(task, fin);
            assignedVm.put(task, vm);
            cost += duration * vm.getCost();
        }
        double makespan = 0.0;
        for (double f : finish.values()) {
            makespan = Math.max(makespan, f);
        }
        return new double[]{makespan, cost};
    }

    protected void commitAssignment(int[] genotype) {
        for (int k = 0; k < taskOrder.size(); k++) {
            taskOrder.get(k).setVmId(vmList.get(genotype[k]).getId());
        }
    }

    // ---------------------------------------------------------------
    // Population, Pareto ranking, density, and the two locust-derived
    // movement operators
    // ---------------------------------------------------------------
    /**
     * Hook for subclasses: return any warm-start genotypes to seed the
     * initial population with before filling the rest randomly.
     * Returns, in order: any genotypes set via setSeedGenotypes() (array-
     * position based, for direct/manual use outside WorkflowSim), then
     * any schedules in CONFIG_SEED_ASSIGNMENTS translated into genotypes
     * via cloudlet-ID lookup (the path used when run through WorkflowSim,
     * e.g. to warm-start from HEFT's/Min-Min's actual schedules).
     * LIWSAMLPlanningAlgorithm further appends ML-biased starting points
     * on top of whatever this returns.
     */
    protected List<int[]> generateSeedGenotypes() {
        List<int[]> seeds = (seedGenotypes != null) ? new ArrayList<>(seedGenotypes) : new ArrayList<>();
        seeds.addAll(buildSeedsFromAssignments());
        return seeds;
    }

    /**
     * Converts each schedule in CONFIG_SEED_ASSIGNMENTS (cloudlet ID -> VM
     * ID) into a genotype matching this instance's own taskOrder and
     * vmList. A schedule is skipped entirely, with a log line, if any
     * task in the current workflow has no entry in it, or if it maps a
     * task to a VM ID not present in the current VM pool -- both would
     * indicate the seed was computed for a different workflow or VM
     * configuration than the one currently running.
     */
    private List<int[]> buildSeedsFromAssignments() {
        List<int[]> result = new ArrayList<>();
        if (CONFIG_SEED_ASSIGNMENTS == null || CONFIG_SEED_ASSIGNMENTS.isEmpty()) {
            return result;
        }
        Map<Integer, Integer> vmIdToIndex = new HashMap<>();
        for (int idx = 0; idx < vmList.size(); idx++) {
            vmIdToIndex.put(vmList.get(idx).getId(), idx);
        }

        for (Map<Integer, Integer> assignment : CONFIG_SEED_ASSIGNMENTS) {
            int[] genotype = new int[taskOrder.size()];
            boolean valid = true;
            for (int k = 0; k < taskOrder.size(); k++) {
                int cloudletId = taskOrder.get(k).getCloudletId();
                Integer vmId = assignment.get(cloudletId);
                if (vmId == null || !vmIdToIndex.containsKey(vmId)) {
                    valid = false;
                    break;
                }
                genotype[k] = vmIdToIndex.get(vmId);
            }
            if (valid) {
                result.add(genotype);
            } else {
                Log.printLine("LIWSA: skipped a seed schedule that did not match "
                    + "the current workflow/VM pool (incomplete cloudlet-ID mapping).");
            }
        }
        return result;
    }

    protected void initializePopulation() {
        population = new ArrayList<>();
        for (int[] g : generateSeedGenotypes()) {
            population.add(g.clone());
        }
        int n = taskOrder.size();
        while (population.size() < populationSize) {
            int[] genotype = new int[n];
            for (int k = 0; k < n; k++) {
                genotype[k] = random.nextInt(vmList.size());
            }
            population.add(genotype);
        }
        if (population.size() > populationSize) {
            population = new ArrayList<>(population.subList(0, populationSize));
        }
        makespans = new double[populationSize];
        costs = new double[populationSize];
        for (int i = 0; i < populationSize; i++) {
            double[] mc = evaluate(population.get(i));
            makespans[i] = mc[0];
            costs[i] = mc[1];
        }
    }

    protected double calibrateTau() {
        List<Double> sample = new ArrayList<>();
        for (int t = 0; t < 200; t++) {
            int a = random.nextInt(populationSize);
            int b = random.nextInt(populationSize);
            if (a != b) {
                sample.add(hamming(population.get(a), population.get(b)));
            }
        }
        if (sample.isEmpty()) {
            return 0.3;
        }
        Collections.sort(sample);
        return sample.get(sample.size() / 2);
    }

    protected double hamming(int[] a, int[] b) {
        int diff = 0;
        for (int k = 0; k < a.length; k++) {
            if (a[k] != b[k]) {
                diff++;
            }
        }
        return (double) diff / a.length;
    }

    protected double kernel(double d) {
        double raw = kernelF * Math.exp(-d / kernelL) - Math.exp(-d);
        return 1.0 / (1.0 + Math.exp(-raw));
    }

    protected boolean dominates(double m1, double c1, double m2, double c2) {
        return (m1 <= m2 && c1 <= c2) && (m1 < m2 || c1 < c2);
    }

    protected List<List<Integer>> nonDominatedSort(int[] frontNumberOut) {
        int n = populationSize;
        int[] domCount = new int[n];
        List<List<Integer>> dominatedBy = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            dominatedBy.add(new ArrayList<Integer>());
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    continue;
                }
                if (dominates(makespans[i], costs[i], makespans[j], costs[j])) {
                    dominatedBy.get(i).add(j);
                } else if (dominates(makespans[j], costs[j], makespans[i], costs[i])) {
                    domCount[i]++;
                }
            }
        }
        List<List<Integer>> fronts = new ArrayList<>();
        List<Integer> current = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            if (domCount[i] == 0) {
                current.add(i);
            }
        }
        int rank = 0;
        while (!current.isEmpty()) {
            for (int i : current) {
                frontNumberOut[i] = rank;
            }
            fronts.add(current);
            List<Integer> next = new ArrayList<>();
            for (int i : current) {
                for (int j : dominatedBy.get(i)) {
                    domCount[j]--;
                    if (domCount[j] == 0) {
                        next.add(j);
                    }
                }
            }
            current = next;
            rank++;
        }
        return fronts;
    }

    protected double localDensity(int i, double tau) {
        int n = populationSize - 1;
        if (n <= 0) {
            return 0.0;
        }
        int count = 0;
        for (int j = 0; j < populationSize; j++) {
            if (j != i && hamming(population.get(i), population.get(j)) < tau) {
                count++;
            }
        }
        return (double) count / n;
    }

    protected int weightedChoice(List<Integer> options, List<Double> weights) {
        double total = 0.0;
        for (double w : weights) {
            total += w;
        }
        if (total <= 0) {
            return options.get(random.nextInt(options.size()));
        }
        double r = random.nextDouble() * total;
        double acc = 0.0;
        for (int k = 0; k < options.size(); k++) {
            acc += weights.get(k);
            if (r <= acc) {
                return options.get(k);
            }
        }
        return options.get(options.size() - 1);
    }

    /**
     * Solitary phase: every other individual casts a signed,
     * distance-weighted vote on each task's VM (positive if better-ranked,
     * negative if worse, zero if on the same Pareto front), then the VM
     * for that task is resampled from a softmax over the accumulated
     * votes. This is the discrete stand-in for summing continuous
     * attraction/repulsion force vectors.
     */
    protected int[] solitaryMove(int i, int[] frontNumber) {
        int n = taskOrder.size();
        int[] child = population.get(i).clone();
        // The signed kernel weight of every voter depends only on (i, j), not on the task k,
        // so it is computed once per individual instead of once per task (result-identical).
        double[] voteWeight = new double[populationSize];
        boolean[] voter = new boolean[populationSize];
        int[] voterList = new int[populationSize];
        int nVoters = 0;
        for (int j = 0; j < populationSize; j++) {
            if (j == i || frontNumber[j] == frontNumber[i]) {
                continue;
            }
            double sign = (frontNumber[j] < frontNumber[i]) ? 1.0 : -1.0;
            double d = hamming(population.get(i), population.get(j));
            voteWeight[j] = sign * kernel(d);
            voter[j] = true;
            voterList[nVoters++] = j;
        }
        for (int k = 0; k < n; k++) {
            Map<Integer, Double> votes = new HashMap<>();
            for (int q = 0; q < nVoters; q++) {
                int j = voterList[q];
                int v = population.get(j)[k];
                votes.put(v, votes.getOrDefault(v, 0.0) + voteWeight[j]);
            }
            if (!votes.isEmpty() && random.nextDouble() < blendProbability) {
                List<Integer> vmsList = new ArrayList<>(votes.keySet());
                double max = Double.NEGATIVE_INFINITY;
                for (int v : vmsList) {
                    max = Math.max(max, votes.get(v));
                }
                List<Double> probs = new ArrayList<>();
                double sum = 0.0;
                for (int v : vmsList) {
                    double e = Math.exp(votes.get(v) - max);
                    probs.add(e);
                    sum += e;
                }
                for (int idx = 0; idx < probs.size(); idx++) {
                    probs.set(idx, probs.get(idx) / sum);
                }
                child[k] = weightedChoice(vmsList, probs);
            }
        }
        return child;
    }

    /**
     * Social phase: roulette-selected attraction toward one member of the
     * current elite (front-1, expanded if too small) set only, then a
     * per-task copy-by-probability toward that single partner.
     */
    protected int[] socialMove(int i, int[] frontNumber, List<Integer> elite) {
        List<Integer> candidates = new ArrayList<>();
        for (int e : elite) {
            if (e != i) {
                candidates.add(e);
            }
        }
        if (candidates.isEmpty()) {
            return population.get(i).clone();
        }
        List<Double> weights = new ArrayList<>();
        for (int e : candidates) {
            double d = hamming(population.get(i), population.get(e));
            weights.add(kernel(d) / (frontNumber[e] + 1));
        }
        int partner = weightedChoice(candidates, weights);
        int[] Y = population.get(partner);
        double dIY = hamming(population.get(i), Y);
        double pCopy = Math.max(0.0, Math.min(1.0, copyAlpha * kernel(dIY)));

        int[] child = population.get(i).clone();
        for (int k = 0; k < child.length; k++) {
            if (random.nextDouble() < pCopy) {
                child[k] = Y[k];
            }
        }
        return child;
    }

    protected void mutate(int[] child) {
        for (int k = 0; k < child.length; k++) {
            if (random.nextDouble() < mutationRate) {
                child[k] = random.nextInt(vmList.size());
            }
        }
    }
}
