package org.workflowsim.planning;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Canonical multi-objective evolutionary baselines (NSGA-II, SPEA2, MOEA/D) on EXACTLY the same problem
 * representation as LIWSA: integer genotype (VM index per task in the fixed topological order), the same
 * insertion-based decoder, the same initial-population construction (including any injected warm-start
 * individuals) and the same objective-evaluation budget. Only the search operators differ.
 *
 * Operators (categorical genes, so real-coded SBX / polynomial mutation do not apply): uniform crossover
 * (probability CONFIG_CROSSOVER_PROB) and per-gene random-reset mutation (rate CONFIG_MUT_RATE, a negative
 * value meaning 1/n). Objective-evaluation budget: initial population + (generationCount x (populationSize-1))
 * evaluations, i.e. 2930 for P=30, G=100, which is what LIWSA itself spends.
 */
public abstract class MOEAHosts extends LIWSASeededPlanningAlgorithm {

    public static double CONFIG_CROSSOVER_PROB = 0.9;
    public static double CONFIG_MUT_RATE = 0.02;

    protected abstract void evolve();

    @Override
    public void run() {
        long searchStartMillis = System.currentTimeMillis();
        runStartNanos = System.nanoTime();
        setupProblem();
        initializePopulation();
        recordCheckpoint(0);
        evolve();
        publish(searchStartMillis);
    }

    protected void publish(long searchStartMillis) {
        int[] fn = new int[populationSize];
        List<List<Integer>> fronts = nonDominatedSort(fn);
        int chosen = bestOf(fronts.get(0));
        LastRunMetrics m = new LastRunMetrics();
        m.chosenMakespan = makespans[chosen];
        m.chosenCost = costs[chosen];
        m.paretoFrontSize = fronts.get(0).size();
        m.paretoFrontPoints = new ArrayList<>();
        for (int i : fronts.get(0)) {
            m.paretoFrontPoints.add(new double[]{makespans[i], costs[i]});
        }
        m.searchWallClockMillis = System.currentTimeMillis() - searchStartMillis;
        m.populationSizeUsed = populationSize;
        m.generationCountUsed = generationCount;
        m.evaluations = evalCount;
        m.archiveFront = new ArrayList<>(archiveFront);
        m.archiveGenotypes = new ArrayList<>(archiveGeno);
        m.trace = trace;
        lastRun = m;
        commitAssignment(population.get(chosen));
    }

    // ------------------------------------------------------------------ shared operators
    protected int[] crossoverMutate(int[] a, int[] b) {
        int n = a.length;
        int[] child = a.clone();
        if (random.nextDouble() < CONFIG_CROSSOVER_PROB) {
            for (int k = 0; k < n; k++) {
                if (random.nextBoolean()) {
                    child[k] = b[k];
                }
            }
        }
        double pm = (CONFIG_MUT_RATE < 0) ? 1.0 / n : CONFIG_MUT_RATE;
        for (int k = 0; k < n; k++) {
            if (random.nextDouble() < pm) {
                child[k] = random.nextInt(vmList.size());
            }
        }
        return child;
    }

    /** Pareto rank (0 = non-dominated) of each point; O(N^2). */
    protected static int[] paretoRanks(double[] m, double[] c) {
        int n = m.length;
        int[] rank = new int[n];
        int[] domCount = new int[n];
        List<List<Integer>> dom = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            dom.add(new ArrayList<>());
        }
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    continue;
                }
                boolean iDomJ = m[i] <= m[j] && c[i] <= c[j] && (m[i] < m[j] || c[i] < c[j]);
                if (iDomJ) {
                    dom.get(i).add(j);
                } else if (m[j] <= m[i] && c[j] <= c[i] && (m[j] < m[i] || c[j] < c[i])) {
                    domCount[i]++;
                }
            }
        }
        List<Integer> cur = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            if (domCount[i] == 0) {
                cur.add(i);
            }
        }
        int r = 0;
        while (!cur.isEmpty()) {
            List<Integer> next = new ArrayList<>();
            for (int i : cur) {
                rank[i] = r;
                for (int j : dom.get(i)) {
                    if (--domCount[j] == 0) {
                        next.add(j);
                    }
                }
            }
            cur = next;
            r++;
        }
        return rank;
    }

    protected void setPopulation(List<int[]> pop, double[] m, double[] c) {
        population = new ArrayList<>(pop);
        makespans = m.clone();
        costs = c.clone();
    }

    // ================================================================== NSGA-II
    public static class NSGA2 extends MOEAHosts {
        @Override
        protected void evolve() {
            int P = populationSize;
            int lambda = P - 1;
            List<int[]> pop = new ArrayList<>(population);
            double[] pm = makespans.clone();
            double[] pc = costs.clone();
            for (int gen = 0; gen < generationCount; gen++) {
                int[] rank = paretoRanks(pm, pc);
                double[] crowd = crowding(pm, pc, rank);
                List<int[]> union = new ArrayList<>(pop);
                double[] um = Arrays.copyOf(pm, P + lambda);
                double[] uc = Arrays.copyOf(pc, P + lambda);
                for (int o = 0; o < lambda; o++) {
                    int a = tournament(rank, crowd, P);
                    int b = tournament(rank, crowd, P);
                    int[] child = crossoverMutate(pop.get(a), pop.get(b));
                    double[] mc = evaluate(child);
                    union.add(child);
                    um[P + o] = mc[0];
                    uc[P + o] = mc[1];
                }
                int[] urank = paretoRanks(um, uc);
                double[] ucrowd = crowding(um, uc, urank);
                Integer[] order = new Integer[P + lambda];
                for (int i = 0; i < order.length; i++) {
                    order[i] = i;
                }
                Arrays.sort(order, (x, y) -> urank[x] != urank[y] ? Integer.compare(urank[x], urank[y])
                        : Double.compare(ucrowd[y], ucrowd[x]));
                List<int[]> np = new ArrayList<>();
                double[] nm = new double[P];
                double[] nc = new double[P];
                for (int i = 0; i < P; i++) {
                    np.add(union.get(order[i]));
                    nm[i] = um[order[i]];
                    nc[i] = uc[order[i]];
                }
                pop = np;
                pm = nm;
                pc = nc;
                setPopulation(pop, pm, pc);
                recordCheckpoint(gen + 1);
            }
            setPopulation(pop, pm, pc);
        }

        private int tournament(int[] rank, double[] crowd, int n) {
            int a = random.nextInt(n);
            int b = random.nextInt(n);
            if (rank[a] != rank[b]) {
                return rank[a] < rank[b] ? a : b;
            }
            if (crowd[a] != crowd[b]) {
                return crowd[a] > crowd[b] ? a : b;
            }
            return random.nextBoolean() ? a : b;
        }

        /** Crowding distance computed within each Pareto rank, objectives scaled by their range. */
        private double[] crowding(double[] m, double[] c, int[] rank) {
            int n = m.length;
            double[] d = new double[n];
            int maxRank = 0;
            for (int r : rank) {
                maxRank = Math.max(maxRank, r);
            }
            for (int r = 0; r <= maxRank; r++) {
                List<Integer> idx = new ArrayList<>();
                for (int i = 0; i < n; i++) {
                    if (rank[i] == r) {
                        idx.add(i);
                    }
                }
                for (int obj = 0; obj < 2; obj++) {
                    final double[] f = (obj == 0) ? m : c;
                    idx.sort((x, y) -> Double.compare(f[x], f[y]));
                    double range = f[idx.get(idx.size() - 1)] - f[idx.get(0)];
                    d[idx.get(0)] = Double.POSITIVE_INFINITY;
                    d[idx.get(idx.size() - 1)] = Double.POSITIVE_INFINITY;
                    if (range <= 0) {
                        continue;
                    }
                    for (int q = 1; q + 1 < idx.size(); q++) {
                        d[idx.get(q)] += (f[idx.get(q + 1)] - f[idx.get(q - 1)]) / range;
                    }
                }
            }
            return d;
        }
    }

    // ================================================================== SPEA2
    public static class SPEA2 extends MOEAHosts {
        @Override
        protected void evolve() {
            int P = populationSize;
            int lambda = P - 1;
            List<int[]> arch = new ArrayList<>();
            double[] am = new double[0];
            double[] ac = new double[0];
            List<int[]> q = new ArrayList<>(population);
            double[] qm = makespans.clone();
            double[] qc = costs.clone();
            for (int gen = 0; gen <= generationCount; gen++) {
                // environmental selection from Q union A
                int nU = q.size() + arch.size();
                List<int[]> union = new ArrayList<>(q);
                union.addAll(arch);
                double[] um = new double[nU];
                double[] uc = new double[nU];
                System.arraycopy(qm, 0, um, 0, q.size());
                System.arraycopy(qc, 0, uc, 0, q.size());
                System.arraycopy(am, 0, um, q.size(), arch.size());
                System.arraycopy(ac, 0, uc, q.size(), arch.size());
                double[] fit = fitness(um, uc);
                List<Integer> sel = new ArrayList<>();
                for (int i = 0; i < nU; i++) {
                    if (fit[i] < 1.0) {
                        sel.add(i);
                    }
                }
                if (sel.size() < P) {
                    List<Integer> rest = new ArrayList<>();
                    for (int i = 0; i < nU; i++) {
                        if (fit[i] >= 1.0) {
                            rest.add(i);
                        }
                    }
                    rest.sort((x, y) -> Double.compare(fit[x], fit[y]));
                    for (int i = 0; sel.size() < P && i < rest.size(); i++) {
                        sel.add(rest.get(i));
                    }
                } else {
                    truncate(sel, um, uc, P);
                }
                arch = new ArrayList<>();
                am = new double[sel.size()];
                ac = new double[sel.size()];
                double[] af = new double[sel.size()];
                for (int i = 0; i < sel.size(); i++) {
                    arch.add(union.get(sel.get(i)));
                    am[i] = um[sel.get(i)];
                    ac[i] = uc[sel.get(i)];
                    af[i] = fit[sel.get(i)];
                }
                setPopulation(arch, am, ac);
                populationSize = arch.size();
                recordCheckpoint(gen);
                if (gen == generationCount) {
                    break;
                }
                q = new ArrayList<>();
                qm = new double[lambda];
                qc = new double[lambda];
                for (int o = 0; o < lambda; o++) {
                    int a = tourn(af);
                    int b = tourn(af);
                    int[] child = crossoverMutate(arch.get(a), arch.get(b));
                    double[] mc = evaluate(child);
                    q.add(child);
                    qm[o] = mc[0];
                    qc[o] = mc[1];
                }
            }
        }

        private int tourn(double[] f) {
            int a = random.nextInt(f.length);
            int b = random.nextInt(f.length);
            return f[a] <= f[b] ? a : b;
        }

        /** SPEA2 fitness F = R + D (lower is better; F < 1 means non-dominated). */
        private double[] fitness(double[] m, double[] c) {
            int n = m.length;
            int[] strength = new int[n];
            boolean[][] dom = new boolean[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i != j && m[i] <= m[j] && c[i] <= c[j] && (m[i] < m[j] || c[i] < c[j])) {
                        dom[i][j] = true;
                        strength[i]++;
                    }
                }
            }
            double mn0 = Double.MAX_VALUE, mx0 = -Double.MAX_VALUE, mn1 = Double.MAX_VALUE, mx1 = -Double.MAX_VALUE;
            for (int i = 0; i < n; i++) {
                mn0 = Math.min(mn0, m[i]);
                mx0 = Math.max(mx0, m[i]);
                mn1 = Math.min(mn1, c[i]);
                mx1 = Math.max(mx1, c[i]);
            }
            double r0 = Math.max(mx0 - mn0, 1e-12);
            double r1 = Math.max(mx1 - mn1, 1e-12);
            int k = Math.max(1, (int) Math.sqrt(n));
            double[] fit = new double[n];
            for (int i = 0; i < n; i++) {
                double raw = 0;
                for (int j = 0; j < n; j++) {
                    if (dom[j][i]) {
                        raw += strength[j];
                    }
                }
                double[] dist = new double[n - 1];
                int q = 0;
                for (int j = 0; j < n; j++) {
                    if (j != i) {
                        double d0 = (m[i] - m[j]) / r0;
                        double d1 = (c[i] - c[j]) / r1;
                        dist[q++] = Math.sqrt(d0 * d0 + d1 * d1);
                    }
                }
                Arrays.sort(dist);
                double sk = dist.length == 0 ? 0 : dist[Math.min(k, dist.length) - 1];
                fit[i] = raw + 1.0 / (sk + 2.0);
            }
            return fit;
        }

        /** Archive truncation: repeatedly drop the member with the lexicographically smallest nearest-neighbour distances. */
        private void truncate(List<Integer> sel, double[] m, double[] c, int target) {
            double mn0 = Double.MAX_VALUE, mx0 = -Double.MAX_VALUE, mn1 = Double.MAX_VALUE, mx1 = -Double.MAX_VALUE;
            for (int i : sel) {
                mn0 = Math.min(mn0, m[i]);
                mx0 = Math.max(mx0, m[i]);
                mn1 = Math.min(mn1, c[i]);
                mx1 = Math.max(mx1, c[i]);
            }
            double r0 = Math.max(mx0 - mn0, 1e-12);
            double r1 = Math.max(mx1 - mn1, 1e-12);
            while (sel.size() > target) {
                int worst = -1;
                double[] worstVec = null;
                for (int a = 0; a < sel.size(); a++) {
                    double[] vec = new double[sel.size() - 1];
                    int q = 0;
                    for (int b = 0; b < sel.size(); b++) {
                        if (a != b) {
                            double d0 = (m[sel.get(a)] - m[sel.get(b)]) / r0;
                            double d1 = (c[sel.get(a)] - c[sel.get(b)]) / r1;
                            vec[q++] = Math.sqrt(d0 * d0 + d1 * d1);
                        }
                    }
                    Arrays.sort(vec);
                    if (worstVec == null || lexLess(vec, worstVec)) {
                        worstVec = vec;
                        worst = a;
                    }
                }
                sel.remove(worst);
            }
        }

        private boolean lexLess(double[] x, double[] y) {
            for (int i = 0; i < x.length; i++) {
                if (x[i] != y[i]) {
                    return x[i] < y[i];
                }
            }
            return false;
        }
    }

    // ================================================================== MOEA/D
    public static class MOEAD extends MOEAHosts {
        public static int CONFIG_NEIGHBOURS = 6;
        public static double CONFIG_DELTA = 0.9;
        public static int CONFIG_NR = 2;
        public static long CONFIG_MAX_EVALS = -1;   // default: initial population + G x (P-1)

        @Override
        protected void evolve() {
            int N = populationSize;
            long budget = (CONFIG_MAX_EVALS > 0) ? CONFIG_MAX_EVALS : (long) N + (long) generationCount * (N - 1);
            double[][] lambda = new double[N][2];
            for (int i = 0; i < N; i++) {
                lambda[i][0] = (N == 1) ? 0.5 : (double) i / (N - 1);
                lambda[i][1] = 1.0 - lambda[i][0];
            }
            int T = Math.min(CONFIG_NEIGHBOURS, N);
            int[][] nb = new int[N][T];
            for (int i = 0; i < N; i++) {
                final int ii = i;
                Integer[] idx = new Integer[N];
                for (int j = 0; j < N; j++) {
                    idx[j] = j;
                }
                Arrays.sort(idx, (x, y) -> Double.compare(Math.abs(lambda[x][0] - lambda[ii][0]), Math.abs(lambda[y][0] - lambda[ii][0])));
                for (int t = 0; t < T; t++) {
                    nb[i][t] = idx[t];
                }
            }
            // random assignment of the initial individuals to sub-problems (standard MOEA/D)
            List<Integer> perm = new ArrayList<>();
            for (int i = 0; i < N; i++) {
                perm.add(i);
            }
            java.util.Collections.shuffle(perm, random);
            List<int[]> x = new ArrayList<>();
            double[] fm = new double[N];
            double[] fc = new double[N];
            for (int i = 0; i < N; i++) {
                x.add(population.get(perm.get(i)));
                fm[i] = makespans[perm.get(i)];
                fc[i] = costs[perm.get(i)];
            }
            double zm = Double.MAX_VALUE, zc = Double.MAX_VALUE, nm = 0, nc = 0;
            for (int i = 0; i < N; i++) {
                zm = Math.min(zm, fm[i]);
                zc = Math.min(zc, fc[i]);
                nm = Math.max(nm, fm[i]);
                nc = Math.max(nc, fc[i]);
            }
            int gen = 0;
            long nextGenEval = N + (N - 1);
            setPopulation(x, fm, fc);
            while (evalCount < budget) {
                for (int i = 0; i < N && evalCount < budget; i++) {
                    int[] pool = (random.nextDouble() < CONFIG_DELTA) ? nb[i] : null;
                    int a = (pool != null) ? pool[random.nextInt(pool.length)] : random.nextInt(N);
                    int b = (pool != null) ? pool[random.nextInt(pool.length)] : random.nextInt(N);
                    int[] child = crossoverMutate(x.get(a), x.get(b));
                    double[] mc = evaluate(child);
                    zm = Math.min(zm, mc[0]);
                    zc = Math.min(zc, mc[1]);
                    nm = Math.max(nm, mc[0]);
                    nc = Math.max(nc, mc[1]);
                    double rm = Math.max(nm - zm, 1e-12);
                    double rc = Math.max(nc - zc, 1e-12);
                    int replaced = 0;
                    List<Integer> order = new ArrayList<>();
                    for (int t = 0; t < (pool != null ? pool.length : N); t++) {
                        order.add(pool != null ? pool[t] : t);
                    }
                    java.util.Collections.shuffle(order, random);
                    for (int j : order) {
                        double gNew = Math.max(lambda[j][0] * (mc[0] - zm) / rm, lambda[j][1] * (mc[1] - zc) / rc);
                        double gOld = Math.max(lambda[j][0] * (fm[j] - zm) / rm, lambda[j][1] * (fc[j] - zc) / rc);
                        if (gNew <= gOld) {
                            x.set(j, child);
                            fm[j] = mc[0];
                            fc[j] = mc[1];
                            if (++replaced >= CONFIG_NR) {
                                break;
                            }
                        }
                    }
                    setPopulation(x, fm, fc);
                    while (evalCount >= nextGenEval && gen < generationCount) {
                        gen++;
                        recordCheckpoint(gen);
                        nextGenEval += (N - 1);
                    }
                }
            }
            setPopulation(x, fm, fc);
        }
    }
}
