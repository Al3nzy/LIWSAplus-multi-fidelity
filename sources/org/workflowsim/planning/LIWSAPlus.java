package org.workflowsim.planning;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;

/**
 * LIWSA+ : the LIWSA recombination moves inside a population-level survival and variation framework.
 *
 *  1. Survival: (mu + lambda) replacement by Pareto rank and crowding distance, so a non-dominated trade-off
 *     solution is never overwritten by an incomparable child (LIWSA's own replacement rule can do that).
 *  2. Output: the anytime archive of every non-dominated solution evaluated (reported next to the final-population front).
 *  3. Variation pool, each child produced by ONE operator:
 *       CANON : uniform crossover + random-reset mutation (canonical),
 *       LIWSA : density-adaptive solitary / social move of the original algorithm (unchanged),
 *       CRIT  : critical-path move: re-assign tasks on the chain that determines the makespan of the parent schedule
 *               (fast VM or a VM already holding a neighbour task, which also removes transfers),
 *       COST  : cost move: move tasks to a VM with a strictly lower price per unit of work.
 *     CRIT and COST use information the decoder already produces (start/finish times, VM prices); they cost no extra evaluations.
 *  4. Online learning of operator usage (adaptive operator selection, probability matching): after every generation the
 *     credit of an operator is the share of its children that entered the archive (1.0) or survived (0.5); selection
 *     probabilities follow the running credit with a floor. This is learning from the run itself, not a metaphor.
 * Every flag can be switched off for ablation (see ExperimentRunner host names LIWSAP.nog / .noa / .nol).
 */
public class LIWSAPlus extends MOEAHosts {

    public static final int OP_CANON = 1, OP_LIWSA = 2, OP_CRIT = 4, OP_COST = 8, OP_EST = 16;
    public static int CONFIG_POOL = 31;
    public static boolean CONFIG_AOS = true;
    public static double CONFIG_AOS_PMIN = 0.10;
    public static double CONFIG_AOS_ALPHA = 0.30;
    public static double CONFIG_MOVE_FRAC = 0.02;      // a guided move re-assigns between 1 and max(1, frac*n) tasks

    /** Surrogate-assisted pre-screening of candidate children (needs the faithful decoder): 0 off, 1 nominal decoder only,
     *  2 nominal decoder corrected by an online ridge model of the expensive objective, 3 random choice (null control). */
    public static int CONFIG_SA_MODE = 0;
    public static int CONFIG_SA_K = 5;                 // candidates generated per child that is actually evaluated
    public static double CONFIG_SA_RIDGE = 0.1;
    /** Selection rule: false = Pareto rank among the candidates only; true = greedy predicted hypervolume improvement over the
     *  archive of truly evaluated solutions (needs calibrated predictions, which is where learning can matter). */
    public static boolean CONFIG_SA_HVI = false;
    /** TCC-style online OLS immigrants: every CONFIG_IMM_PERIOD generations CONFIG_IMM_COUNT genotypes sampled (softmax) from a
     *  per-(task, VM) OLS model of the evaluated objectives replace children, so the evaluation budget is unchanged. */
    public static boolean CONFIG_IMM = false;
    public static int CONFIG_IMM_COUNT = 4;
    public static int CONFIG_IMM_PERIOD = 10;
    public static double CONFIG_IMM_TEMP = 0.5;

    /** Diagnostic mode: every candidate is also evaluated truly (not counted) to measure how well nominal and corrected predictions rank them. */
    public static boolean CONFIG_DIAG = false;
    public static double[] lastDiag = null;                // {rho nominal makespan, rho corrected makespan, rho nominal cost, rho corrected cost, overlap nominal, overlap corrected}
    private double[] diagSum = new double[6];
    private int diagN = 0;

    private OnlineRidge ridge;
    private OnlineRidge ridgeRich;
    private BoostedTrees[] gbm = new BoostedTrees[2];
    private final List<double[]> gbmX = new ArrayList<>();
    private final List<double[]> gbmY = new ArrayList<>();
    private double refM, refC;
    private boolean refSet;
    /** Ridge penalty of the rich-feature linear surrogate. */
    public static double CONFIG_SA_RIDGE_RICH = 1.0;
    /** Hedged selection: alternate between the nominal and the learned ranking so that the extremes of the front stay reachable. */
    public static boolean CONFIG_SA_MIX = false;
    public static int CONFIG_GBM_TREES = 60;
    public static int CONFIG_GBM_DEPTH = 3;
    private double[][][] pairFeat;
    private int[][] pairCount;
    private double[][] pairSumM, pairSumC;
    private double totM, totC;
    private int totN;
    private double[][] pairCoef;                       // [2][9]

    private int[][] childIdx;
    private double[] cpi;                              // price per unit of work of each VM (rate / mips)
    private int[][] vmBySpeed;                         // per task: VMs sorted by execution time (fastest first)
    private final Map<int[], int[]> chainCache = new IdentityHashMap<>();

    // tracked decoder buffers
    private double[] tStart, tFin;
    private int[][] evTask;
    private double[][] evS, evF;
    private int[] evN;

    @Override
    protected void evolve() {
        int P = populationSize;
        int lambda = P - 1;
        int n = nTasks;
        prepare();
        double tau = calibrateTau();
        List<Integer> ops = new ArrayList<>();
        for (int b = 0; b < 5; b++) {
            if ((CONFIG_POOL & (1 << b)) != 0) {
                ops.add(1 << b);
            }
        }
        int K = ops.size();
        double[] q = new double[K];
        double[] prob = new double[K];
        Arrays.fill(q, 1.0);
        Arrays.fill(prob, 1.0 / K);

        List<int[]> pop = new ArrayList<>(population);
        double[] pm = makespans.clone();
        double[] pc = costs.clone();
        initLearning();
        for (int i = 0; i < P; i++) {
            learn(pop.get(i), new double[]{pm[i], pc[i]}, null);
        }
        for (int gen = 0; gen < generationCount; gen++) {
            setPopulation(pop, pm, pc);
            int[] rank = paretoRanks(pm, pc);
            double[] crowd = crowding(pm, pc, rank);
            List<Integer> elite = new ArrayList<>();
            for (int i = 0; i < P; i++) {
                if (rank[i] == 0) {
                    elite.add(i);
                }
            }
            for (int r = 1; elite.size() < 3 && r < P; r++) {
                for (int i = 0; i < P; i++) {
                    if (rank[i] == r) {
                        elite.add(i);
                    }
                }
            }
            List<int[]> union = new ArrayList<>(pop);
            double[] um = Arrays.copyOf(pm, P + lambda);
            double[] uc = Arrays.copyOf(pc, P + lambda);
            int immCount = (CONFIG_IMM && gen > 0 && gen % CONFIG_IMM_PERIOD == 0) ? Math.min(CONFIG_IMM_COUNT, lambda - 1) : 0;
            int nChildren = lambda - immCount;
            int Ksa = (CONFIG_SA_MODE == 0) ? 1 : CONFIG_SA_K;
            int nCand = Ksa * nChildren;
            int[][] cand = new int[nCand][];
            int[] candOp = new int[nCand];
            for (int c = 0; c < nCand; c++) {
                int oi = pick(prob);
                int op = ops.get(oi);
                candOp[c] = oi;
                int a = tournament(rank, crowd, P);
                int b = tournament(rank, crowd, P);
                int[] child;
                if (op == OP_CANON) {
                    child = crossoverMutate(pop.get(a), pop.get(b));
                } else if (op == OP_LIWSA) {
                    double density = localDensity(a, tau);
                    double pSocial = 0.5 * ((double) gen / Math.max(generationCount, 1)) + 0.5 * density;
                    child = (random.nextDouble() > pSocial) ? solitaryMove(a, rank) : socialMove(a, rank, elite);
                    mutate(child);
                } else if (op == OP_CRIT) {
                    child = criticalMove(pop.get(a));
                } else if (op == OP_EST) {
                    child = estimateMove(pop.get(a));
                } else {
                    child = costMove(pop.get(a));
                }
                cand[c] = child;
            }
            double[][] candFeat = new double[nCand][];
            int[] chosen = new int[nChildren];
            if (Ksa == 1) {
                for (int o = 0; o < nChildren; o++) {
                    chosen[o] = o;
                }
            } else {
                chosen = selectCandidates(cand, candFeat, nChildren);
            }
            int[] opOf = new int[nChildren];
            boolean[] added = new boolean[nChildren];
            for (int o = 0; o < nChildren; o++) {
                int c = chosen[o];
                int[] child = cand[c];
                opOf[o] = candOp[c];
                double[] mc = evaluate(child);
                added[o] = lastAdded;
                union.add(child);
                um[P + o] = mc[0];
                uc[P + o] = mc[1];
                learn(child, mc, candFeat[c]);
            }
            if (immCount > 0) {
                fitPairModel();
                for (int i = 0; i < immCount; i++) {
                    int[] g = immigrant(i);
                    double[] mc = evaluate(g);
                    union.add(g);
                    um[P + nChildren + i] = mc[0];
                    uc[P + nChildren + i] = mc[1];
                    learn(g, mc, null);
                }
            }
            // (mu + lambda) survival by rank and crowding
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
            boolean[] survived = new boolean[P + lambda];
            for (int i = 0; i < P; i++) {
                np.add(union.get(order[i]));
                nm[i] = um[order[i]];
                nc[i] = uc[order[i]];
                survived[order[i]] = true;
            }
            // online update of operator probabilities
            if (CONFIG_AOS && K > 1) {
                double[] sum = new double[K];
                int[] cnt = new int[K];
                for (int o = 0; o < nChildren; o++) {
                    double r = added[o] ? 1.0 : (survived[P + o] ? 0.5 : 0.0);
                    sum[opOf[o]] += r;
                    cnt[opOf[o]]++;
                }
                double tot = 0;
                for (int k = 0; k < K; k++) {
                    if (cnt[k] > 0) {
                        q[k] = (1 - CONFIG_AOS_ALPHA) * q[k] + CONFIG_AOS_ALPHA * (sum[k] / cnt[k]);
                    }
                    tot += q[k];
                }
                for (int k = 0; k < K; k++) {
                    prob[k] = (tot > 0) ? CONFIG_AOS_PMIN + (1 - K * CONFIG_AOS_PMIN) * q[k] / tot : 1.0 / K;
                }
            }
            pop = np;
            pm = nm;
            pc = nc;
            setPopulation(pop, pm, pc);
            recordCheckpoint(gen + 1);
        }
        setPopulation(pop, pm, pc);
    }

    // ------------------------------------------------------------------ learned components
    private boolean learning() {
        return CONFIG_SA_MODE >= 2 || CONFIG_IMM;
    }

    private void initLearning() {
        ridge = new OnlineRidge(FaithfulModel.N_SURR_FEATURES, 2, CONFIG_SA_RIDGE);
        ridgeRich = new OnlineRidge(richDim(), 2, CONFIG_SA_RIDGE_RICH);
        gbmX.clear();
        gbmY.clear();
        refSet = false;
        totM = totC = 0;
        totN = 0;
        if ((CONFIG_SA_MODE > 0) && faithful == null) {
            throw new IllegalStateException("surrogate-assisted modes need the faithful decoder (.f)");
        }
        int n = nTasks;
        int m = nVms;
        pairCount = new int[n][m];
        pairSumM = new double[n][m];
        pairSumC = new double[n][m];
        pairCoef = new double[2][9];
        // 9 per-(task, VM) features of the TCC predictor
        int[] depth = new int[n];
        double maxLen = 1, maxDepth = 1, maxCh = 1, maxPa = 1, maxMips = 1, maxRate = 1, slow = Double.MAX_VALUE;
        for (int k = 0; k < n; k++) {
            int d = 0;
            for (int p : parentIdx[k]) {
                d = Math.max(d, depth[p] + 1);
            }
            depth[k] = d;
            maxDepth = Math.max(maxDepth, d);
            maxLen = Math.max(maxLen, dur[k][0] * vmList.get(0).getMips());
            maxCh = Math.max(maxCh, childIdx[k].length);
            maxPa = Math.max(maxPa, parentIdx[k].length);
        }
        for (int v = 0; v < m; v++) {
            maxMips = Math.max(maxMips, vmList.get(v).getMips());
            maxRate = Math.max(maxRate, vmRate[v]);
            slow = Math.min(slow, vmList.get(v).getMips());
        }
        double maxDur = maxLen / Math.max(slow, 1.0);
        double maxPc = Math.max(maxDur * maxRate, 1e-9);
        pairFeat = new double[n][m][];
        for (int k = 0; k < n; k++) {
            double len = dur[k][0] * vmList.get(0).getMips();
            for (int v = 0; v < m; v++) {
                double d = dur[k][v];
                pairFeat[k][v] = new double[]{len / maxLen, depth[k] / maxDepth, childIdx[k].length / maxCh, parentIdx[k].length / maxPa,
                    vmList.get(v).getMips() / maxMips, vmRate[v] / maxRate, d / maxDur, d * vmRate[v] / maxPc, 1.0};
            }
        }
    }

    /** Adds one expensively evaluated schedule to the surrogate and to the pair-level OLS statistics. */
    private void learn(int[] g, double[] mc, double[] feat) {
        if (!learning()) {
            return;
        }
        if (CONFIG_SA_MODE >= 2 && CONFIG_SA_MODE != 3) {
            double[] f = feat;
            if (f == null) {
                f = faithful.features(g);
                nominalDecodes++;
            }
            double[] y = {Math.log(Math.max(mc[0], 1e-9) / Math.max(f[0], 1e-9)), Math.log(Math.max(mc[1], 1e-9) / Math.max(f[1], 1e-9))};
            if (CONFIG_SA_MODE == 2) {
                ridge.add(Arrays.copyOfRange(f, 2, 2 + FaithfulModel.N_SURR_FEATURES), y);
            } else if (CONFIG_SA_MODE == 4) {
                ridgeRich.add(richX(f), y);
            } else if (CONFIG_SA_MODE == 5) {
                gbmX.add(richX(f));
                gbmY.add(y);
            }
        }
        if (CONFIG_IMM) {
            for (int k = 0; k < g.length; k++) {
                pairCount[k][g[k]]++;
                pairSumM[k][g[k]] += mc[0];
                pairSumC[k][g[k]] += mc[1];
            }
            totM += mc[0];
            totC += mc[1];
            totN++;
        }
    }

    // ------------------------------------------------------------------ richer surrogates (modes 4 and 5)
    private int richDim() {
        return 1 + 9 + 2 + nVms + 5;
    }

    /**
     * Input of the rich surrogates: structural descriptors, log nominal objectives relative to the first schedule, the busy-time share
     * of every machine, and five scenario-exposure features computed from the machine loads and the speed multipliers of the fixed
     * search scenarios (bottleneck-inflated load: mean, tail mean, spread; inflated cost: mean, tail mean).
     */
    private double[] richX(double[] f) {
        if (!refSet) {
            refM = Math.max(f[0], 1e-9);
            refC = Math.max(f[1], 1e-9);
            refSet = true;
        }
        int m = nVms;
        double ms = Math.max(f[0], 1e-9);
        double[] x = new double[richDim()];
        int p = 0;
        x[p++] = 1.0;
        for (int i = 3; i < 12; i++) {
            x[p++] = f[i];
        }
        x[p++] = Math.log(ms / refM);
        x[p++] = Math.log(Math.max(f[1], 1e-9) / refC);
        for (int v = 0; v < m; v++) {
            x[p++] = f[12 + v];
        }
        int S = robMips == null ? 0 : robMips.length;
        if (S > 0) {
            double[] b = new double[S], cc = new double[S];
            double base = 0;
            for (int v = 0; v < m; v++) {
                base += vmRate[v] * f[12 + v] * ms;
            }
            base = Math.max(base, 1e-9);
            for (int s = 0; s < S; s++) {
                double mx = 0, cs = 0;
                for (int v = 0; v < m; v++) {
                    double inf = f[12 + v] * ms / robMips[s][v];
                    mx = Math.max(mx, inf);
                    cs += vmRate[v] * inf;
                }
                b[s] = mx / ms;
                cc[s] = cs / base;
            }
            int tail = Math.max(1, (int) Math.ceil((1.0 - (CONFIG_ROBUST_CVAR > 0 ? CONFIG_ROBUST_CVAR : 0.8)) * S));
            double[] sb = b.clone(), sc = cc.clone();
            Arrays.sort(sb);
            Arrays.sort(sc);
            double mb = 0, tb = 0, mc = 0, tc = 0;
            for (int s = 0; s < S; s++) {
                mb += b[s];
                mc += cc[s];
            }
            for (int s = S - tail; s < S; s++) {
                tb += sb[s];
                tc += sc[s];
            }
            mb /= S;
            mc /= S;
            double var = 0;
            for (int s = 0; s < S; s++) {
                var += (b[s] - mb) * (b[s] - mb);
            }
            x[p++] = mb;
            x[p++] = tb / tail;
            x[p++] = Math.sqrt(var / S);
            x[p++] = mc;
            x[p++] = tc / tail;
        }
        return x;
    }

    /** Refits the active surrogate on everything evaluated so far; false while there is too little data (or for the nominal and random modes). */
    private boolean fitModel() {
        if (CONFIG_SA_MODE == 2 && ridge.size() >= 12) {
            ridge.fit();
            return true;
        }
        if (CONFIG_SA_MODE == 4 && ridgeRich.size() >= 12) {
            ridgeRich.fit();
            return true;
        }
        if (CONFIG_SA_MODE == 5 && gbmX.size() >= 12) {
            double[][] xs = gbmX.toArray(new double[0][]);
            for (int t = 0; t < 2; t++) {
                double[] y = new double[xs.length];
                for (int i = 0; i < y.length; i++) {
                    y[i] = gbmY.get(i)[t];
                }
                gbm[t] = new BoostedTrees(CONFIG_GBM_TREES, CONFIG_GBM_DEPTH, 0.1, 5, 16, 0.8, 1234L + t);
                gbm[t].fit(xs, y);
            }
            return true;
        }
        return false;
    }

    /** Predicted log ratio of the expensive to the nominal objectives for the schedule with nominal features f. */
    private double[] corr(double[] f) {
        if (CONFIG_SA_MODE == 2) {
            double[] x = Arrays.copyOfRange(f, 2, 2 + FaithfulModel.N_SURR_FEATURES);
            return new double[]{ridge.predict(0, x), ridge.predict(1, x)};
        }
        double[] x = richX(f);
        if (CONFIG_SA_MODE == 4) {
            return new double[]{ridgeRich.predict(0, x), ridgeRich.predict(1, x)};
        }
        return new double[]{gbm[0].predict(x), gbm[1].predict(x)};
    }

    /** Chooses nChildren of the candidates by predicted Pareto rank and crowding (candidates only); computes nominal features. */
    private int[] selectCandidates(int[][] cand, double[][] candFeat, int nChildren) {
        int nc = cand.length;
        double[] pm = new double[nc];
        double[] pc = new double[nc];
        if (CONFIG_SA_MODE == 3) {
            Integer[] idx = new Integer[nc];
            for (int i = 0; i < nc; i++) {
                idx[i] = i;
            }
            java.util.Collections.shuffle(Arrays.asList(idx), random);
            int[] out = new int[nChildren];
            for (int i = 0; i < nChildren; i++) {
                out[i] = idx[i];
            }
            return out;
        }
        boolean useModel = fitModel();
        for (int c = 0; c < nc; c++) {
            double[] f = faithful.features(cand[c]);
            nominalDecodes++;
            candFeat[c] = f;
            pm[c] = f[0];
            pc[c] = f[1];
            if (useModel) {
                double[] lr = corr(f);
                pm[c] *= Math.exp(lr[0]);
                pc[c] *= Math.exp(lr[1]);
            }
        }
        if (CONFIG_SA_HVI && !archiveFront.isEmpty()) {
            return selectByHvi(pm, pc, nChildren);
        }
        if (CONFIG_SA_MIX && useModel) {
            double[] nm0 = new double[nc], nc0 = new double[nc];
            for (int c = 0; c < nc; c++) {
                nm0[c] = candFeat[c][0];
                nc0[c] = candFeat[c][1];
            }
            int[] byNom = topByPareto(nm0, nc0, nc);
            int[] byMod = topByPareto(pm, pc, nc);
            boolean[] taken = new boolean[nc];
            int[] out = new int[nChildren];
            int got = 0, a = 0, b = 0;
            while (got < nChildren) {
                while (a < nc && taken[byNom[a]]) {
                    a++;
                }
                if (got < nChildren && a < nc) {
                    taken[byNom[a]] = true;
                    out[got++] = byNom[a];
                }
                while (b < nc && taken[byMod[b]]) {
                    b++;
                }
                if (got < nChildren && b < nc) {
                    taken[byMod[b]] = true;
                    out[got++] = byMod[b];
                }
            }
            return out;
        }
        if (CONFIG_DIAG) {
            diagnose(cand, candFeat, nChildren);
        }
        int[] rk = paretoRanks(pm, pc);
        double[] cw = crowding(pm, pc, rk);
        Integer[] order = new Integer[nc];
        for (int i = 0; i < nc; i++) {
            order[i] = i;
        }
        Arrays.sort(order, (x, y) -> rk[x] != rk[y] ? Integer.compare(rk[x], rk[y]) : Double.compare(cw[y], cw[x]));
        int[] out = new int[nChildren];
        for (int i = 0; i < nChildren; i++) {
            out[i] = order[i];
        }
        return out;
    }

    private int[] topByPareto(double[] a, double[] b, int k) {
        int[] rk = paretoRanks(a, b);
        double[] cw = crowding(a, b, rk);
        Integer[] order = new Integer[a.length];
        for (int i = 0; i < a.length; i++) {
            order[i] = i;
        }
        Arrays.sort(order, (x, y) -> rk[x] != rk[y] ? Integer.compare(rk[x], rk[y]) : Double.compare(cw[y], cw[x]));
        int[] out = new int[k];
        for (int i = 0; i < k; i++) {
            out[i] = order[i];
        }
        return out;
    }

    private void diagnose(int[][] cand, double[][] candFeat, int nChildren) {
        int nc = cand.length;
        double[] tm = new double[nc], tc = new double[nc], nm = new double[nc], nc2 = new double[nc], cm = new double[nc], cc = new double[nc];
        boolean model = fitModel();
        for (int c = 0; c < nc; c++) {
            double[] t = decode(cand[c]);
            tm[c] = t[0];
            tc[c] = t[1];
            nm[c] = candFeat[c][0];
            nc2[c] = candFeat[c][1];
            cm[c] = nm[c];
            cc[c] = nc2[c];
            if (model) {
                double[] lr = corr(candFeat[c]);
                cm[c] *= Math.exp(lr[0]);
                cc[c] *= Math.exp(lr[1]);
            }
        }
        int[] truth = topByPareto(tm, tc, nChildren);
        java.util.Set<Integer> ts = new java.util.HashSet<>();
        for (int i : truth) {
            ts.add(i);
        }
        int[] sn = topByPareto(nm, nc2, nChildren);
        int[] sc = topByPareto(cm, cc, nChildren);
        int on = 0, oc = 0;
        for (int i : sn) {
            if (ts.contains(i)) {
                on++;
            }
        }
        for (int i : sc) {
            if (ts.contains(i)) {
                oc++;
            }
        }
        double[] v = {spearman(nm, tm), spearman(cm, tm), spearman(nc2, tc), spearman(cc, tc), on / (double) nChildren, oc / (double) nChildren};
        for (int i = 0; i < 6; i++) {
            if (!Double.isNaN(v[i])) {
                diagSum[i] += v[i];
            }
        }
        diagN++;
        lastDiag = new double[6];
        for (int i = 0; i < 6; i++) {
            lastDiag[i] = diagSum[i] / diagN;
        }
    }

    private static double hv(List<double[]> pts, double r0, double r1) {
        List<double[]> a = new ArrayList<>();
        for (double[] p : pts) {
            if (p[0] < r0 && p[1] < r1) {
                a.add(p);
            }
        }
        a.sort((x, y) -> Double.compare(x[0], y[0]));
        double v = 0, best = r1;
        for (double[] p : a) {
            if (p[1] < best) {
                v += (r0 - p[0]) * (best - p[1]);
                best = p[1];
            }
        }
        return v;
    }

    /** Greedy selection of the candidates with the largest predicted hypervolume improvement over the evaluated archive. */
    private int[] selectByHvi(double[] pm, double[] pc, int nChildren) {
        int nc = pm.length;
        double r0 = 0, r1 = 0;
        List<double[]> front = new ArrayList<>(archiveFront);
        for (double[] p : front) {
            r0 = Math.max(r0, p[0]);
            r1 = Math.max(r1, p[1]);
        }
        r0 *= 1.2;
        r1 *= 1.2;
        boolean[] used = new boolean[nc];
        int[] out = new int[nChildren];
        double cur = hv(front, r0, r1);
        // fallback ordering for rounds where no candidate improves the hypervolume: distance-like score (sum of normalised excess)
        for (int i = 0; i < nChildren; i++) {
            int best = -1;
            double bestGain = 0;
            for (int c = 0; c < nc; c++) {
                if (used[c]) {
                    continue;
                }
                List<double[]> t = new ArrayList<>(front);
                t.add(new double[]{pm[c], pc[c]});
                double gain = hv(t, r0, r1) - cur;
                if (gain > bestGain) {
                    bestGain = gain;
                    best = c;
                }
            }
            if (best < 0) {
                double bs = Double.MAX_VALUE;
                for (int c = 0; c < nc; c++) {
                    if (used[c]) {
                        continue;
                    }
                    double sc = 0;
                    for (double[] p : front) {
                        sc = Math.max(sc, 0);
                    }
                    double d = Double.MAX_VALUE;
                    for (double[] p : front) {
                        double dx = Math.max(pm[c] - p[0], 0) / r0, dy = Math.max(pc[c] - p[1], 0) / r1;
                        d = Math.min(d, Math.max(dx, dy));
                    }
                    if (d < bs) {
                        bs = d;
                        best = c;
                    }
                }
            }
            used[best] = true;
            out[i] = best;
            front.add(new double[]{pm[best], pc[best]});
            cur = hv(front, r0, r1);
        }
        return out;
    }

    /** Refits the two per-(task, VM) OLS models from all evaluated schedules (normal equations from sufficient statistics). */
    private void fitPairModel() {
        if (totN == 0) {
            return;
        }
        int p = 9;
        double[][] xtx = new double[p][p];
        double[][] xty = new double[2][p];
        double mm = totM / totN;
        double mc = totC / totN;
        for (int k = 0; k < nTasks; k++) {
            for (int v = 0; v < nVms; v++) {
                int c = pairCount[k][v];
                if (c == 0) {
                    continue;
                }
                double[] x = pairFeat[k][v];
                for (int i = 0; i < p; i++) {
                    for (int j = 0; j < p; j++) {
                        xtx[i][j] += c * x[i] * x[j];
                    }
                    xty[0][i] += x[i] * pairSumM[k][v] / mm;
                    xty[1][i] += x[i] * pairSumC[k][v] / mc;
                }
            }
        }
        for (int t = 0; t < 2; t++) {
            double[][] aug = new double[p][p + 1];
            for (int i = 0; i < p; i++) {
                for (int j = 0; j < p; j++) {
                    aug[i][j] = xtx[i][j] + (i == j ? 1e-9 : 0.0);
                }
                aug[i][p] = xty[t][i];
            }
            pairCoef[t] = OnlineRidge.solve(aug, p);
        }
    }

    private static final double[] IMM_WEIGHTS = {0.9, 0.7, 0.5, 0.3};

    /** One immigrant: per task, a softmax sample over machines of the negative predicted (weighted, mean-normalised) objective. */
    private int[] immigrant(int idx) {
        double wM = IMM_WEIGHTS[idx % IMM_WEIGHTS.length];
        int[] g = new int[nTasks];
        double[] sc = new double[nVms];
        for (int k = 0; k < nTasks; k++) {
            double best = -Double.MAX_VALUE;
            for (int v = 0; v < nVms; v++) {
                double[] x = pairFeat[k][v];
                double pmv = 0, pcv = 0;
                for (int i = 0; i < 9; i++) {
                    pmv += pairCoef[0][i] * x[i];
                    pcv += pairCoef[1][i] * x[i];
                }
                sc[v] = -(wM * pmv + (1 - wM) * pcv);
                best = Math.max(best, sc[v]);
            }
            double sum = 0;
            double[] pr = new double[nVms];
            for (int v = 0; v < nVms; v++) {
                pr[v] = Math.exp((sc[v] - best) / Math.max(CONFIG_IMM_TEMP, 1e-6));
                sum += pr[v];
            }
            double r = random.nextDouble() * sum;
            double cum = 0;
            int pick = nVms - 1;
            for (int v = 0; v < nVms; v++) {
                cum += pr[v];
                if (r <= cum) {
                    pick = v;
                    break;
                }
            }
            g[k] = pick;
        }
        return g;
    }

    // ------------------------------------------------------------------ helpers
    private void prepare() {
        int n = nTasks;
        int m = nVms;
        List<List<Integer>> kids = new ArrayList<>();
        for (int k = 0; k < n; k++) {
            kids.add(new ArrayList<Integer>());
        }
        for (int k = 0; k < n; k++) {
            for (int p : parentIdx[k]) {
                kids.get(p).add(k);
            }
        }
        childIdx = new int[n][];
        for (int k = 0; k < n; k++) {
            childIdx[k] = new int[kids.get(k).size()];
            for (int i = 0; i < childIdx[k].length; i++) {
                childIdx[k][i] = kids.get(k).get(i);
            }
        }
        cpi = new double[m];
        for (int v = 0; v < m; v++) {
            cpi[v] = vmRate[v] / vmList.get(v).getMips();
        }
        vmBySpeed = new int[n][m];
        for (int k = 0; k < n; k++) {
            final int kk = k;
            Integer[] idx = new Integer[m];
            for (int v = 0; v < m; v++) {
                idx[v] = v;
            }
            Arrays.sort(idx, (x, y) -> Double.compare(dur[kk][x], dur[kk][y]));
            for (int v = 0; v < m; v++) {
                vmBySpeed[k][v] = idx[v];
            }
        }
        tStart = new double[n];
        tFin = new double[n];
        evTask = new int[m][n + 1];
        evS = new double[m][n + 1];
        evF = new double[m][n + 1];
        evN = new int[m];
    }

    protected int pick(double[] prob) {
        double r = random.nextDouble();
        double acc = 0;
        for (int k = 0; k < prob.length; k++) {
            acc += prob[k];
            if (r < acc) {
                return k;
            }
        }
        return prob.length - 1;
    }

    protected int tournament(int[] rank, double[] crowd, int n) {
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

    protected double[] crowding(double[] m, double[] c, int[] rank) {
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
                for (int i = 1; i + 1 < idx.size(); i++) {
                    d[idx.get(i)] += (f[idx.get(i + 1)] - f[idx.get(i - 1)]) / range;
                }
            }
        }
        return d;
    }

    protected int numChanges() {
        int cap = Math.max(1, (int) (CONFIG_MOVE_FRAC * nTasks));
        return 1 + random.nextInt(cap);
    }

    /** Re-assign tasks on the makespan-determining chain of the parent schedule. */
    private int[] criticalMove(int[] parent) {
        int[] chain = chainOf(parent);
        if (chain.length == 0) {
            return crossoverMutate(parent, parent);
        }
        int[] child = parent.clone();
        int changes = numChanges();
        int m = nVms;
        for (int t = 0; t < changes; t++) {
            int k = chain[random.nextInt(chain.length)];
            int cur = child[k];
            int v = cur;
            if (random.nextBoolean()) {
                int top = Math.max(2, m / 4);
                v = vmBySpeed[k][random.nextInt(top)];
            } else {
                int np = parentIdx[k].length;
                int nc = childIdx[k].length;
                if (np + nc > 0) {
                    int pick = random.nextInt(np + nc);
                    v = child[pick < np ? parentIdx[k][pick] : childIdx[k][pick - np]];
                }
            }
            if (v == cur) {
                v = random.nextInt(m);
            }
            child[k] = v;
        }
        return child;
    }

    /**
     * Estimate-guided move (needs the faithful model): for tasks on the critical chain or at random, choose the VM that minimises a
     * random time/cost weighting of the estimated occupied time (execution plus fetching inputs produced on other VMs); the best VM
     * is taken with probability 0.7, otherwise one of the three best, so repeated moves explore different trade-offs.
     */
    private int[] estimateMove(int[] parent) {
        if (faithful == null) {
            return costMove(parent);
        }
        int[] chain = chainOf(parent);
        int[] child = parent.clone();
        int changes = numChanges();
        int m = nVms;
        double[] est = new double[m];
        double[] score = new double[m];
        for (int t = 0; t < changes; t++) {
            int k = (chain.length > 0 && random.nextBoolean()) ? chain[random.nextInt(chain.length)] : random.nextInt(nTasks);
            double w = random.nextDouble();
            double minE = Double.MAX_VALUE, minC = Double.MAX_VALUE;
            for (int v = 0; v < m; v++) {
                est[v] = faithful.estimate(k, v, child);
                minE = Math.min(minE, est[v]);
                minC = Math.min(minC, est[v] * vmRate[v]);
            }
            for (int v = 0; v < m; v++) {
                score[v] = w * est[v] / minE + (1 - w) * est[v] * vmRate[v] / minC;
            }
            Integer[] idx = new Integer[m];
            for (int v = 0; v < m; v++) {
                idx[v] = v;
            }
            final double[] sc = score;
            Arrays.sort(idx, (x, y) -> Double.compare(sc[x], sc[y]));
            child[k] = (random.nextDouble() < 0.7) ? idx[0] : idx[random.nextInt(Math.min(3, m))];
        }
        return child;
    }

    /** Move tasks to a VM with a strictly lower price per unit of work. */
    private int[] costMove(int[] parent) {
        int[] child = parent.clone();
        int changes = numChanges();
        int m = nVms;
        int moved = 0;
        for (int tries = 0; tries < 8 * changes && moved < changes; tries++) {
            int k = random.nextInt(nTasks);
            int cur = child[k];
            int[] better = new int[m];
            int nb = 0;
            for (int v = 0; v < m; v++) {
                if (cpi[v] < cpi[cur] - 1e-12) {
                    better[nb++] = v;
                }
            }
            if (nb == 0) {
                continue;
            }
            child[k] = better[random.nextInt(nb)];
            moved++;
        }
        if (moved == 0) {
            return crossoverMutate(parent, parent);
        }
        return child;
    }

    /** Chain of tasks that binds the makespan (data dependency or VM queue), by tracing binding predecessors. */
    protected int[] chainOf(int[] g) {
        int[] cached = chainCache.get(g);
        if (cached != null) {
            return cached;
        }
        if (faithful != null) {
            faithful.decode(g);
            int[] fc = faithful.criticalChain();
            chainCache.put(g, fc);
            return fc;
        }
        int n = nTasks;
        Arrays.fill(evN, 0);
        for (int k = 0; k < n; k++) {
            int v = g[k];
            double ready = 0;
            for (int q = 0; q < parentIdx[k].length; q++) {
                int pk = parentIdx[k][q];
                double pf = tFin[pk] + (g[pk] != v ? parentTc[k][q] : 0.0);
                if (pf > ready) {
                    ready = pf;
                }
            }
            double d = dur[k][v];
            int cnt = evN[v];
            double start = Math.max(ready, cnt > 0 ? evF[v][cnt - 1] : 0.0);
            int pos = cnt;
            // identical gap rule family: earliest feasible gap
            double t = ready;
            for (int i = 0; i < cnt; i++) {
                if (t + d <= evS[v][i]) {
                    pos = i;
                    start = t;
                    break;
                }
                t = Math.max(t, evF[v][i]);
                if (i == cnt - 1) {
                    start = t;
                }
            }
            for (int i = cnt; i > pos; i--) {
                evS[v][i] = evS[v][i - 1];
                evF[v][i] = evF[v][i - 1];
                evTask[v][i] = evTask[v][i - 1];
            }
            evS[v][pos] = start;
            evF[v][pos] = start + d;
            evTask[v][pos] = k;
            evN[v] = cnt + 1;
            tStart[k] = start;
            tFin[k] = start + d;
        }
        int cur = 0;
        for (int k = 1; k < n; k++) {
            if (tFin[k] > tFin[cur]) {
                cur = k;
            }
        }
        List<Integer> chain = new ArrayList<>();
        chain.add(cur);
        for (int guard = 0; guard < n; guard++) {
            int pred = -1;
            double bestArr = -1;
            for (int q = 0; q < parentIdx[cur].length; q++) {
                int pk = parentIdx[cur][q];
                double arr = tFin[pk] + (g[pk] != g[cur] ? parentTc[cur][q] : 0.0);
                if (arr >= tStart[cur] - 1e-9 && arr > bestArr) {
                    bestArr = arr;
                    pred = pk;
                }
            }
            if (pred < 0) {
                int v = g[cur];
                for (int i = 1; i < evN[v]; i++) {
                    if (evTask[v][i] == cur) {
                        int qk = evTask[v][i - 1];
                        if (tFin[qk] >= tStart[cur] - 1e-9) {
                            pred = qk;
                        }
                        break;
                    }
                }
            }
            if (pred < 0 || tStart[cur] <= 1e-12) {
                break;
            }
            chain.add(pred);
            cur = pred;
        }
        int[] out = new int[chain.size()];
        for (int i = 0; i < out.length; i++) {
            out[i] = chain.get(i);
        }
        chainCache.put(g, out);
        return out;
    }
}
