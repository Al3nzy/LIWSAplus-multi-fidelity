package org.workflowsim.planning;

import java.util.ArrayList;
import java.util.List;

/**
 * LIWSA host that injects, next to the HEFT and Min-Min schedules, a sweep of weighted list-scheduling seeds
 * (the 'plus' start). Every algorithm of the study that uses this start receives the same seeds, so they cannot
 * favour one algorithm. The seeds use only the problem data and cost no search evaluations.
 */
public class LIWSASeededPlanningAlgorithm extends LIWSAPlanningAlgorithm {

    /** Number of weighted list-scheduling (HEFT-style, time/cost trade-off sweep) seeds injected; 0 = none. */
    public static int CONFIG_SWEEP_SEEDS = 0;

    @Override
    protected List<int[]> generateSeedGenotypes() {
        List<int[]> seeds = super.generateSeedGenotypes();
        if (CONFIG_SWEEP_SEEDS > 0) {
            seeds.addAll(sweepSeeds(CONFIG_SWEEP_SEEDS));
        }
        return seeds;
    }


    // ------------------------------------------------------------------ weighted list-scheduling seeds
    /**
     * Deterministic constructive seeds: tasks are taken in HEFT upward-rank order and each is placed on the VM that
     * minimises  alpha * normalisedEarliestFinish + (1 - alpha) * normalisedCost  (insertion-based, so idle gaps are
     * reused), with alpha swept from 1 (pure HEFT-style) to 0 (cost only). Uses only the problem data, no search
     * evaluations. Every host receives the same seeds in the 'plus' start, so they cannot favour one algorithm.
     */
    protected List<int[]> sweepSeeds(int count) {
        int n = nTasks;
        int m = nVms;
        List<List<int[]>> kids = new ArrayList<>();
        for (int k = 0; k < n; k++) {
            kids.add(new ArrayList<int[]>());
        }
        double[] avgDur = new double[n];
        for (int k = 0; k < n; k++) {
            double sum = 0;
            for (int v = 0; v < m; v++) {
                sum += dur[k][v];
            }
            avgDur[k] = sum / m;
            for (int q = 0; q < parentIdx[k].length; q++) {
                kids.get(parentIdx[k][q]).add(new int[]{k, q});
            }
        }
        double[] rank = new double[n];
        for (int k = n - 1; k >= 0; k--) {
            double best = 0;
            for (int[] c : kids.get(k)) {
                best = Math.max(best, parentTc[c[0]][c[1]] + rank[c[0]]);
            }
            rank[k] = avgDur[k] + best;
        }
        Integer[] order = new Integer[n];
        for (int k = 0; k < n; k++) {
            order[k] = k;
        }
        java.util.Arrays.sort(order, (x, y) -> rank[x] != rank[y] ? Double.compare(rank[y], rank[x]) : Integer.compare(x, y));

        List<int[]> out = new ArrayList<>();
        for (int s = 0; s < count; s++) {
            double alpha = (count == 1) ? 1.0 : 1.0 - (double) s / (count - 1);
            double[][] evS = new double[m][n + 1];
            double[][] evF = new double[m][n + 1];
            int[] cnt = new int[m];
            double[] fin = new double[n];
            int[] g = new int[n];
            double[] eft = new double[m];
            double[] cst = new double[m];
            double[] st = new double[m];
            int[] at = new int[m];
            for (int idx = 0; idx < n; idx++) {
                int k = order[idx];
                double minE = Double.MAX_VALUE, maxE = -Double.MAX_VALUE, minC = Double.MAX_VALUE, maxC = -Double.MAX_VALUE;
                for (int v = 0; v < m; v++) {
                    double ready = 0;
                    for (int q = 0; q < parentIdx[k].length; q++) {
                        int pk = parentIdx[k][q];
                        double arr = fin[pk] + (g[pk] != v ? parentTc[k][q] : 0.0);
                        ready = Math.max(ready, arr);
                    }
                    double t = ready;
                    int pos = cnt[v];
                    for (int i = 0; i < cnt[v]; i++) {
                        if (t + dur[k][v] <= evS[v][i]) {
                            pos = i;
                            break;
                        }
                        t = Math.max(t, evF[v][i]);
                    }
                    st[v] = t;
                    at[v] = pos;
                    eft[v] = t + dur[k][v];
                    cst[v] = dur[k][v] * vmRate[v];
                    minE = Math.min(minE, eft[v]);
                    maxE = Math.max(maxE, eft[v]);
                    minC = Math.min(minC, cst[v]);
                    maxC = Math.max(maxC, cst[v]);
                }
                int bestV = 0;
                double bestScore = Double.MAX_VALUE;
                for (int v = 0; v < m; v++) {
                    double e = (eft[v] - minE) / (maxE - minE + 1e-12);
                    double c = (cst[v] - minC) / (maxC - minC + 1e-12);
                    double score = alpha * e + (1 - alpha) * c;
                    if (score < bestScore - 1e-12) {
                        bestScore = score;
                        bestV = v;
                    }
                }
                g[k] = bestV;
                fin[k] = eft[bestV];
                int pos = at[bestV];
                for (int i = cnt[bestV]; i > pos; i--) {
                    evS[bestV][i] = evS[bestV][i - 1];
                    evF[bestV][i] = evF[bestV][i - 1];
                }
                evS[bestV][pos] = st[bestV];
                evF[bestV][pos] = eft[bestV];
                cnt[bestV]++;
            }
            out.add(g);
        }
        return out;
    }

    static double spearman(double[] a, double[] b) {
        double[] ra = ranks(a), rb = ranks(b);
        double ma = 0, mb = 0;
        for (int i = 0; i < ra.length; i++) {
            ma += ra[i];
            mb += rb[i];
        }
        ma /= ra.length;
        mb /= rb.length;
        double num = 0, da = 0, db = 0;
        for (int i = 0; i < ra.length; i++) {
            num += (ra[i] - ma) * (rb[i] - mb);
            da += (ra[i] - ma) * (ra[i] - ma);
            db += (rb[i] - mb) * (rb[i] - mb);
        }
        return (da == 0 || db == 0) ? Double.NaN : num / Math.sqrt(da * db);
    }

    private static double[] ranks(double[] v) {
        Integer[] idx = new Integer[v.length];
        for (int i = 0; i < v.length; i++) {
            idx[i] = i;
        }
        java.util.Arrays.sort(idx, (x, y) -> Double.compare(v[x], v[y]));
        double[] r = new double[v.length];
        int i = 0;
        while (i < idx.length) {
            int j = i;
            while (j + 1 < idx.length && v[idx[j + 1]] == v[idx[i]]) {
                j++;
            }
            double avg = (i + j) / 2.0 + 1.0;
            for (int q = i; q <= j; q++) {
                r[idx[q]] = avg;
            }
            i = j + 1;
        }
        return r;
    }
}
