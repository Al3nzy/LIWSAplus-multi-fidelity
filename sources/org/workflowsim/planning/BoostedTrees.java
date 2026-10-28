package org.workflowsim.planning;

import java.util.Arrays;
import java.util.Random;

/**
 * Histogram gradient-boosted regression trees (squared loss), refitted from scratch on all evaluated schedules.
 * Small and dependency free: shallow trees, quantile bins, row subsampling with a fixed seed (deterministic).
 */
public class BoostedTrees {
    private final int nTrees, depth, minLeaf, nBins;
    private final double lr, subsample;
    private final long seed;
    private double base;
    private int[][] feat;        // [tree][node]
    private double[][] thr, val; // split threshold / leaf value
    private int[][] left, right;

    public BoostedTrees(int nTrees, int depth, double lr, int minLeaf, int nBins, double subsample, long seed) {
        this.nTrees = nTrees;
        this.depth = depth;
        this.lr = lr;
        this.minLeaf = minLeaf;
        this.nBins = nBins;
        this.subsample = subsample;
        this.seed = seed;
    }

    public void fit(double[][] x, double[] y) {
        int n = y.length, d = x[0].length;
        double[][] edges = new double[d][];
        int[][] bin = new int[n][d];
        for (int j = 0; j < d; j++) {
            double[] col = new double[n];
            for (int i = 0; i < n; i++) {
                col[i] = x[i][j];
            }
            Arrays.sort(col);
            double[] e = new double[nBins - 1];
            int ne = 0;
            for (int b = 1; b < nBins; b++) {
                double q = col[Math.min(n - 1, (int) ((long) b * n / nBins))];
                if (ne == 0 || q > e[ne - 1]) {
                    e[ne++] = q;
                }
            }
            edges[j] = Arrays.copyOf(e, ne);
            for (int i = 0; i < n; i++) {
                int b = 0;
                while (b < ne && x[i][j] > edges[j][b]) {
                    b++;
                }
                bin[i][j] = b;
            }
        }
        double mean = 0;
        for (double v : y) {
            mean += v;
        }
        base = mean / n;
        double[] pred = new double[n];
        Arrays.fill(pred, base);
        feat = new int[nTrees][];
        thr = new double[nTrees][];
        val = new double[nTrees][];
        left = new int[nTrees][];
        right = new int[nTrees][];
        Random rnd = new Random(seed + n);
        int maxNodes = (1 << (depth + 1)) - 1;
        for (int t = 0; t < nTrees; t++) {
            double[] r = new double[n];
            for (int i = 0; i < n; i++) {
                r[i] = y[i] - pred[i];
            }
            int[] rows = new int[n];
            int nr = 0;
            for (int i = 0; i < n; i++) {
                if (rnd.nextDouble() < subsample) {
                    rows[nr++] = i;
                }
            }
            if (nr < 2 * minLeaf) {
                nr = n;
                for (int i = 0; i < n; i++) {
                    rows[i] = i;
                }
            }
            feat[t] = new int[maxNodes];
            thr[t] = new double[maxNodes];
            val[t] = new double[maxNodes];
            left[t] = new int[maxNodes];
            right[t] = new int[maxNodes];
            Arrays.fill(feat[t], -1);
            int[] cnt = {1};
            grow(t, 0, Arrays.copyOf(rows, nr), 0, r, bin, edges, cnt);
            for (int i = 0; i < n; i++) {
                pred[i] += predictTree(t, x[i]);
            }
        }
    }

    private void grow(int t, int node, int[] idx, int level, double[] r, int[][] bin, double[][] edges, int[] cnt) {
        double sum = 0;
        for (int i : idx) {
            sum += r[i];
        }
        val[t][node] = lr * sum / Math.max(1, idx.length);
        if (level >= depth || idx.length < 2 * minLeaf) {
            return;
        }
        int d = edges.length;
        double bestGain = 1e-12;
        int bf = -1, bb = -1;
        for (int j = 0; j < d; j++) {
            int nb = edges[j].length + 1;
            if (nb < 2) {
                continue;
            }
            double[] s = new double[nb];
            int[] c = new int[nb];
            for (int i : idx) {
                s[bin[i][j]] += r[i];
                c[bin[i][j]]++;
            }
            double sl = 0;
            int cl = 0;
            for (int b = 0; b < nb - 1; b++) {
                sl += s[b];
                cl += c[b];
                int cr = idx.length - cl;
                if (cl < minLeaf || cr < minLeaf) {
                    continue;
                }
                double sr = sum - sl;
                double gain = sl * sl / cl + sr * sr / cr - sum * sum / idx.length;
                if (gain > bestGain) {
                    bestGain = gain;
                    bf = j;
                    bb = b;
                }
            }
        }
        if (bf < 0) {
            return;
        }
        int[] li = new int[idx.length], ri = new int[idx.length];
        int nl = 0, nr2 = 0;
        for (int i : idx) {
            if (bin[i][bf] <= bb) {
                li[nl++] = i;
            } else {
                ri[nr2++] = i;
            }
        }
        feat[t][node] = bf;
        thr[t][node] = edges[bf][bb];
        left[t][node] = cnt[0]++;
        right[t][node] = cnt[0]++;
        grow(t, left[t][node], Arrays.copyOf(li, nl), level + 1, r, bin, edges, cnt);
        grow(t, right[t][node], Arrays.copyOf(ri, nr2), level + 1, r, bin, edges, cnt);
    }

    private double predictTree(int t, double[] x) {
        int node = 0;
        while (feat[t][node] >= 0) {
            node = x[feat[t][node]] <= thr[t][node] ? left[t][node] : right[t][node];
        }
        return val[t][node];
    }

    public double predict(double[] x) {
        double p = base;
        for (int t = 0; t < nTrees; t++) {
            p += predictTree(t, x);
        }
        return p;
    }
}
