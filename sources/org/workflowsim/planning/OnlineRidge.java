package org.workflowsim.planning;

/**
 * Online ridge regression with closed-form refits from sufficient statistics (X'X and X'y are accumulated, so the model over ALL
 * data seen so far is recovered exactly at every refit, as in the TCC predictor). Two targets share the same design matrix.
 */
public class OnlineRidge {
    private final int d;
    private final double lambda;
    private final double[][] xtx;
    private final double[][] xty;     // [target][d]
    private final double[][] beta;
    private int n = 0;

    public OnlineRidge(int d, int targets, double lambda) {
        this.d = d;
        this.lambda = lambda;
        xtx = new double[d][d];
        xty = new double[targets][d];
        beta = new double[targets][d];
    }

    public int size() {
        return n;
    }

    public void add(double[] x, double[] y) {
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                xtx[i][j] += x[i] * x[j];
            }
            for (int t = 0; t < y.length; t++) {
                xty[t][i] += x[i] * y[t];
            }
        }
        n++;
    }

    public void fit() {
        for (int t = 0; t < xty.length; t++) {
            double[][] a = new double[d][d + 1];
            for (int i = 0; i < d; i++) {
                for (int j = 0; j < d; j++) {
                    a[i][j] = xtx[i][j] + (i == j && i > 0 ? lambda : 0.0);   // the bias term is not penalised
                }
                a[i][d] = xty[t][i];
            }
            beta[t] = solve(a, d);
        }
    }

    public double predict(int t, double[] x) {
        double v = 0;
        for (int i = 0; i < d; i++) {
            v += beta[t][i] * x[i];
        }
        return v;
    }

    /** Gaussian elimination with partial pivoting on an augmented p x (p+1) system. */
    static double[] solve(double[][] aug, int p) {
        for (int col = 0; col < p; col++) {
            int piv = col;
            for (int r = col + 1; r < p; r++) {
                if (Math.abs(aug[r][col]) > Math.abs(aug[piv][col])) {
                    piv = r;
                }
            }
            double[] tmp = aug[col];
            aug[col] = aug[piv];
            aug[piv] = tmp;
            if (Math.abs(aug[col][col]) < 1e-12) {
                continue;
            }
            for (int r = col + 1; r < p; r++) {
                double f = aug[r][col] / aug[col][col];
                for (int c = col; c <= p; c++) {
                    aug[r][c] -= f * aug[col][c];
                }
            }
        }
        double[] b = new double[p];
        for (int i = p - 1; i >= 0; i--) {
            if (Math.abs(aug[i][i]) < 1e-12) {
                b[i] = 0.0;
                continue;
            }
            double sum = aug[i][p];
            for (int j = i + 1; j < p; j++) {
                sum -= aug[i][j] * b[j];
            }
            b[i] = sum / aug[i][i];
        }
        return b;
    }
}
