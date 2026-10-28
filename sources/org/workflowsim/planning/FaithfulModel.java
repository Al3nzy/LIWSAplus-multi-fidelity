package org.workflowsim.planning;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.workflowsim.CondorVM;
import org.workflowsim.FileItem;
import org.workflowsim.Task;
import org.workflowsim.utils.Parameters.FileType;

/**
 * Event-driven re-implementation of the rules WorkflowSim applies to a statically scheduled workflow, so that the planning
 * level objective equals what the simulator measures. The analytic decoder differs from the simulator in four ways:
 *  1. a job's stage-in time is added to its length, so the VM is OCCUPIED (and billed) while inputs are fetched, instead of the
 *     transfer overlapping with other work;
 *  2. input files are fetched per VM and cached there: a file already on the VM (produced there or fetched earlier) costs nothing,
 *     and workflow input files are charged the first time a VM needs them;
 *  3. a VM executes its jobs strictly first-come-first-served in release order (the static scheduler only submits to an idle VM),
 *     not by insertion into the earliest idle gap in topological order;
 *  4. a stage-in job (110 MI on VM 0) runs first and releases the root tasks.
 * A job started on VM v at time t runs for ((long)(len + mips_v * stage)) / mips_v seconds (CloudSim truncates the extra length to
 * a long); cost is that occupied time times the VM price.
 */
public class FaithfulModel {
    final int n, m;
    final double[] len;
    final double[] mips, rate, bw;
    final int[][] parents, children, inFiles, outFiles;
    final double[] fsize;
    final int[] rootOrder;
    final int[] producer;                 // task index that produces each file (-1 for workflow inputs)
    final double t0;

    // work buffers
    private final int[] indeg;
    private final int[] fmask;
    private final int[][] queue;
    private final int[] qh, qt;
    private final boolean[] busy;
    private final double[] heapT;
    private final int[] heapK, heapS;
    private int heapN, seqCounter;
    private final int[] taskVm;
    // schedule of the last decode (for guided variation operators)
    final double[] st, fin;
    final int[] relParent, vmPrev, lastOnVm;

    public FaithfulModel(List<Task> taskOrder, List<CondorVM> vms) {
        n = taskOrder.size();
        m = vms.size();
        Map<Task, Integer> pos = new HashMap<>();
        for (int k = 0; k < n; k++) {
            pos.put(taskOrder.get(k), k);
        }
        len = new double[n];
        parents = new int[n][];
        children = new int[n][];
        inFiles = new int[n][];
        outFiles = new int[n][];
        Map<String, Integer> fid = new HashMap<>();
        List<Double> sizes = new ArrayList<>();
        for (int k = 0; k < n; k++) {
            Task t = taskOrder.get(k);
            len[k] = t.getCloudletLength();
            parents[k] = new int[t.getParentList().size()];
            for (int i = 0; i < parents[k].length; i++) {
                parents[k][i] = pos.get((Task) t.getParentList().get(i));
            }
            children[k] = new int[t.getChildList().size()];
            for (int i = 0; i < children[k].length; i++) {
                children[k][i] = pos.get((Task) t.getChildList().get(i));
            }
            List<Integer> ins = new ArrayList<>();
            List<Integer> outs = new ArrayList<>();
            for (FileItem f : t.getFileList()) {
                Integer id = fid.get(f.getName());
                if (id == null) {
                    id = sizes.size();
                    fid.put(f.getName(), id);
                    sizes.add(f.getSize());
                }
                if (f.getType() == FileType.OUTPUT) {
                    outs.add(id);
                } else if (f.isRealInputFile(t.getFileList())) {
                    ins.add(id);
                }
            }
            inFiles[k] = ins.stream().mapToInt(Integer::intValue).toArray();
            outFiles[k] = outs.stream().mapToInt(Integer::intValue).toArray();
        }
        producer = new int[sizes.size()];
        java.util.Arrays.fill(producer, -1);
        for (int k = 0; k < n; k++) {
            for (int f : outFiles[k]) {
                producer[f] = k;
            }
        }
        fsize = new double[sizes.size()];
        for (int i = 0; i < fsize.length; i++) {
            fsize[i] = sizes.get(i);
        }
        mips = new double[m];
        rate = new double[m];
        bw = new double[m];
        for (int v = 0; v < m; v++) {
            mips[v] = vms.get(v).getMips();
            rate[v] = vms.get(v).getCost();
            bw[v] = vms.get(v).getBw();
        }
        List<Integer> roots = new ArrayList<>();
        for (int k = 0; k < n; k++) {
            if (parents[k].length == 0) {
                roots.add(k);
            }
        }
        final List<Task> to = taskOrder;
        roots.sort((x, y) -> Integer.compare(to.get(x).getCloudletId(), to.get(y).getCloudletId()));
        rootOrder = roots.stream().mapToInt(Integer::intValue).toArray();
        t0 = 110.0 / mips[0] + 0.1;      // stage-in job on VM 0, plus CloudSim's 0.1 s minimum time between events
        indeg = new int[n];
        fmask = new int[fsize.length];
        queue = new int[m][n];
        qh = new int[m];
        qt = new int[m];
        busy = new boolean[m];
        heapT = new double[n + 1];
        heapK = new int[n + 1];
        heapS = new int[n + 1];
        taskVm = new int[n];
        st = new double[n];
        fin = new double[n];
        relParent = new int[n];
        vmPrev = new int[n];
        lastOnVm = new int[m];
        vmBusy = new double[m];
    }

    private void push(double t, int k) {
        int i = heapN++;
        heapT[i] = t;
        heapK[i] = k;
        heapS[i] = seqCounter++;
        while (i > 0) {
            int p = (i - 1) / 2;
            if (less(i, p)) {
                swap(i, p);
                i = p;
            } else {
                break;
            }
        }
    }

    private boolean less(int a, int b) {
        return heapT[a] < heapT[b] || (heapT[a] == heapT[b] && heapS[a] < heapS[b]);
    }

    private void swap(int a, int b) {
        double t = heapT[a]; heapT[a] = heapT[b]; heapT[b] = t;
        int k = heapK[a]; heapK[a] = heapK[b]; heapK[b] = k;
        int s = heapS[a]; heapS[a] = heapS[b]; heapS[b] = s;
    }

    private void pop() {
        heapN--;
        if (heapN > 0) {
            heapT[0] = heapT[heapN];
            heapK[0] = heapK[heapN];
            heapS[0] = heapS[heapN];
            int i = 0;
            while (true) {
                int l = 2 * i + 1, r = l + 1, s = i;
                if (l < heapN && less(l, s)) {
                    s = l;
                }
                if (r < heapN && less(r, s)) {
                    s = r;
                }
                if (s == i) {
                    break;
                }
                swap(i, s);
                i = s;
            }
        }
    }

    private double cost;
    private double[] curLen, curMips;
    private double stageSum;
    private final double[] vmBusy;

    private void dispatch(double t) {
        for (int v = 0; v < m; v++) {
            if (!busy[v] && qh[v] < qt[v]) {
                int k = queue[v][qh[v]++];
                busy[v] = true;
                st[k] = t;
                vmPrev[k] = lastOnVm[v];
                lastOnVm[v] = k;
                double stage = 0.0;
                for (int f : inFiles[k]) {
                    if ((fmask[f] & (1 << v)) == 0) {
                        stage += fsize[f] / 1000000.0 / bw[v];
                        fmask[f] |= (1 << v);
                    }
                }
                long length = (long) (curLen[k] + curMips[v] * stage);
                double dur = length / curMips[v];
                cost += dur * rate[v];
                stageSum += stage;
                vmBusy[v] += dur;
                fin[k] = t + dur;
                push(t + dur, k);
            }
        }
    }

    /**
     * Same as decode(g) but with every task length multiplied by lenMult[k] and every VM speed by mipsMult[v]: used to test how a
     * schedule planned on nominal values behaves when runtime estimates are wrong or VMs do not deliver their peak MIPS.
     */
    public double[] decodePerturbed(int[] g, double[] lenMult, double[] mipsMult) {
        double[] pl = new double[n];
        double[] pm = new double[m];
        for (int k = 0; k < n; k++) {
            pl[k] = len[k] * lenMult[k];
        }
        for (int v = 0; v < m; v++) {
            pm[v] = mips[v] * mipsMult[v];
        }
        curLen = pl;
        curMips = pm;
        double[] r = run(g);
        curLen = len;
        curMips = mips;
        return r;
    }

    /** @return {makespan, cost} as WorkflowSim would measure them for this task-to-VM assignment. */
    public double[] decode(int[] g) {
        curLen = len;
        curMips = mips;
        return run(g);
    }

    private double[] run(int[] g) {
        java.util.Arrays.fill(fmask, 0);
        java.util.Arrays.fill(busy, false);
        java.util.Arrays.fill(qh, 0);
        java.util.Arrays.fill(qt, 0);
        java.util.Arrays.fill(lastOnVm, -1);
        java.util.Arrays.fill(vmBusy, 0.0);
        stageSum = 0.0;
        java.util.Arrays.fill(relParent, -1);
        heapN = 0;
        seqCounter = 0;
        cost = 0.0;
        for (int k = 0; k < n; k++) {
            indeg[k] = parents[k].length;
            taskVm[k] = g[k];
        }
        for (int k : rootOrder) {
            queue[g[k]][qt[g[k]]++] = k;
        }
        double t = t0;
        dispatch(t);
        double makespan = 0.0;
        while (heapN > 0) {
            t = heapT[0];
            int k = heapK[0];
            pop();
            int v = g[k];
            busy[v] = false;
            makespan = Math.max(makespan, t);
            for (int f : outFiles[k]) {
                fmask[f] |= (1 << v);
            }
            for (int c : children[k]) {
                if (--indeg[c] == 0) {
                    relParent[c] = k;
                    queue[g[c]][qt[g[c]]++] = c;
                }
            }
            dispatch(t);
        }
        return new double[]{makespan, cost};
    }

    /** Estimated occupied time of task k on VM v given the assignment g: execution plus fetching inputs produced elsewhere. */
    public double estimate(int k, int v, int[] g) {
        double stage = 0.0;
        for (int f : inFiles[k]) {
            int pr = producer[f];
            if (pr >= 0 && g[pr] != v) {
                stage += fsize[f] / 1000000.0 / bw[v];
            }
        }
        return len[k] / mips[v] + stage;
    }

    /** Chain of tasks that binds the makespan of the LAST decoded schedule (VM queue or data release), latest task first. */
    public int[] criticalChain() {
        int cur = 0;
        for (int k = 1; k < n; k++) {
            if (fin[k] > fin[cur]) {
                cur = k;
            }
        }
        List<Integer> chain = new ArrayList<>();
        chain.add(cur);
        for (int guard = 0; guard < n; guard++) {
            double release = (relParent[cur] >= 0) ? fin[relParent[cur]] : t0;
            int pred;
            if (st[cur] > release + 1e-9 && vmPrev[cur] >= 0) {
                pred = vmPrev[cur];          // started late because the VM was busy
            } else {
                pred = relParent[cur];       // started when its last parent finished
            }
            if (pred < 0) {
                break;
            }
            chain.add(pred);
            cur = pred;
        }
        int[] out = new int[chain.size()];
        for (int i = 0; i < out.length; i++) {
            out[i] = chain.get(i);
        }
        return out;
    }

    /** Number of tasks on the makespan-binding chain of the last decoded schedule. */
    private int chainLength() {
        return chainStats(null)[0];
    }

    /** @return {chain length, tasks of the chain on its most used machine} (the second entry needs g). */
    private int[] chainStats(int[] g) {
        int cur = 0;
        for (int k = 1; k < n; k++) {
            if (fin[k] > fin[cur]) {
                cur = k;
            }
        }
        int len = 1;
        int[] onVm = (g == null) ? null : new int[m];
        if (onVm != null) {
            onVm[g[cur]]++;
        }
        for (int guard = 0; guard < n; guard++) {
            double release = (relParent[cur] >= 0) ? fin[relParent[cur]] : t0;
            int pred = (st[cur] > release + 1e-9 && vmPrev[cur] >= 0) ? vmPrev[cur] : relParent[cur];
            if (pred < 0) {
                break;
            }
            len++;
            cur = pred;
            if (onVm != null) {
                onVm[g[cur]]++;
            }
        }
        int mx = 0;
        if (onVm != null) {
            for (int c : onVm) {
                mx = Math.max(mx, c);
            }
        }
        return new int[]{len, mx};
    }

    public static final int N_SURR_FEATURES = 10;

    /**
     * Nominal decode plus cheap structural descriptors of the schedule, used as inputs of the online surrogate.
     * @return {nominal makespan, nominal cost, f0..f6} with f0 = 1 (bias), critical-chain share, share of machines used,
     *         bottleneck-machine utilisation, mean utilisation, share of dependencies crossing machines, staging share of busy time
     */
    public double[] features(int[] g) {
        double[] r = decode(g);
        double sumBusy = 0, maxBusy = 0;
        int used = 0;
        for (int v = 0; v < m; v++) {
            sumBusy += vmBusy[v];
            maxBusy = Math.max(maxBusy, vmBusy[v]);
            if (vmBusy[v] > 0) {
                used++;
            }
        }
        int edges = 0, cross = 0;
        for (int k = 0; k < n; k++) {
            for (int p : parents[k]) {
                edges++;
                if (g[p] != g[k]) {
                    cross++;
                }
            }
        }
        double ms = Math.max(r[0], 1e-9);
        int[] cs = chainStats(g);
        double herf = 0, slowBusy = 0;
        double[] sp = mips.clone();
        java.util.Arrays.sort(sp);
        double median = sp[m / 2];
        for (int v = 0; v < m; v++) {
            double sh = sumBusy > 0 ? vmBusy[v] / sumBusy : 0;
            herf += sh * sh;
            if (mips[v] < median) {
                slowBusy += vmBusy[v];
            }
        }
        double[] out = new double[12 + m];
        double[] head = {r[0], r[1], 1.0, cs[0] / (double) n, used / (double) m, maxBusy / ms, sumBusy / (ms * m),
                edges > 0 ? cross / (double) edges : 0.0, sumBusy > 0 ? stageSum / sumBusy : 0.0,
                herf, cs[0] > 0 ? cs[1] / (double) cs[0] : 0.0, sumBusy > 0 ? slowBusy / sumBusy : 0.0};
        System.arraycopy(head, 0, out, 0, 12);
        for (int v = 0; v < m; v++) {
            out[12 + v] = vmBusy[v] / ms;      // busy-time share of every machine (input of the richer surrogates)
        }
        return out;
    }
}
