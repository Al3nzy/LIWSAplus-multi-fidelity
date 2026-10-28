package org.workflowsim.examples.experiments;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import org.cloudbus.cloudsim.CloudletSchedulerSpaceShared;
import org.workflowsim.CondorVM;
import org.workflowsim.Task;
import org.workflowsim.WorkflowParser;
import org.workflowsim.utils.ClusteringParameters;
import org.workflowsim.utils.OverheadParameters;
import org.workflowsim.utils.Parameters;
import org.workflowsim.utils.ReplicaCatalog;

/**
 * Runs planning algorithms OUTSIDE the discrete-event simulator, on exactly the task graph and the
 * 16-VM pool used by the benchmark (4 x {250,500,1000,2000} MIPS at {0.15,0.30,0.60,0.90} per second,
 * 160 Mbit/s). The planners only ever consume (task list, VM list), so this is result-identical to the
 * planner call inside WorkflowPlanner.processPlanning() and avoids the simulator's per-run overhead.
 */
public class StandaloneHarness {

    public static final double[][] VM_TYPES = {
        { 250.0, 160.0, 0.15, 512, 10000, 4},
        { 500.0, 160.0, 0.30, 512, 10000, 4},
        {1000.0, 160.0, 0.60, 512, 10000, 4},
        {2000.0, 160.0, 0.90, 512, 10000, 4},
    };

    /** Price schemes (per second) for the four speed tiers: original, strictly linear, convex speed premium. */
    static {
        String sc = System.getProperty("pricing", "original");
        double[] rates = sc.equals("linear") ? new double[]{0.15, 0.30, 0.60, 1.20}
                : sc.equals("premium") ? new double[]{0.15, 0.35, 0.85, 2.00}
                : new double[]{0.15, 0.30, 0.60, 0.90};
        for (int i = 0; i < 4; i++) {
            VM_TYPES[i][2] = rates[i];
        }
    }

    public static List<Task> loadTasks(String daxPath) {
        OverheadParameters op = new OverheadParameters(0, null, null, null, null, 0);
        ClusteringParameters cp = new ClusteringParameters(0, 0, ClusteringParameters.ClusteringMethod.NONE, null);
        Parameters.init(16, daxPath, null, null, op, cp,
                Parameters.SchedulingAlgorithm.STATIC, Parameters.PlanningAlgorithm.INVALID, null, 0);
        ReplicaCatalog.init(ReplicaCatalog.FileSystem.LOCAL);
        WorkflowParser parser = new WorkflowParser(0);
        parser.parse();
        return new ArrayList<>(parser.getTaskList());
    }

    public static List<CondorVM> makeVms() {
        LinkedList<CondorVM> list = new LinkedList<>();
        int vmId = 0;
        for (double[] type : VM_TYPES) {
            for (int k = 0; k < (int) type[5]; k++) {
                list.add(new CondorVM(vmId++, 0, type[0], 1, (int) type[3], (long) type[1], (long) type[4],
                        "Xen", type[2], 0.0, 0.0, 0.0, new CloudletSchedulerSpaceShared()));
            }
        }
        return list;
    }
}
