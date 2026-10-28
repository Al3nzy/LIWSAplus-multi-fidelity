package org.workflowsim.examples.experiments;

import java.io.File;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.LinkedList;
import java.util.List;
import org.cloudbus.cloudsim.CloudletSchedulerSpaceShared;
import org.cloudbus.cloudsim.DatacenterCharacteristics;
import org.cloudbus.cloudsim.HarddriveStorage;
import org.cloudbus.cloudsim.Host;
import org.cloudbus.cloudsim.Pe;
import org.cloudbus.cloudsim.Storage;
import org.cloudbus.cloudsim.VmAllocationPolicySimple;
import org.cloudbus.cloudsim.VmSchedulerTimeShared;
import org.cloudbus.cloudsim.core.CloudSim;
import org.cloudbus.cloudsim.provisioners.BwProvisionerSimple;
import org.cloudbus.cloudsim.provisioners.PeProvisionerSimple;
import org.cloudbus.cloudsim.provisioners.RamProvisionerSimple;
import org.workflowsim.CondorVM;
import org.workflowsim.Job;
import org.workflowsim.WorkflowDatacenter;
import org.workflowsim.WorkflowEngine;
import org.workflowsim.WorkflowPlanner;
import org.workflowsim.examples.planning.RunMetricsCalculator;
import org.workflowsim.utils.ClusteringParameters;
import org.workflowsim.utils.OverheadParameters;
import org.workflowsim.utils.Parameters;
import org.workflowsim.utils.ReplicaCatalog;

/**
 * Runs HEFT and Min-Min once per workflow inside the real WorkflowSim simulation, and writes (a) their
 * cloudlet-id -> VM-id assignments (the warm-start seeds shared by every algorithm in every experiment) and
 * (b) their simulator-measured makespan and cost. Seeds are cached so that all conditions use identical seeds.
 * Usage: SeedBuilder <outDir> <wf1,wf2,...>
 */
public class SeedBuilder {
    public static void main(String[] args) throws Exception {
        String out = args[0];
        new File(out).mkdirs();
        try (PrintWriter summary = new PrintWriter(out + "/heuristics_sim.csv")) {
            summary.println("workflow,algorithm,sim_makespan,sim_cost");
            for (String wf : args[1].split(",")) {
                for (String alg : new String[]{"HEFT", "MINMIN"}) {
                    run(out, wf, alg, summary);
                }
            }
        }
    }

    /** Replays a genotype inside the full WorkflowSim simulation; returns {makespan, cost} as measured by the simulator. */
    static double[] simulateForced(String wf, int[] genotype) throws Exception {
        org.workflowsim.planning.LIWSAPlanningAlgorithm.CONFIG_FORCED_GENOTYPE = genotype;
        try {
            OverheadParameters op = new OverheadParameters(0, null, null, null, null, 0);
            ClusteringParameters cp = new ClusteringParameters(0, 0, ClusteringParameters.ClusteringMethod.NONE, null);
            Parameters.init(16, "config/dax/" + wf + ".xml", null, null, op, cp,
                    Parameters.SchedulingAlgorithm.STATIC, Parameters.PlanningAlgorithm.LIWSA, null, 0);
            Parameters.setCostModel(Parameters.CostModel.VM);
            ReplicaCatalog.init(ReplicaCatalog.FileSystem.LOCAL);
            CloudSim.init(1, Calendar.getInstance(), false);
            WorkflowDatacenter dc = createDatacenter("DC_forced");
            WorkflowPlanner planner = new WorkflowPlanner("planner_forced", 1);
            WorkflowEngine engine = planner.getWorkflowEngine();
            List<CondorVM> vms = createVMs(engine.getSchedulerId(0));
            engine.submitVmList(vms, 0);
            engine.bindSchedulerDatacenter(dc.getId(), 0);
            CloudSim.startSimulation();
            List<Job> jobs = engine.getJobsReceivedList();
            CloudSim.stopSimulation();
            RunMetricsCalculator.Result m = RunMetricsCalculator.compute(jobs, vms.size(), 2000.0);
            return new double[]{m.makespan, m.cost};
        } finally {
            org.workflowsim.planning.LIWSAPlanningAlgorithm.CONFIG_FORCED_GENOTYPE = null;
        }
    }

    private static void run(String out, String wf, String alg, PrintWriter summary) throws Exception {
        String dax = "config/dax/" + wf + ".xml";
        Parameters.PlanningAlgorithm pa = alg.equals("HEFT") ? Parameters.PlanningAlgorithm.HEFT : Parameters.PlanningAlgorithm.INVALID;
        Parameters.SchedulingAlgorithm sa = alg.equals("HEFT") ? Parameters.SchedulingAlgorithm.STATIC : Parameters.SchedulingAlgorithm.MINMIN;
        OverheadParameters op = new OverheadParameters(0, null, null, null, null, 0);
        ClusteringParameters cp = new ClusteringParameters(0, 0, ClusteringParameters.ClusteringMethod.NONE, null);
        Parameters.init(16, dax, null, null, op, cp, sa, pa, null, 0);
        Parameters.setCostModel(Parameters.CostModel.VM);
        ReplicaCatalog.init(ReplicaCatalog.FileSystem.LOCAL);
        CloudSim.init(1, Calendar.getInstance(), false);
        WorkflowDatacenter dc = createDatacenter("DC_" + alg + "_" + wf);
        WorkflowPlanner planner = new WorkflowPlanner("planner_" + alg + "_" + wf, 1);
        WorkflowEngine engine = planner.getWorkflowEngine();
        List<CondorVM> vms = createVMs(engine.getSchedulerId(0));
        engine.submitVmList(vms, 0);
        engine.bindSchedulerDatacenter(dc.getId(), 0);
        CloudSim.startSimulation();
        List<Job> jobs = engine.getJobsReceivedList();
        CloudSim.stopSimulation();
        RunMetricsCalculator.Result m = RunMetricsCalculator.compute(jobs, vms.size(), 2000.0);
        try (PrintWriter pw = new PrintWriter(out + "/" + wf + "_" + alg + ".txt")) {
            for (Job job : jobs) {
                if (job.getClassType() == Parameters.ClassType.STAGE_IN.value) {
                    continue;
                }
                for (org.workflowsim.Task t : job.getTaskList()) {
                    pw.println(t.getCloudletId() + " " + job.getVmId());
                }
            }
        }
        summary.println(wf + "," + alg + "," + m.makespan + "," + m.cost);
        summary.flush();
        System.out.println(wf + " " + alg + " sim makespan=" + m.makespan + " cost=" + m.cost);
    }

    static List<CondorVM> createVMs(int userId) {
        LinkedList<CondorVM> list = new LinkedList<>();
        int vmId = 0;
        for (double[] type : StandaloneHarness.VM_TYPES) {
            for (int k = 0; k < (int) type[5]; k++) {
                list.add(new CondorVM(vmId++, userId, type[0], 1, (int) type[3], (long) type[1], (long) type[4],
                        "Xen", type[2], 0.0, 0.0, 0.0, new CloudletSchedulerSpaceShared()));
            }
        }
        return list;
    }

    static WorkflowDatacenter createDatacenter(String name) {
        int numHosts = 4;
        List<Host> hostList = new ArrayList<>();
        for (int h = 0; h < numHosts; h++) {
            List<Pe> peList = new ArrayList<>();
            for (int p = 0; p < 16; p++) {
                peList.add(new Pe(p, new PeProvisionerSimple(2000)));
            }
            hostList.add(new Host(h, new RamProvisionerSimple(8192), new BwProvisionerSimple(100000),
                    1000000L, peList, new VmSchedulerTimeShared(peList)));
        }
        DatacenterCharacteristics dc = new DatacenterCharacteristics("x86", "Linux", "Xen", hostList, 10.0, 3.0, 0.05, 0.1, 0.1);
        WorkflowDatacenter datacenter = null;
        try {
            HarddriveStorage s = new HarddriveStorage(name, 1e12);
            s.setMaxTransferRate(100);
            LinkedList<Storage> sl = new LinkedList<>();
            sl.add(s);
            datacenter = new WorkflowDatacenter(name, dc, new VmAllocationPolicySimple(hostList), sl, 0);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return datacenter;
    }
}
