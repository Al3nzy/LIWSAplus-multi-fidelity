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
package org.workflowsim.examples.planning;

import java.util.List;
import org.workflowsim.Job;

/**
 * Computes the metrics that apply identically to every planning algorithm
 * (HEFT, Min-Min, MLEAO, LIWSA, LIWSA-ML alike), directly from the
 * simulator's actual job results rather than from any algorithm's
 * internal state. Using one shared implementation across all four
 * example/benchmark drivers means these numbers can never drift out of
 * sync between files.
 */
public class RunMetricsCalculator {

    public static class Result {
        public double makespan;
        public double cost;
        /** Fraction in [0,1]: mean, across all VMs in the pool, of (busy time / makespan). */
        public double avgUtilization;
        /** Jain's fairness index over per-VM busy time: 1.0 = perfectly even load,
         *  1/numVMs = all load on a single VM. */
        public double fairnessIndex;
        /** (workflow total task length / fastest available VM's MIPS) / makespan. */
        public double speedup;
    }

    public static Result compute(List<Job> jobs, int numVMs, double fastestMips) {
        Result r = new Result();
        double[] vmBusyTime = new double[numVMs];
        double totalTaskLength = 0.0;

        for (Job job : jobs) {
            if (job.getClassType() == org.workflowsim.utils.Parameters.ClassType.STAGE_IN.value) {
                continue;
            }
            r.makespan = Math.max(r.makespan, job.getFinishTime());
            r.cost += job.getActualCPUTime() * job.getCostPerSec();
            int vmId = job.getVmId();
            if (vmId >= 0 && vmId < numVMs) {
                vmBusyTime[vmId] += job.getActualCPUTime();
            }
            totalTaskLength += job.getCloudletTotalLength();
        }

        double busySum = 0, busySqSum = 0;
        for (double b : vmBusyTime) { busySum += b; busySqSum += b * b; }
        r.avgUtilization = (r.makespan > 0) ? (busySum / numVMs) / r.makespan : 0.0;
        r.fairnessIndex = (busySqSum > 0) ? (busySum * busySum) / (numVMs * busySqSum) : 0.0;

        double sequentialTime = totalTaskLength / fastestMips;
        r.speedup = (r.makespan > 0) ? sequentialTime / r.makespan : 0.0;

        return r;
    }
}
