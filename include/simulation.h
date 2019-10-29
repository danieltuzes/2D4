// 
// simulation.h : contains the function declarations for simulation.cpp

/*
# 0.6
* got rid off StressProtocol's calcExtStress and getExtStress; introduced extStress
* calculateSpeeds is removed, calculateSpeedsAtStress is introduced instead
* private member of Simulations are eliminated as they were only used in the class once or twice
* rearranged code: more interesting parts come first
* some new name convention inside the code

# 0.5
* calculateSpeedsAtStresses calculates the speeds at two different external stress values, eliminates 1 speed calculation
* unused function declarations with Dislocation classes are removed
* iteration loop unrolled

# 0.4
* Burgers' vector value is checked and used on a bas of bool, no multiplication
* stepStages are moved to run()
* first run code part in stepStageI is moved out to run()
* prints out start time

# 0.3
USE_POINT_DEFECTS is checked to eliminate codes on point defects

# 0.2
* Dislocations are stored without Burgers' vector, it is deduced from dislocation ID
* unused python bindings are removed

# 0.1
First version tracked source
*/

#ifndef SDDDST_CORE_SIMULATION_H
#define SDDDST_CORE_SIMULATION_H

#define VERSION_simulation 0.6

#include "dislocation.h"
#include "precision_handler.h"
#include "simulation_data.h"
#include "stress_protocol.h"

#include <memory>
#include <sstream>
#include <vector>

namespace sdddstCore {

    class Simulation
    {
    public:
        Simulation(std::shared_ptr<SimulationData> sD);

        void calculateXError();
        void calculateSparseFormForJacobian();
        void solveEQSys();

        double calculateOrderParameter(const std::vector<double>& speeds) const;

        double getElement(int j, int si, int ei) const;

        double getSimTime() const;

        void run();

        // with DislwoB: dislocation without Burger's vector
        void integrate(double stepsize, std::vector<DislwoB>& newDislocation, const std::vector<DislwoB>& old, bool useSpeed2, bool calculateInitSpeed, StressProtocolStepType origin, StressProtocolStepType end);

        // calculates the forces (therefore, the speed too) between all d-d and d-p (d: dislocation, p: fixed point defect) at a given external stress
        void calculateSpeedsAtStress(const std::vector<DislwoB>& dis, std::vector<double>& forces, double extStress) const;

        // calculates the forces (therefore, the speed too) between all d-d and d-p (d: dislocation, p: fixed point defect) at two given external stresses
        void calculateSpeedsAtStresses(const std::vector<DislwoB>& dis, std::vector<double>& forces_A, std::vector<double>& forces_B, double extStress_A, double extStress_B) const;

        // calculates the forces (therefore, the speed too) between all d-d and d-p (d: dislocation, p: fixed point defect) if dislocations are at dis and the force will be calculated from simTime
        void calculateSpeedsAtTime(const std::vector<DislwoB>& dis, std::vector<double>& forces, double simTime) const;
        
        // calculates the forces (therefore, the speed too) between all d-d and d-p (d: dislocation, p: fixed point defect) if dislocations are at dis and the external forces will be calculated from simTime_A and simTime_B
        void calculateSpeedsAtTimes(const std::vector<DislwoB>& dis, std::vector<double>& forces_A, std::vector<double>& forces_B, double simTime_A, double simTime_B) const; 

        // calculates the g vector; modifies only g
        void calculateG(double stepsize, const std::vector<DislwoB>& newDislocation, const std::vector<DislwoB>& old, bool useSpeed2, bool calculateInitSpeed, StressProtocolStepType origin, StressProtocolStepType end) const;

        double calculateStrainIncrement(const std::vector<DislwoB>& old, const std::vector<DislwoB>& newD) const;

        void calcJacobian(double stepsize, const std::vector<DislwoB>& dislocs);

    private:
        std::shared_ptr<SimulationData> sD;
        std::unique_ptr<PrecisionHandler> pH;
    };
}

#endif
