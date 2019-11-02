// 
// simulation.h : contains the function declarations for simulation.cpp

/*
# 0.7
* calcJacobian returns the total number of non0 values in Ax, useful for debugging purposes
* getElement got new interface taking arguments of i, j; old getElement became private
* calcJacobianAndSpeedsAtTime calculates the Jacobian and the speeds parallel
* calcJacobianAndSpeedsAtTimes calculates the speeds at two different external stress values
* calcJacobianFromPrev calculates the Jacobian for the first small step from the one large step
* calcG and integrate are removed, because in each cases everything needs to be calculated a bit differently
* functions are commented

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

#define VERSION_simulation 0.7

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
        // constructor: saves the sD pointer, creates pH object, sets output format
        Simulation(std::shared_ptr<SimulationData> sD);

        // starts the simulation, all required data are present; modifies the state of the system and creates files
        void run();

    private:
        std::shared_ptr<SimulationData> sD;
        std::unique_ptr<PrecisionHandler> pH;

        // calculates the difference of the large and two small steps and store it for all dislocs and saves the (largest relative to toleranceAndError[ID].first)
        void calculateXError();
        
        // modifies only Symbolic, constant for Ap, Ai, Ax; doesn't read dVec, indexes, speeds or positions
        void calculateSparseFormForJacobian();

        // solves A * Δx = g for Δx
        void solveEQSys();

        double calculateOrderParameter(const std::vector<double>& speeds) const;

        double getElement(int j, int si, int ei) const;

        // returns the i,j element of the dense matrix A
        double getElement(int i, int j) const;
        
        // calculates the forces (therefore, the speed too) between all d-d and d-p (d: dislocation, p: fixed point defect) at two given external stresses
        void calculateSpeedsAtStresses(const std::vector<DislwoB>& dis, std::vector<double>& forces_A, std::vector<double>& forces_B, double extStress_A, double extStress_B) const;

        // calculates the forces (therefore, the speed too) between all d-d and d-p (d: dislocation, p: fixed point defect) if dislocations are at dis and the force will be calculated from simTime
        void calculateSpeedsAtTime(const std::vector<DislwoB>& dis, std::vector<double>& forces, double simTime) const;

        // calculates the forces (therefore, the speed too) between all d-d and d-p (d: dislocation, p: fixed point defect) if dislocations are at dis and the external forces will be calculated from simTime_A and simTime_B
        void calculateSpeedsAtTimes(const std::vector<DislwoB>& dis, std::vector<double>& forces_A, std::vector<double>& forces_B, double simTime_A, double simTime_B) const;

        // calculates the Burgers' vector weighted sum of the displacements
        double calculateStrainIncrement(const std::vector<DislwoB>& old, const std::vector<DislwoB>& newD) const;

        /**
        @brief calcJacobian:    calculates the Jacobian matrix containing the field derivatives multiplied with stepsize; modifies Ai, Ax, Ap, indexes, dVec
        @param stepsize:        how large time step should be made
        @param dislocs:         the actual positions of the dislocations
        @return int:            totalElementCounter, the total number of nonezero elements in the matrix J_{i,j}^k
        */
        int calcJacobian(double stepsize, const std::vector<DislwoB>& dislocs);

        // like calcJacobianAndSpeedsAtTime but calculates the forces at two given time point where the forces can be different due to the load protocol
        int calcJacobianAndSpeedsAtTimes(double stepsize, const std::vector<DislwoB>& dislocs, std::vector<double>& forces_A, std::vector<double>& forces_B, double simTime_A, double simTime_B);

        // like calcJacobian but also calculates the forces at the given time
        int calcJacobianAndSpeedsAtTime(double stepsize, const std::vector<DislwoB>& dislocs, std::vector<double>& forces, double simTime);

        // Calculates the new Jacobian J_{i,j}^k from the previous one by halfing the non-diagonal elements and also the weights
        void calcJacobianFromPrev();
    };
}

#endif
