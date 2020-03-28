// 
// simulation.h : contains the function declarations for simulation.cpp

/*
# 2.0
Bugfix: if first writout simulation time was too small, an error occured (i guess this was the problem)

# 1.9
Bugfix: if the first step was unsuccessful and logarithmic config writeout was requested, the program stepped into a state with 0 stepSize

# 1.8
* stepSizeBeforeWriteout is implemented: after writeout, the new stepSize is loaded from the suggested stepSizeBeforeWriteout
* subConfigTimesType == 'b' is implemented: exponential writeout frequency
* weight expects a char variable for the weight function

# 1.7
* improved debug tools:
    * they are moved to here from simulation_data (so that getElement can be used)
    * printOut is a template
    * Ax_nonSparse_to_file prints out Ax in a nonsparse way
* bool interrupt() collects all the interrupt conditions
* calcJacobianAndSpeedsAtTimes calculates faster bc
    * subSum calculated in the same loop (1 loop is eliminated)
    * exp2 is used instead of exp
* calcJacobianAndSpeedsFromPrev calculates the stresses properly even when the force is not monotonous
* logfile: contains stepSize too (no need to calculate), and cols are numbered
* stepSize is a bit closer to the requested value if sub-config-times is defined
* if heavisideCutoff, after cut-off, no values will be stored
* weight is introduced to be able modify the implicit part's weight more easily

# 1.6
stage IV checks if pH->getMaxErrorRatioSqr() < 1, as it was before, but now getMaxErrorRatioSqr() returns maxErrorToleranceRatioSq / dipolePrecisitySq

# 1.5
* remainder is removed and normalize does a loop
* strain is calculated again in the right time: before normalization
* memory for Ax and Ai (re)allocated only when necessary

# 1.4
* std::remainder is added when modifying x with Δx
* solveEQSys are relabelled

# 1.3
std::isfinite is useless with ffast-math; ranges are used instead of std::isfinite, it is implemented in debugging functions isFinite and solveEQSys too

# 1.2
solveEQSys contains label to make error messages easier to understand; calcGSolveAndUpdate also takes one therefore

# 1.1
 * bugfix: sD->subConfigDelay >= sD->subconfigDistanceCounter was checked but <= was expected
 * if subConfigTimes modifies sD->stepSize, sD->subConfigDelay is increased to trigger file writeout next time
 * checks for umfpack_di_solve's return value: if it says that the matrix is singular, nonfinite x values are zeroed out

# 1.0
subConfigTimes is implemented to print subconfigs at give times

# 0.9
* calcGSolveAndUpdate is implemented, it simplifies the code
* unused functions are removed
* code reformatted and rearranged; pragma regions added
* energies are calculated a bit differently bc speeds at the end of the 2nd small step is not available at the end of the 2nd large step.

# 0.8
Bugfix: calcJacobianFromPrev didn't calculate the new speed values at the end of the first small step, possibly lead to higher number of unsuccessful step if stress wasn't constant. calcJacobianAndSpeedsFromPrev does the right job

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

#define VERSION_simulation 1.9

#include "dislocation.h"
#include "precision_handler.h"
#include "simulation_data.h"
#include "stress_protocol.h"
#include "sim_utils.h"

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

#pragma region Directly called from run

        // calculates the forces (therefore, the speed too) between all d-d and d-p (d: dislocation, p: fixed point defect) if dislocations are at dis and the force will be calculated from simTime
        void calculateSpeedsAtTime(const std::vector<DislwoB>& dis, std::vector<double>& forces, double simTime) const;

        /**
        @brief calcJacobianAndSpeedsAtTimes:    calculates the Jacobian matrix containing the field derivatives multiplied with stepsize; modifies Ai, Ax, Ap, indexes, dVec; also calculates the forces at two different time
        @param stepsize:                        how large time step should be made
        @param dislocs:                         the actual positions of the dislocations
        @param forces_A:                        the estimated speeds of the particles at simTime_A
        @param forces_B:                        the other estimated speeds of the particles at simTime_B
        @param simTime_A:                       the which time should the external stress be evaluated to estimate the forces
        @param simTime_B:                       the which other time should the external stress be evaluated to estimate the forces
        @return int:                            totalElementCounter, the total number of nonezero elements in the matrix J_{i,j}^k
        */
        int calcJacobianAndSpeedsAtTimes(double stepsize, const std::vector<DislwoB>& dislocs, std::vector<double>& forces_A, std::vector<double>& forces_B, double simTime_A, double simTime_B);

        /**
        @brief calcJacobian:    like calcJacobianAndSpeedsAtTimes: calculates the Jacobian matrix containing the field derivatives multiplied with stepsize; modifies Ai, Ax, Ap, indexes, dVec; also calculates the force but at only 1 time point
        @param stepsize:        how large time step should be made
        @param dislocs:         the actual positions of the dislocations
        @param forces:          the speeds of the particles
        @param simTime:         the value of the external force will be calcualted from tim
        @return int:            totalElementCounter, the total number of nonezero elements in the matrix J_{i,j}^k
        */
        int calcJacobianAndSpeedsAtTime(double stepsize, const std::vector<DislwoB>& dislocs, std::vector<double>& forces, double simTime);

        // Calculates the new Jacobian J_{i,j}^k from the previous one by halfing the non-diagonal elements and also the weights
        // baseTime: the the time the Jacobian was in the correct state
        // newTime:  the new time point, where the new configuration should be predicted
        // oldTime:  the old time point, where the last prediction was made with the Jacobian
        void calcJacobianAndSpeedsFromPrev(double baseTime, double newTime, double oldTime);

        /**
        @brief calcGSolveAndUpdate: calculates the new g vector, solves the linear equations and updates the dislocation positions
        @param new_disloc:          the container of the target dislocation arrangement
        @param old_config const:    the actual dislocation configuration
        @param stepSize:            the time size of the step
        @param endSpeed const:      the estimated speeds at the end of the step
        @param initSpeedconst:      the speeds at the beginning of the step
        */
        void calcGSolveAndUpdate(std::vector<DislwoB>& new_disloc, const std::vector<DislwoB>& old_config, double stepSize, const std::vector<double>& endSpeed, const std::vector<double>& initSpeed, std::string label);

#pragma endregion

#pragma region Related to the sparse matrix handling and solving

        // modifies only Symbolic, constant for Ap, Ai, Ax; doesn't read dVec, indexes, speeds or positions
        void calculateSparseFormForJacobian();

        // solves A * Δx = g for Δx
        void solveEQSys(std::string label);

        // return the jth line element that is in between si and ei indices
        double getElement(int j, int si, int ei) const;

        // returns the i,j element of the dense matrix A
        double getElement(int i, int j) const;

#pragma endregion

#pragma region Required for analysis and outputs

        // calculates the difference of the large and two small steps and store it for all dislocs and saves the (largest relative to toleranceAndError[ID].first)
        void calculateXError();

        // calcualtes the deformation speed: Burgers' vector signed sum of the speeds
        double calculateOrderParameter(const std::vector<double>& speeds) const;

        // calculates the Burgers' vector weighted sum of the displacements
        double calculateStrainIncrement(const std::vector<DislwoB>& old, const std::vector<DislwoB>& newD) const;

        // returns if there is any reason to stop the simulation and prints out the reason
        bool interrupt() const;
#pragma endregion

#ifdef DEBUG_VERSION

        // prints out the selected container to fname        
        template <typename T>
        void printOut(std::string fname, const std::vector<T>& m_vector) const;

        // prints out size number of elements from the selected array's values to fname        
        template<typename T>
        void printOut(std::string fname, T* array, int size) const;

        // prints out the whole container for vectors and nz number of elements from dynamically allocated arrays to file container name + fname
        void printAll(std::string fname, unsigned int nz) const;

        // checks if all values in the container are finite
        bool isFinite(std::vector<double> m_vector, double lb, double ub) const;

        // checks if all x coordinate values in the container are finite
        bool isFinite(std::vector<DislwoB> disl, double lb, double ub) const;

        // checks if the first nz number of elements in the array are finite
        bool isFinite(double* m_array, size_t size, double lb, double ub) const;

        // checks if all the containers and arrays up to nz number of elements contain only finite values and
        // print out results to labeled filenames, label should match ^[\w,\s-]+
        bool isAllFinite(size_t nz, std::string label) const;

        void Ax_nonSparse_to_file(std::string fname) const;

#endif
    };

#ifdef DEBUG_VERSION
    template void Simulation::printOut(std::string, double*, int) const;
    template void Simulation::printOut(std::string, int*, int) const;
    template void Simulation::printOut(std::string, const std::vector<double>&) const;
    template void Simulation::printOut(std::string, const std::vector<DislwoB>&) const;
#endif
}

#endif
