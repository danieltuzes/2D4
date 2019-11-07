// 
// simulation_data.h : contains the function declaration for simulation_data.cpp, project_parser.h, simulation.h

/*
# 0.6
subConfigTimes is added

# 0.5
Eliminated ic and macros according to constants.h v0.5

# 0.4
* disl_sorted were twice as large as needed, values were copied there two times

# 0.3
* Dislocations are stored without their Burgers' vector, old dislocation code is removed from the source
* Unused phython binding is removed

# 0.2
Memory allocation is checked in Release mode too

# 0.1
The first version tracked file
*/

#ifndef SDDDST_CORE_SIMULATION_DATA_H
#define SDDDST_CORE_SIMULATION_DATA_H

#define VERSION_simulation_data 0.6

#include "dislocation.h"
#include "point_defect.h"
#include "Fields/Field.h"
#include "stress_protocol.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#include <vector>
#include <map>
#include <memory>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace sdddstCore {

    class SimulationData
    {
    public:
        /**
         * @brief SimulationData can be created with giving the file path for the dislocation and the point defect
         * data files
         * @param dislocationDataFilePath
         * @param pointDefectDataFilePath
         */
        SimulationData(const std::string& dislocationConfigurationPath, const std::string& pointDefectDataFilePath);

#pragma region utility functions

        /// Data file handling utilities
        void readDislocationDataFromFile(std::string dislocationDataFilePath);
        void writeDislocationDataToFile(std::string dislocationDataFilePath) const;

        void readPointDefectDataFromFile(std::string pointDefectDataFilePath);
        void writePointDefectDataToFile(std::string pointDefectDataFilePath) const;

        /// Other utilities
        void initSimulationVariables();
        void updateCutOff();

        // returns the Burgers' vector type based on the index value ID; first half: +1; second half: -1
        int b(unsigned int ID) const;

        // returns true if Burgers' vector for the IDth dislocation is positive; false otherwise
        bool is_pos_b(unsigned int ID) const;

#pragma endregion

#pragma region debugging function tools

        // prints out totalElementCounter number of elements from Ax and all elements from dVec to file fname + ".txt"; helps debugging
        void printAxD(std::string fname, unsigned int totalElementCounter) const;

        // check if the container is finite
        bool isFinite(std::vector<double> m_vector);

        // checks if the dislocation has finite x coordinate
        bool isFinite(std::vector<DislwoB> disl);

        // check if the dynamically allocated array is finite
        bool isFinite(double* m_array, size_t size);

        // checks for all container if the are finite
        bool isAllFinite(size_t nz, std::string label);

#pragma endregion

#pragma region data fields

        // Valid dislocation position data -> state of the simulation at simTime; the sorted dislocations, Burger's vector is not needed
        std::vector<DislwoB> disl_sorted;

        // the order of the dislocations, useful for writing out in the original order
        std::vector<unsigned int> disl_order;

        //The positions of the fix points
        std::vector<PointDefect> points;

        // the g vector from the calculations
        std::vector<double> g;

        // Used to store initial speeds for the big step and for the first small step
        std::vector<double> initSpeed;

        // Used to store initial speeds for the second small step
        std::vector<double> initSpeed2;

        // Stores speed during the NR iteration for the big step and for the first small step
        std::vector<double> speed;

        // Stores speed during the NR iteration for the second small step
        std::vector<double> speed2;

        // Stores the d values for the integration scheme
        std::vector<double> dVec;

        // The dislocation data after the big step
        std::vector<DislwoB> bigStep_sorted;

        // The dislocation data after the first small step
        std::vector<DislwoB> firstSmall_sorted;

        // The dislocation data after the second small step
        std::vector<DislwoB> secondSmall_sorted;

        // UMFPack specified sparse format stored Jacobian
        int* Ap;        // index values i ∈ [Ap[j], Ap[j+1]) is used to determine the row values by Ai[i] for which A_{i,j} is non zero
        int* Ai;        // from Ap[j] to Ap[j+1] it stores the row index for the non zero elements in A_{i,j}
        double* Ax;     // The values of the sparse matrix, row-column order
        double* x;      // for which the linear equations will be solved; Ax * Δx = g for Δx

        // UMFPack required variables
        double* null;
        void* Symbolic, * Numeric;

        // Diagonal indexes in the Jacobian
        std::vector<int> indexes;

#pragma endregion

#pragma region variables

        // The used interaction field
        Field tau;

        // The value of the cut off multiplier read from input "cutoff-multiplier"
        double cutOffMultiplier;

        // The value of the cutoff = cutOffMultiplier / sqrt(dc)
        double cutOff;

        // The value of the cutoff^2
        double cutOffSqr;

        // The value of the 1/(cutoff^2)
        double onePerCutOffSqr;

        // Precisity of the simulation, set from position-precision
        double prec;

        // Count of the point defects in the system
        unsigned int pc;

        // Count of the dislocations in the system
        unsigned int dc;

        // Simulation time limit. After it is reached there should be no more calculations
        double timeLimit;

        // The current step size of the simulation
        double stepSize;

        // The current time in the simulation
        double simTime;

        // Scaling factor for point defect interaction strength calculation
        double KASQR;

        // Interaction strength between a point defect and a dislocation
        double A;

        // Number of the successfully finished steps
        size_t succesfulSteps;

        // Number of the unsuccessfully finished steps
        size_t failedSteps;

        // If strain is calculated, the total accumulated strain during the simulation
        double totalAccumulatedStrainIncrease;

        // True, if a simulation limit is set on the total accumulated strain
        bool isStrainIncreaseLimit;

        // The total accumulated strain increase limit if set
        double totalAccumulatedStrainIncreaseLimit;

        // True, if there is an upper limit set for the step size
        bool isMaxStepSizeLimit;

        // The upper limit of a step size if it is set
        double maxStepSizeLimit;

        // True if a simulation time limit is set for the simulation
        bool isTimeLimit;

        // True if a step count limit is set for the simulation
        bool isStepCountLimit;

        // The step count limit if set
        unsigned int stepCountLimit;

        // True if the strain should be calculated during the simulation
        bool calculateStrainDuringSimulation;

        // True if the order parameter should be calculated during the simulation
        bool orderParameterCalculationIsOn;

        // The standard log entries will be written into this stream
        std::ofstream standardOutputLog;

        // The path to the initial dislocation configuration
        std::string startDislocationConfigurationPath;

        // The path to the final dislocation configuration
        std::string endDislocationConfigurationPath;

        // External stress can be applied to the simulation with a specified protocol
        std::unique_ptr<StressProtocol> externalStressProtocol;

        // True if avalanches should be counted for limit
        bool countAvalanches;

        // This speed threshold is used if avalanche counting is needed
        double avalancheSpeedThreshold;

        // How many avalanches should be recorded before stop
        unsigned int avalancheTriggerLimit;

        // The current count of avalanches
        unsigned int avalancheCount;

        // True if avalanche limit is used and speed is above threshold
        bool inAvalanche;

        // True if subconfigs should be saved
        bool isSaveSubConfigs;

        // The path (and fname prefix) where the sub configs should be saved
        std::string subConfigPath;

        // The number of successful steps between two sub config output
        unsigned int subConfigDelay;

        // The number of successful steps between two sub config output during avalanches if avalanche detection is on
        unsigned int subConfigDelayDuringAvalanche;

        // The number of elapsed steps since the last subconfig written
        unsigned int subconfigDistanceCounter;

        // subconfigs must be written out at every simTime / subConfigTimes elment of N
        double subConfigTimes;

        // What kind of stress state should be used
        StressProtocolStepType currentStressStateType;

        // Cutoff multiplier changing threshold
        double speedThresholdForCutoffChange;

        bool isSpeedThresholdForCutoffChange;
#pragma endregion

    private:
    };

}

#endif
