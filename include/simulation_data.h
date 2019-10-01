/*
 * SDDDST Simple Discrete Dislocation Dynamics Toolkit
 * Copyright (C) 2015-2019  Gábor Péterffy <peterffy95@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA
 */

#ifndef SDDDST_CORE_SIMULATION_DATA_H
#define SDDDST_CORE_SIMULATION_DATA_H

#include "dislocation.h"
#include "point_defect.h"
#include "Fields/Field.h"
#include "StressProtocols/stress_protocol.h"

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

#ifdef BUILD_PYTHON_BINDINGS
#include <boost/python.hpp>
#endif

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

        ///////////////////
        /// UTILITIES
        ///

        /// Data file handling utilities
        void readDislocationDataFromFile(std::string dislocationDataFilePath);
        void writeDislocationDataToFile(std::string dislocationDataFilePath) const;
        void readPointDefectDataFromFile(std::string pointDefectDataFilePath);
        void writePointDefectDataToFile(std::string pointDefectDataFilePath) const;

        /// Other utilities
        void initSimulationVariables();
        void updateCutOff();

        //////////////////
        /// DATA FIELDS
        ///

        //Valid dislocation position data -> state of the simulation at simTime
        std::vector<Dislocation> dislocations;

        // the sorted dislocations, Burger's vector is not needed
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

        // The value of the cut off multiplier
        double cutOffMultiplier;

        // The value of the cutoff
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

        // Count of the iterations during the NR
        unsigned int ic;

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

        // The dislocation data after the big step
        std::vector<Dislocation> bigStep;
        std::vector<DislwoB> oBigStep;

        // The dislocation data after the first small step
        std::vector<Dislocation> firstSmall;
        std::vector<DislwoB> oFirstSmall;

        // The dislocation data after the second small step
        std::vector<Dislocation> secondSmall;
        std::vector<DislwoB> oSecondSmall;

        // The used interaction field
        std::unique_ptr<Field> tau;

        // UMFPack specified sparse format stored Jacobian
        int* Ap;
        int* Ai;
        double* Ax;
        // Result data
        double* x;

        // UMFPack required variables
        double* null;
        void* Symbolic, * Numeric;

        // Diagonal indexes in the Jacobian
        std::vector<int> indexes;

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

        // The upper limit of a step size if is set
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

        // The path where the sub configs should be saved
        std::string subConfigPath;

        // The number of successful steps between two sub config output
        unsigned int subConfigDelay;

        // The number of successful steps between two sub config output during avalanches if avalanche detection is on
        unsigned int subConfigDelayDuringAvalanche;

        // The number of elapsed steps since the last subconfig written
        unsigned int subconfigDistanceCounter;

        // What kind of stress state should be used
        StressProtocolStepType currentStressStateType;

        // Cutoff multiplier changing threshold
        double speedThresholdForCutoffChange;

        bool isSpeedThresholdForCutoffChange;

#ifdef BUILD_PYTHON_BINDINGS

        Field const& getField();
        void setField(boost::python::object field);

        StressProtocol const& getStressProtocol();
        void setStressProtocol(boost::python::object protocol);

        // deleteDislocationCountRelatedData free memory what was allocated for dislocation related data
        // void deleteDislocationCountRelatedData(); */

#endif

    private:
    };

}

#endif
