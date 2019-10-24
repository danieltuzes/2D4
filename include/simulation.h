
// simulation.h : contains the function declarations for simulation.cpp

/*
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

#define VERSION_simulation 0.3

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

        void stepStageI();
        void stepStageII();
        void stepStageIII();

        // with Dislocation: dislocation with Burgers' vector
        void integrate(double stepsize, std::vector<Dislocation>& newDislocation, const std::vector<Dislocation>& old, bool useSpeed2, bool calculateInitSpeed, StressProtocolStepType origin, StressProtocolStepType end);

        // calculates the forces (therefore, the speed too) between all d-d and d-p (d: dislocation, p: fixed point defect)
        void calculateSpeeds(const std::vector<Dislocation>& dis, std::vector<double>& res) const;
        void calculateG(double stepsize, const std::vector<Dislocation>& newDislocation, const std::vector<Dislocation>& old, bool useSpeed2, bool calculateInitSpeed, bool useInitSpeedForFirstStep, StressProtocolStepType origin, StressProtocolStepType end) const;
        void calculateJacobian(double stepsize, const std::vector<Dislocation>& data);

        double calculateStrainIncrement(const std::vector<Dislocation>& old, const std::vector<Dislocation>& newD) const;

        // with DislwoB: dislocation without Burger's vector
        void integrate(double stepsize, std::vector<DislwoB>& newDislocation, const std::vector<DislwoB>& old, bool useSpeed2, bool calculateInitSpeed, StressProtocolStepType origin, StressProtocolStepType end);

        // calculates the forces (therefore, the speed too) between all d-d and d-p (d: dislocation, p: fixed point defect)
        void calculateSpeeds(const std::vector<DislwoB>& dis, std::vector<double>& res) const;
        void calculateG(double stepsize, const std::vector<DislwoB>& newDislocation, const std::vector<DislwoB>& old, bool useSpeed2, bool calculateInitSpeed, bool useInitSpeedForFirstStep, StressProtocolStepType origin, StressProtocolStepType end) const;

        double calculateStrainIncrement(const std::vector<DislwoB>& old, const std::vector<DislwoB>& newD) const;

        void calculateJacobian(double stepsize, const std::vector<DislwoB>& data);

    private:
        double lastWriteTimeFinished;
        double startTime;
        bool initSpeedCalculationIsNeeded;
        double energy;

        std::shared_ptr<SimulationData> sD;
        std::unique_ptr<PrecisionHandler> pH;
    };
}

#endif
