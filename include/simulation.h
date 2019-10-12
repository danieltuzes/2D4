
// simulation.h : contains the function declarations for simulation.cpp

#ifndef SDDDST_CORE_SIMULATION_H
#define SDDDST_CORE_SIMULATION_H

#define VERSION_simulation 0.1

#include "dislocation.h"
#include "precision_handler.h"
#include "simulation_data.h"
#include "stress_protocol.h"

#include <memory>
#include <sstream>
#include <vector>

#ifdef BUILD_PYTHON_BINDINGS
#include <boost/python.hpp>
#endif

namespace sdddstCore {

    class Simulation
    {
    public:
        Simulation(std::shared_ptr<SimulationData> sD);
        ~Simulation();

        void integrate(double stepsize, std::vector<Dislocation>& newDislocation, const std::vector<Dislocation>& old, bool useSpeed2, bool calculateInitSpeed, StressProtocolStepType origin, StressProtocolStepType end);

        // calculates the forces (therefore, the speed too) between all d-d and d-p (d: dislocation, p: fixed point defect)
        void calculateSpeeds(const std::vector<Dislocation>& dis, std::vector<double>& res) const;
        void calculateG(double stepsize, const std::vector<Dislocation>& newDislocation, const std::vector<Dislocation>& old, bool useSpeed2, bool calculateInitSpeed, bool useInitSpeedForFirstStep, StressProtocolStepType origin, StressProtocolStepType end) const;
        void calculateJacobian(double stepsize, const std::vector<Dislocation>& data);
        void calculateXError();

        void calculateSparseFormForJacobian();
        void solveEQSys();

        double calculateOrderParameter(const std::vector<double>& speeds) const;
        double calculateStrainIncrement(const std::vector<Dislocation>& old, const std::vector<Dislocation>& newD) const;

        double getElement(int j, int si, int ei) const;

        double getSimTime() const;

        void run();

        void stepStageI();
        void stepStageII();
        void stepStageIII();

#ifdef BUILD_PYTHON_BINDINGS
        // const std::vector<Dislocation>& getStoredDislocationData();
        static Simulation* create(boost::python::object simulationData);
#endif

    private:
        bool succesfulStep;
        double lastWriteTimeFinished;
        double startTime;
        bool initSpeedCalculationIsNeeded;
        bool firstStepRequest;
        double energy;

        std::shared_ptr<SimulationData> sD;
        std::unique_ptr<PrecisionHandler> pH;
    };

}

#endif
