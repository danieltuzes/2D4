
// stress_protocol.h : contains the function declaration for stress_protocol.cpp, also included in project_parser.cpp

#ifndef SDDDST_CORE_STRESS_PROTOCOL_H
#define SDDDST_CORE_STRESS_PROTOCOL_H

#define VERSION_stress_protocol 0.1

#include "dislocation.h"

#include <memory>
#include <vector>
#include <cmath>

namespace sdddstCore {

    enum class StressProtocolStepType
    {
        Original,
        EndOfBigStep,
        EndOfFirstSmallStep,
        EndOfSecondSmallStep
    };

    /**
     * @brief The StressProtocol class is the base class for all external stress protocol classes,
     * all parameters for calculation is used from a SimulationData object or should be preset
     */
    class StressProtocol
    {
    public:
        StressProtocol(double initExtStress);
        StressProtocol();

        /**
         * @brief calculateStress calculates the stress value for the given situation
         * @param simulationTime
         * @param dislocations
         * @param type
         */
        virtual void calcExtStress(double simulationTime, StressProtocolStepType type);

        /**
         * @brief getStress returns with the stress value at the given situation
         * @param type
         * @return
         */
        virtual double getExtStress(StressProtocolStepType type) const;

    protected:
        double initExtStress;
    };

    // the external stress is increasing linearly with time without limitation
    class FixedRateProtocol : public StressProtocol
    {
    public:
        FixedRateProtocol(double initExtStress, double stressRate);
        FixedRateProtocol();

        virtual void calcExtStress(double simulationTime, StressProtocolStepType type);
        virtual double getExtStress(StressProtocolStepType type) const;

    private:
        double m_rate;
        double stressValues[4];
    };

    // the external stress is periodic in time with a given time period, linearly increasing and decreasing in a triangular shaped function centered to initExtStress
    class CyclicLoadProtocol : public StressProtocol
    {
    public:
        CyclicLoadProtocol(double initExtStress, double stressRate, double period);
        CyclicLoadProtocol();

        virtual void calcExtStress(double simulationTime, StressProtocolStepType type);
        virtual double getExtStress(StressProtocolStepType type) const;

    private:
        double m_rate;
        double m_timePeriod;
        double stressValues[4];
    };
}

#endif
