// 
// stress_protocol.h : contains the function declaration for stress_protocol.cpp, also included in project_parser.cpp

/*
# 0.3
* IncreasingCyclicLoadProtocol has been added; CyclicLoadProtocol was already here but not version-tracked into 0.2, sorry
* CyclicLoadProtocol can take the initExtStress as a constant bias
* CyclicLoadProtocol is derived from FixedRateProtocol
* private members are protected only from now on, so that no need to define the same vars

# 0.2
FixedRateProtocol has been dramatically simplified. extStress is the only non CTOR/DTOR function returning the stress value at a given simulation time

# 0.1 
First version tracked source
*/

#ifndef SDDDST_CORE_STRESS_PROTOCOL_H
#define SDDDST_CORE_STRESS_PROTOCOL_H

#define VERSION_stress_protocol 0.3

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
        
        // returns the external stress at a given time
        virtual double extStress(double simulationTime) const;

    protected:
        double initExtStress;
    };

    // the external stress is increasing linearly with time without limitation
    class FixedRateProtocol : public StressProtocol
    {
    public:
        FixedRateProtocol(double initExtStress, double stressRate);
        FixedRateProtocol();

        // returns the external stress at a given time
        virtual double extStress(double simulationTime) const;

    protected:
        double m_rate;
        double stressValues[4];
    };

    // the external stress is periodic in time with a given time period, linearly increasing and decreasing in a triangular shaped function centered to initExtStress
    class CyclicLoadProtocol : public FixedRateProtocol
    {
    public:
        CyclicLoadProtocol(double initExtStress, double stressRate, double period);
        CyclicLoadProtocol();

        // returns the external stress at a given time
        virtual double extStress(double simulationTime) const;

    protected:
        double m_timePeriod;
    };

    // the external stress is periodic in time with a given time period
    // linearly increasing and decreasing in a triangular shaped function centered to initExtStress
    // and the time period is linearly increasing with time (with a given rate it also means larger amplitude)
    class IncreasingCyclicLoadProtocol : public CyclicLoadProtocol
    {
    public:
        IncreasingCyclicLoadProtocol(double initExtStress, double stressRate, double period, double amplitudeRate);
        IncreasingCyclicLoadProtocol();

        // returns the external stress at a given time
        virtual double extStress(double simulationTime) const;

    protected:
        double m_amplitudeRate;
    };
}

#endif
