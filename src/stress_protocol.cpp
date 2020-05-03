
// stress_protocol.cpp : this file contains the different methods how the external stress should be chosen for a given configuration at a given time

#include "stress_protocol.h"

#pragma region StressProtocol

sdddstCore::StressProtocol::StressProtocol(double initExtStress) : initExtStress(initExtStress) {}

sdddstCore::StressProtocol::StressProtocol() : StressProtocol(0) {}

// returs the external stress value as a function of the argument simTime
double sdddstCore::StressProtocol::extStress(double) const
{
    return 0;
}

#pragma endregion

#pragma region FixedRateProtocol

sdddstCore::FixedRateProtocol::FixedRateProtocol(double initExtStress, double stressRate) :
    StressProtocol(initExtStress),
    m_rate(stressRate),
    stressValues{ initExtStress, initExtStress, initExtStress, initExtStress } {}

sdddstCore::FixedRateProtocol::FixedRateProtocol() :
    FixedRateProtocol(0, 0) {}

double sdddstCore::FixedRateProtocol::extStress(double simulationTime) const
{
    return simulationTime * m_rate + initExtStress;
}


#pragma endregion

#pragma region CyclicLoadProtocol

sdddstCore::CyclicLoadProtocol::CyclicLoadProtocol(double initExtStress, double stressRate, double timePeriod) :
    FixedRateProtocol(initExtStress, stressRate),
    m_timePeriod(timePeriod) {}

sdddstCore::CyclicLoadProtocol::CyclicLoadProtocol() :
    CyclicLoadProtocol(0, 0, 0) {}

double sdddstCore::CyclicLoadProtocol::extStress(double simulationTime) const
{
    double periodicTime = simulationTime - std::floor((simulationTime / m_timePeriod + 1. / 4)) * m_timePeriod;

    double stress = m_rate * periodicTime + initExtStress; // if the derivative is positive
    if (periodicTime > m_timePeriod / 4)
        stress = m_rate * m_timePeriod / 2 - stress + initExtStress; // otherwise, if it is negative

    return stress;
}

#pragma endregion


#pragma region IncreasingCyclicLoadProtocol

sdddstCore::IncreasingCyclicLoadProtocol::IncreasingCyclicLoadProtocol(double initExtStress, double stressRate, double timePeriod, double amplitudeRate) :
    CyclicLoadProtocol(initExtStress, stressRate, timePeriod),
    m_amplitudeRate(amplitudeRate) {}

sdddstCore::IncreasingCyclicLoadProtocol::IncreasingCyclicLoadProtocol() :
    IncreasingCyclicLoadProtocol(0, 0, 0, 0) {}

double sdddstCore::IncreasingCyclicLoadProtocol::extStress(double simulationTime) const
{
    return CyclicLoadProtocol::extStress(simulationTime) * m_amplitudeRate * simulationTime;
}

#pragma endregion