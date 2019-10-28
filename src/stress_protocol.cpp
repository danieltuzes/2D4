
// stress_protocol.cpp : this file contains the different methods how the external stress should be chosen for a given configuration at a given time

#include "stress_protocol.h"

#pragma region StressProtocol

sdddstCore::StressProtocol::StressProtocol(double initExtStress) : initExtStress(initExtStress) {}

sdddstCore::StressProtocol::StressProtocol() : StressProtocol(0) {}

double sdddstCore::StressProtocol::extStress(double simulationTime) const
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

#pragma region PeriodicLoadProtocol

sdddstCore::CyclicLoadProtocol::CyclicLoadProtocol(double initExtStress, double stressRate, double timePeriod) :
    StressProtocol(initExtStress),
    m_rate(stressRate),
    m_timePeriod(timePeriod),
    stressValues{ initExtStress, initExtStress, initExtStress, initExtStress } {}

sdddstCore::CyclicLoadProtocol::CyclicLoadProtocol() :
    CyclicLoadProtocol(0, 0, 0) {}

double sdddstCore::CyclicLoadProtocol::extStress(double simulationTime) const
{
    double periodicTime = simulationTime - std::floor((simulationTime / m_timePeriod + 1. / 4)) * m_timePeriod;

    double stress = m_rate * periodicTime; // if the derivative is positive
    if (periodicTime > m_timePeriod / 4)
        stress = m_rate * m_timePeriod / 2 - stress; // otherwise, if it is negative

    return stress;
}

#pragma endregion