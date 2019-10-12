
// stress_protocol.cpp : this file contains the different methods how the external stress should be chosen for a given configuration at a given time

#include "stress_protocol.h"

#pragma region StressProtocol

sdddstCore::StressProtocol::StressProtocol(double initExtStress) : initExtStress(initExtStress) {}

sdddstCore::StressProtocol::StressProtocol() : StressProtocol(0) {}

void sdddstCore::StressProtocol::calcExtStress(double, StressProtocolStepType) {}

double sdddstCore::StressProtocol::getExtStress(StressProtocolStepType) const
{
    return initExtStress;
}

#pragma endregion

#pragma region FixedRateProtocol

sdddstCore::FixedRateProtocol::FixedRateProtocol(double initExtStress, double stressRate) :
    StressProtocol(initExtStress),
    m_rate(stressRate),
    stressValues{ initExtStress, initExtStress, initExtStress, initExtStress } {}

sdddstCore::FixedRateProtocol::FixedRateProtocol() :
    FixedRateProtocol(0, 0) {}

void sdddstCore::FixedRateProtocol::calcExtStress(double simulationTime, StressProtocolStepType type)
{
    double value = simulationTime * m_rate + initExtStress;
    switch (type)
    {
    case StressProtocolStepType::Original:
        stressValues[0] = value;
        break;
    case StressProtocolStepType::EndOfBigStep:
        stressValues[1] = value;
        break;
    case StressProtocolStepType::EndOfFirstSmallStep:
        stressValues[2] = value;
        break;
    case StressProtocolStepType::EndOfSecondSmallStep:
        stressValues[3] = value;
        break;
    }
}

double sdddstCore::FixedRateProtocol::getExtStress(StressProtocolStepType type) const
{
    switch (type)
    {
    case StressProtocolStepType::Original:
        return stressValues[0];
    case StressProtocolStepType::EndOfBigStep:
        return stressValues[1];
    case StressProtocolStepType::EndOfFirstSmallStep:
        return stressValues[2];
    case StressProtocolStepType::EndOfSecondSmallStep:
        return stressValues[3];
    }
    return 0;
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

void sdddstCore::CyclicLoadProtocol::calcExtStress(double simulationTime, StressProtocolStepType type)
{
    double periodicTime = simulationTime - std::floor((simulationTime / m_timePeriod + 1. / 4)) * m_timePeriod;

    double value = m_rate * periodicTime;
    if (periodicTime > m_timePeriod / 4)
        value = m_rate * m_timePeriod / 2 - value;

    switch (type)
    {
    case StressProtocolStepType::Original:
        stressValues[0] = value;
        break;
    case StressProtocolStepType::EndOfBigStep:
        stressValues[1] = value;
        break;
    case StressProtocolStepType::EndOfFirstSmallStep:
        stressValues[2] = value;
        break;
    case StressProtocolStepType::EndOfSecondSmallStep:
        stressValues[3] = value;
        break;
    }
}

double sdddstCore::CyclicLoadProtocol::getExtStress(StressProtocolStepType type) const
{
    switch (type)
    {
    case StressProtocolStepType::Original:
        return stressValues[0];
    case StressProtocolStepType::EndOfBigStep:
        return stressValues[1];
    case StressProtocolStepType::EndOfFirstSmallStep:
        return stressValues[2];
    case StressProtocolStepType::EndOfSecondSmallStep:
        return stressValues[3];
    }
    return 0;
}

#pragma endregion