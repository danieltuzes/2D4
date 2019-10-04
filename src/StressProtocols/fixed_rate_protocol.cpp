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

#include "StressProtocols/fixed_rate_protocol.h"

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
    double periodicTime = simulationTime - floor((simulationTime / m_timePeriod + 1. / 4)) * m_timePeriod;

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