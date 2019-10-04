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

#include "simulation_data.h"
#include "constants.h"
#include "StressProtocols/stress_protocol.h"


#ifdef BUILD_PYTHON_BINDINGS
#include "field_wrapper.h"
#include "periodic_shear_stress_ELTE_wrapper.h"
#include "analytic_field_wrapper.h"
#include "stress_protocol_wrapper.h"
#endif

using namespace sdddstCore;

SimulationData::SimulationData(const std::string& startDislocationConfigurationPath, const std::string& fixpointsDataFilePath) :
    cutOffMultiplier(DEFAULT_CUTOFF_MULTIPLIER),
    cutOff(DEFAULT_CUTOFF),
    cutOffSqr(cutOff* cutOff),
    onePerCutOffSqr(1. / cutOffSqr),
    prec(DEFAULT_PRECISION),
    pc(0),
    dc(0),
    ic(DEFAULT_ITERATION_COUNT),
    timeLimit(DEFAULT_TIME_LIMIT),
    stepSize(DEFAULT_STEP_SIZE),
    simTime(DEFAULT_SIM_TIME),
    KASQR(DEFAULT_KASQR),
    A(DEFAULT_A),
    tau(nullptr),
    Ap(nullptr),
    Ai(nullptr),
    Ax(nullptr),
    x(nullptr),
    null(NULL),
    Symbolic(nullptr),
    Numeric(nullptr),
    succesfulSteps(0),
    failedSteps(0),
    totalAccumulatedStrainIncrease(0),
    isStrainIncreaseLimit(false),
    totalAccumulatedStrainIncreaseLimit(0),
    isMaxStepSizeLimit(false),
    maxStepSizeLimit(0),
    isTimeLimit(false),
    isStepCountLimit(false),
    stepCountLimit(0),
    calculateStrainDuringSimulation(false),
    orderParameterCalculationIsOn(false),
    standardOutputLog(),
    startDislocationConfigurationPath(startDislocationConfigurationPath),
    endDislocationConfigurationPath(""),
    externalStressProtocol(nullptr),
    countAvalanches(false),
    avalancheSpeedThreshold(0),
    avalancheTriggerLimit(0),
    avalancheCount(0),
    inAvalanche(false),
    isSaveSubConfigs(false),
    subConfigPath(""),
    subConfigDelay(0),
    subConfigDelayDuringAvalanche(0),
    subconfigDistanceCounter(0),
    currentStressStateType(StressProtocolStepType::Original),
    speedThresholdForCutoffChange(0),
    isSpeedThresholdForCutoffChange(false)
{
    readDislocationDataFromFile(startDislocationConfigurationPath);
    readPointDefectDataFromFile(fixpointsDataFilePath);
    initSimulationVariables();
}

void SimulationData::readDislocationDataFromFile(std::string dislocationDataFilePath)
{
    std::ifstream ifile(dislocationDataFilePath);
    if (!ifile)
    {
        std::cerr << "Cannot open dislocation data file " << dislocationDataFilePath << ". Program terminates." << std::endl;
        assert(ifile.is_open() && "Cannot open the data file to read!"); // if debug mode is active, program will print this too; good to show whether the program is in debug mode
        exit(-1);
    }

    std::vector<DislwId> dislocs_w_id;
    double sum_b = 0; // the sum of the Burger's vector, it must be 0
    // Reading in dislocations
    for (double x; ifile >> x;)
    {
        double y, b;
        if (!(ifile >> y && ifile >> b))
        {
            std::cerr << "Error in " << dislocationDataFilePath << ". Cannot read in the y coordinate and Burger's vector for an x coordinate with value ~ " << x << ". Program terminates." << std::endl;
            exit(-1);
        }

        if (fabs(b - rint(b)) > 1e-5)
        {
            std::cerr << "Error in " << dislocationDataFilePath << ". Burger's vector supposed to be an integer, -1 or 1, but value " << b << " is found. Program terminates." << std::endl;
            exit(-1);
        }

        dislocs_w_id.emplace_back(x, y, b, dislocs_w_id.size());
        sum_b += b;
        dc++;
    }
    
    if (sum_b)
    {
        std::cerr << "Error in " << dislocationDataFilePath << ". The sum of the Burger's vector supposed to be 0, but it is " << sum_b << ". Program terminates." << std::endl;
        exit(-1);
    }

    // order dislocations
    std::sort(dislocs_w_id.begin(), dislocs_w_id.end(), [](const DislwId& l, const DislwId& r) {return l.y + l.b > r.y + r.b; });
    for (const auto& disloc_w_id : dislocs_w_id)
        disl_sorted.emplace_back(disloc_w_id.x, disloc_w_id.y);

    disl_order.resize(dislocs_w_id.size());
    for (size_t i = 0; i < dislocs_w_id.size(); ++i)
    {
        dislocations.emplace_back(dislocs_w_id[i].x, dislocs_w_id[i].y, dislocs_w_id[i].b);
        disl_order[dislocs_w_id[i].id] = i;
    }

    // update memory usage according to dislocation count
    g.resize(dc);
    initSpeed.resize(dc);
    initSpeed2.resize(dc);
    speed.resize(dc);
    speed2.resize(dc);
    dVec.resize(dc);
    bigStep.resize(dc);
    firstSmall.resize(dc);
    secondSmall.resize(dc);
    Ap = (int*)calloc(size_t(dc) + 1, sizeof(int));
    Ai = (int*)calloc(size_t(dc) * dc, sizeof(int));
    Ax = (double*)calloc(size_t(dc) * dc, sizeof(double));
    x = (double*)calloc(dc, sizeof(double));
    assert(Ap && "Memory allocation for Ap failed!");
    assert(Ai && "Memory allocation for Ai failed!");
    assert(Ax && "Memory allocation for Ax failed!");
    assert(x && "Memory allocation for x failed!");
    indexes.resize(dc);
}

void SimulationData::writeDislocationDataToFile(std::string dislocationDataFilePath) const
{
    std::ofstream ofile(dislocationDataFilePath);
    if (!ofile)
    {
        std::cerr << "Warning: the program was unable to create the file " << dislocationDataFilePath << ".\n";
        return;
    }

    ofile << std::setprecision(14);
    for (auto id : disl_order)
        ofile << dislocations[id].x << "\t"
        << dislocations[id].y << "\t"
        << dislocations[id].b << "\n";
}

void SimulationData::readPointDefectDataFromFile(std::string pointDefectDataFilePath)
{
    if (pointDefectDataFilePath.empty())
    {
        pc = 0;
        return;
    }
    pc = 0;
    points.resize(0);
    std::ifstream in(pointDefectDataFilePath);
    assert(in.is_open() && "Cannot open fixpoints data file!"); // why? miért csak debugban?

    // Iterating through the file
    while (!in.eof())
    {
        std::string data;
        in >> data;
        if (data == "")
        {
            break;
        }
        PointDefect tmp;
        tmp.x = std::stod(data);
        in >> tmp.y;
        points.push_back(tmp);
        pc++;
    }
}

void SimulationData::writePointDefectDataToFile(std::string pointDefectDataFilePath) const
{
    std::ofstream out(pointDefectDataFilePath);
    assert(out.is_open() && "Cannot open the data file to write!"); // why? miért csak debugban?
    out << std::scientific << std::setprecision(16);
    for (auto& i : points)
        out << i.x << " " << i.y << "\n";
}

void SimulationData::initSimulationVariables()
{
    updateCutOff();
}

void SimulationData::updateCutOff()
{
    double multiplier = cutOffMultiplier;
    if (cutOffMultiplier >= sqrt(2 * dc) / 12)
        multiplier = 1e20;

    cutOff = multiplier / sqrt(dc);
    cutOffSqr = cutOff * cutOff;
    onePerCutOffSqr = 1 / cutOffSqr;
}

#ifdef BUILD_PYTHON_BINDINGS

Field const& SimulationData::getField()
{
    return *tau;
}

void SimulationData::setField(boost::python::object field)
{
    boost::python::extract<PySdddstCore::PyField&> x(field);
    if (x.check())
    {
        PySdddstCore::PyField& tmp = x();
        tau.reset(tmp.release());
    }
}

const StressProtocol& SimulationData::getStressProtocol()
{
    return *externalStressProtocol;
}

void SimulationData::setStressProtocol(boost::python::object protocol)
{
    boost::python::extract<PySdddstCore::PyStressProtocol&> x(protocol);
    if (x.check())
    {
        PySdddstCore::PyStressProtocol& tmp = x();
        externalStressProtocol.reset(tmp.release());
    }
}
/*
void SimulationData::deleteDislocationCountRelatedData()
{
    dc = 0;
    dislocations.resize(0);
    g.resize(0);
    initSpeed.resize(0);
    initSpeed2.resize(0);
    speed.resize(0);
    speed2.resize(0);
    dVec.resize(0);
    bigStep.resize(0);
    firstSmall.resize(0);
    secondSmall.resize(0);
    free(Ap);
    Ap = nullptr;
    free(Ai);
    Ai = nullptr;
    free(Ax);
    Ax = nullptr;
    free(x);
    x = nullptr;
    indexes.resize(0);
}*/
#endif