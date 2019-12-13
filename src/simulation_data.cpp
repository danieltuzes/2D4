// 
// simulation_data.cpp : contains the function definitions for simulation_data.h

#include "simulation_data.h"
#include "constants.h"
#include "stress_protocol.h"

#include <type_traits>

using namespace sdddstCore;

SimulationData::SimulationData(const std::string& startDislocationConfigurationPath, const std::string& fixpointsDataFilePath) :
    Ap(nullptr),
    Ai(nullptr),
    Ax(nullptr),
    x(nullptr),
    null(NULL),
    Symbolic(nullptr),
    Numeric(nullptr),
    cutOffMultiplier(0),
    cutOff(0),
    cutOffSqr(cutOff* cutOff),
    onePerCutOffSqr(1. / cutOffSqr),
    heavisideCutoff(false),
    currentStorageSize(0),
    prec(0),
    weightFunc('c'),
    pc(0),
    dc(0),
    timeLimit(0),
    stepSize(0),
    initStepSize(0),
    stepSizeBeforeWriteout(0),
    simTime(0),
    KASQR(DEFAULT_KASQR),
    A(DEFAULT_A),
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
    subConfigTimes(0),
    subConfigTimesType('a'),
    currentStressStateType(StressProtocolStepType::Original),
    speedThresholdForCutoffChange(0),
    isSpeedThresholdForCutoffChange(false)
{
    readDislocationDataFromFile(startDislocationConfigurationPath);
    readPointDefectDataFromFile(fixpointsDataFilePath);
    initSimulationVariables();
}

#pragma region utility functions
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
            std::cerr << "Error in " << dislocationDataFilePath << ". Burgers' vector supposed to be an integer, -1 or 1, but value " << b << " is found. Program terminates." << std::endl;
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

    disl_order.resize(dc);
    for (size_t i = 0; i < dc; ++i)
    {
        disl_sorted.emplace_back(dislocs_w_id[i].x, dislocs_w_id[i].y);
        disl_order[dislocs_w_id[i].id] = i;
    }

    // update memory usage according to dislocation count
    g.resize(dc);
    initSpeed.resize(dc);
    initSpeed2.resize(dc);
    speed.resize(dc);
    speed2.resize(dc);
    dVec.resize(dc);
    bigStep_sorted.resize(dc);
    firstSmall_sorted.resize(dc);
    secondSmall_sorted.resize(dc);
    
    currentStorageSize = dc;
    Ap = (int*)calloc(size_t(dc) + 1, sizeof(int));
    Ai = (int*)calloc(dc, sizeof(int));
    Ax = (double*)calloc(dc, sizeof(double));
    x = (double*)calloc(dc, sizeof(double));
    if (Ap == NULL)
        throw std::runtime_error("Memory allocation for Ap failed. Program terminates.");

    if (Ai == NULL)
        throw std::runtime_error("Memory allocation for Ai failed. Program terminates.");

    if (Ax == NULL)
        throw std::runtime_error("Memory allocation for Ax failed. Program terminates.");

    if (x == NULL)
        throw std::runtime_error("Memory allocation for x failed. Program terminates.");

    indexes.resize(dc);
}

// increases the reserved size of Ax and Ai by dc and increases currentStorageSize by dc
void SimulationData::increaseCurrentStorageSize(int lastUsedSize)
{
    currentStorageSize += dc;
    double* tmp_Ax = (double*)(realloc(Ax, currentStorageSize * sizeof(double)));
    int* tmp_Ai = (int*)(realloc(Ai, currentStorageSize * sizeof(int)));
    if (tmp_Ax == nullptr || tmp_Ai == nullptr)
    {
        std::cerr << "Cannot realloc memory for Ax and Ai. Program terminates." << std::endl;
        exit(-4);
    }
    Ax = tmp_Ax;
    Ai = tmp_Ai;

    for (unsigned int i = lastUsedSize; i < currentStorageSize; i++)
    {
        Ax[i] = 0;
        Ai[i] = 0;
    }
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
        ofile
        << disl_sorted[id] << "\t"
        << b(id) << "\n";
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
    if (!heavisideCutoff && cutOffMultiplier >= sqrt(2 * dc) / 12)
        multiplier = 1e20;

    cutOff = multiplier / sqrt(dc);
    cutOffSqr = cutOff * cutOff;
    onePerCutOffSqr = 1 / cutOffSqr;
}

// returns the Burgers' vector type based on the index value ID; first half: +1; second half: -1
int SimulationData::b(unsigned int ID) const
{
    if (ID < dc / 2)
        return 1;
    else return -1;
}

// returns true if Burgers' vector for the IDth dislocation is positive; false otherwise
bool SimulationData::is_pos_b(unsigned int ID) const
{
    return ID < dc / 2;
}

#pragma endregion