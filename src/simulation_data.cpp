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
    currentStorageSize(0),
    prec(0),
    pc(0),
    dc(0),
    timeLimit(0),
    stepSize(0),
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
        ofile << disl_sorted[id].x << "\t"
        << disl_sorted[id].y << "\t"
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
    if (cutOffMultiplier >= sqrt(2 * dc) / 12)
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

#ifdef DEBUG_VERSION

void SimulationData::printOut(std::string fname, const std::vector<double>& m_vector) const
{
    std::ofstream of(fname, std::ios_base::app);
    if (!of)
        std::cerr << "Cannot create or append to file " << fname << std::endl;
    else
    {
        for (size_t i = 0; i < m_vector.size(); ++i)
            of << m_vector[i] << "\n";
        std::cout << fname << " is created for debugging." << std::endl;
    }
}

// prints out the selected container's x values to fname        
void SimulationData::printOut(std::string fname, const std::vector<DislwoB>& m_vector) const
{
    std::ofstream of(fname, std::ios_base::app);
    if (!of)
        std::cerr << "Cannot create or append to file " << fname << std::endl;
    else
    {
        for (size_t i = 0; i < m_vector.size(); ++i)
            of << m_vector[i].x << "\n";
        std::cout << fname << " is created for debugging." << std::endl;
    }
}

// prints out size number of elements from the selected array's values to fname        
void SimulationData::printOut(std::string fname, double* array, int size) const
{
    std::ofstream of(fname, std::ios_base::app);
    if (!of)
        std::cerr << "Cannot create or append to file " << fname << std::endl;
    else
    {
        for (int i = 0; i < size; ++i)
            of << array[i] << "\n";
        std::cout << fname << " is created for debugging." << std::endl;
    }
}

// prints out the whole container for vectors and nz number of elements from dynamically allocated arrays to file container name + fname
void SimulationData::printAll(std::string fname, unsigned int nz) const
{
    printOut("speed" + fname, speed);
    printOut("speed2" + fname, speed2);
    printOut("initSpeed" + fname, initSpeed);
    printOut("initSpeed2" + fname, initSpeed2);
    printOut("disl_sorted" + fname, disl_sorted);
    printOut("firstSmall_sorted" + fname, firstSmall_sorted);
    printOut("secondSmall_sorted" + fname, secondSmall_sorted);
    printOut("bigStep_sorted" + fname, bigStep_sorted);
    printOut("dVec" + fname, dVec);
    printOut("g" + fname, g);
    printOut("Ax" + fname, Ax, nz);
    printOut("x" + fname, x, dc);
}

// checks if all values in the container are finite
bool SimulationData::isFinite(std::vector<double> m_vector, double lb, double ub)
{
    for (size_t i = 0; i < m_vector.size(); ++i)
    {
        if (!((lb <= m_vector[i]) && (m_vector[i] <= ub)))
        {
            std::cerr << "value at [" << i << "] is " << m_vector[i] << " is unexpected\n";
            return false;
        }
    }
    return true;
}

// checks if all x coordinate values in the container are finite
bool SimulationData::isFinite(std::vector<DislwoB> m_vector, double lb, double ub)
{
    for (size_t i = 0; i < m_vector.size(); ++i)
    {
        if (!((lb <= m_vector[i].x) && (m_vector[i].x <= ub)))
        {
            std::cerr << "value at [" << i << "] is " << m_vector[i].x << " is unexpected\n";
            return false;
        }
    }
    return true;
}

// checks if the first nz number of elements in the array are finite
bool SimulationData::isFinite(double* m_array, size_t nz, double lb, double ub)
{
    std::vector<double> m_vector(m_array, m_array + nz);
    return isFinite(m_vector, lb, ub);
}

// checks if all the containers and arrays up to nz number of elements contain only finite values and
// print out results to labeled filenames, label should match ^[\w,\s-]+
bool SimulationData::isAllFinite(size_t nz, std::string label)
{
    if (simTime < 8.433e-04)
        return true;

    bool allfinite = true;
    if (!isFinite(speed, -1E18, 1E18))
    {
        std::cerr << "speed can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(speed2, -1E18, 1E18))
    {
        std::cerr << "speed2 can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(initSpeed, -1E18, 1E18))
    {
        std::cerr << "initSpeed can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(initSpeed2, -1E18, 1E18))
    {
        std::cerr << "initSpeed2 can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(disl_sorted, -1, 1))
    {
        std::cerr << "disl_sorted can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(firstSmall_sorted, -1, 1))
    {
        std::cerr << "firstSmall_sorted can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(secondSmall_sorted, -1, 1))
    {
        std::cerr << "secondSmall_sorted can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(bigStep_sorted, -1, 1))
    {
        std::cerr << "bigStep_sorted can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(dVec, 0, 1))
    {
        std::cerr << "dVec can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(g, -1E18, 1E18))
    {
        std::cerr << "g can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(Ax, nz, -(int)dc, dc))
    {
        std::cerr << "Ax can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(x, dc, -1, 1))
    {
        std::cerr << "x can lead to nonfinite result\n";
        allfinite = false;
    }

    if (!allfinite)
    {
        std::cerr
            << "The error was found at label "
            << label << "\t"
            << succesfulSteps << "\t"
            << failedSteps << "\t"
            << stepSize << "\t"
            << simTime << "\n"
            << "Debugging files will be creatred."
            << std::endl;

        std::stringstream ss;
        ss << label << "_" << succesfulSteps << "_" << failedSteps << ".txt";
        printAll(ss.str(), nz);
        std::cerr << std::endl;
    }
    return allfinite;
}

#endif