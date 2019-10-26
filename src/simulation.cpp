//
// simulation.cpp : contains the function definitions from simulation.h; time evolve the dislocation system, choose the stepsize adaptively by calculating the error

#define _USE_MATH_DEFINES
#include "simulation.h"
#include <suitesparse/umfpack.h>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <math.h>
#include <ctime>
#include <chrono>

#pragma region utilities
 // transforms n into the range [-0.5:0.5)
inline void normalize(double& n)
{
    while (n < -0.5) // it shouldn't be smaller than -0.9999999, "while" could be changed to "if", but it isn't faster
        n += 1;

    while (n >= 0.5) // it shouldn't be larger than 0.99999999, "while" could be changed to "if", but it isn't faster
        n -= 1;
}

inline double X(double x)
{
    return sin(2 * M_PI * x) * 0.5 / M_PI;
}

inline double X2(double x)
{
    return (1 - cos(2 * M_PI * x)) * 0.5 / M_PI / M_PI;
}


// calculates the absolute value square of a vector as in mathematics
double absvalsq(const std::vector<double>& input)
{
    return std::accumulate(input.begin(), input.end(), 0., [](double a, double b) {return a + b * b; });
}

// returns wall time in ms
double get_wall_time()
{
    // From https://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows

    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_start_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(t_start);
    auto t_start_se = t_start_ms.time_since_epoch();

    double time_in_ms = static_cast<double>(t_start_se.count());

    return time_in_ms / 1000.;
}

#pragma endregion

using namespace sdddstCore;

Simulation::Simulation(std::shared_ptr<SimulationData> sD) :
    lastWriteTimeFinished(0),
    startTime(0),
    initSpeedCalculationIsNeeded(true),
    energy(0),
    sD(sD),
    pH(new PrecisionHandler)
{
    // Format setting
    sD->standardOutputLog << std::scientific << std::setprecision(16);

    pH->setMinPrecisity(sD->prec);
    pH->setSize(sD->dc);
}

double Simulation::getElement(int j, int si, int ei) const
{
    int len = ei - si;
    if (len > 1)
    {
        int tmp = len / 2;
        double a;

        if (sD->Ai[si + tmp] > j)
        {
            a = getElement(j, si, si + tmp);
            if (a != 0)
                return a;
        }
        else
        {
            a = getElement(j, si + tmp, ei);
            if (a != 0)
                return a;
        }
    }
    else if (sD->Ai[si] == j)
        return sD->Ax[si];

    return 0;
}

double Simulation::getSimTime() const
{
    return sD->simTime;
}

void Simulation::calculateSparseFormForJacobian()
{
    (void)umfpack_di_symbolic(sD->dc, sD->dc, sD->Ap, sD->Ai, sD->Ax, &(sD->Symbolic), sD->null, sD->null);
    (void)umfpack_di_numeric(sD->Ap, sD->Ai, sD->Ax, sD->Symbolic, &(sD->Numeric), sD->null, sD->null);
    umfpack_di_free_symbolic(&(sD->Symbolic));
}

void Simulation::solveEQSys()
{
    (void)umfpack_di_solve(UMFPACK_A, sD->Ap, sD->Ai, sD->Ax, sD->x, sD->g.data(), sD->Numeric, sD->null, sD->null);
}

double Simulation::calculateOrderParameter(const std::vector<double>& speeds) const
{
    double orderParameter = 0;
    for (size_t i = 0; i < sD->dc; i++)
        orderParameter += sD->b(i) * speeds[i];

    return orderParameter;
}

void Simulation::run()
{
    {
        std::time_t start_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::cout << "Simulation started on " << std::ctime(&start_time) << std::endl;
    }
    startTime = get_wall_time();
    lastWriteTimeFinished = get_wall_time();
    sD->externalStressProtocol->calcExtStress(sD->simTime, StressProtocolStepType::Original);
    calculateSpeeds(sD->disl_sorted, sD->initSpeed);
    initSpeedCalculationIsNeeded = false; // if initSpeed corresponds to the actual state of the dislocations; false after a successful step
    double sumAvgSp = std::accumulate(sD->initSpeed.begin(), sD->initSpeed.end(), 0., [](double a, double b) {return a + fabs(b); }) / sD->dc;
    double vsquare = std::accumulate(sD->initSpeed.begin(), sD->initSpeed.end(), 0., [](double a, double b) {return a + b * b; });

    // First two line in the log file
    sD->standardOutputLog << "# simtime\tsuccessfullsteps\tfailedsteps\tmaxErrorRatioSqr\tsumAvgSp\tcutOff\torder parameter\texternal stress\tstageI - III time\tstrain\tvsquare\tenergy\twall_time_elapsed" << std::endl;
    sD->standardOutputLog << sD->simTime << "\t"
        << sD->succesfulSteps << "\t"
        << sD->failedSteps << "\t"
        << 0 << "\t"
        << sumAvgSp << "\t"
        << sD->cutOff << "\t"
        << "-" << "\t"
        << sD->externalStressProtocol->getExtStress(sD->currentStressStateType) << "\t"
        << "-" << "\t"
        << sD->totalAccumulatedStrainIncrease << "\t"
        << vsquare << "\t"
        << energy << "\t"
        << 0 << std::endl;

    while (
        ((sD->isTimeLimit && sD->simTime < sD->timeLimit) || !sD->isTimeLimit) &&
        ((sD->isStrainIncreaseLimit && sD->totalAccumulatedStrainIncrease < sD->totalAccumulatedStrainIncreaseLimit) || !sD->isStrainIncreaseLimit) &&
        ((sD->isStepCountLimit && sD->succesfulSteps < sD->stepCountLimit) || !sD->isStepCountLimit) &&
        ((sD->countAvalanches && sD->avalancheCount < sD->avalancheTriggerLimit) || !sD->countAvalanches)
        )
    {
        /////////////////////////////////
        /// step stage I
        /////////////////////////////////
        {
            sD->currentStressStateType = StressProtocolStepType::Original;
            integrate(sD->stepSize, sD->bigStep_sorted, sD->disl_sorted, false, initSpeedCalculationIsNeeded, StressProtocolStepType::Original, StressProtocolStepType::EndOfBigStep);
        }

        /////////////////////////////////
        /// step stageII
        /////////////////////////////////
        {
            integrate(0.5 * sD->stepSize, sD->firstSmall_sorted, sD->disl_sorted, false, false, StressProtocolStepType::Original, StressProtocolStepType::EndOfFirstSmallStep);
        }

        /////////////////////////////////
        /// step stage III
        /////////////////////////////////
        {
            integrate(0.5 * sD->stepSize, sD->secondSmall_sorted, sD->firstSmall_sorted, true, true, StressProtocolStepType::EndOfFirstSmallStep, StressProtocolStepType::EndOfSecondSmallStep);

            double vsquare2 = absvalsq(sD->initSpeed2);

            double energyThisStep = (absvalsq(sD->initSpeed) + vsquare2) * 0.5 * sD->stepSize * 0.5;

            calculateXError();

            /// Precision related error handling
            if (pH->getMaxErrorRatioSqr() < 1)
            {
                initSpeedCalculationIsNeeded = true;

                if (sD->calculateStrainDuringSimulation)
                {
                    sD->totalAccumulatedStrainIncrease += calculateStrainIncrement(sD->disl_sorted, sD->firstSmall_sorted);
                    sD->totalAccumulatedStrainIncrease += calculateStrainIncrement(sD->firstSmall_sorted, sD->secondSmall_sorted);
                }

                sD->disl_sorted.swap(sD->secondSmall_sorted);
                for (size_t i = 0; i < sD->dc; i++)
                    normalize(sD->disl_sorted[i].x);

                sD->simTime += sD->stepSize;
                sD->succesfulSteps++;

                double orderParameter = 0;
                if (sD->orderParameterCalculationIsOn)
                    orderParameter = calculateOrderParameter(sD->speed);

                sD->currentStressStateType = StressProtocolStepType::Original;
                sD->externalStressProtocol->calcExtStress(sD->simTime, StressProtocolStepType::Original);
                calculateSpeeds(sD->disl_sorted, sD->initSpeed);
                initSpeedCalculationIsNeeded = false;
                double sumAvgSp = std::accumulate(sD->initSpeed.begin(), sD->initSpeed.end(), 0., [](double a, double b) {return a + fabs(b); }) / sD->dc;

                if (sD->countAvalanches)
                {
                    if (sD->inAvalanche && sumAvgSp < sD->avalancheSpeedThreshold)
                    {
                        sD->avalancheCount++;
                        sD->inAvalanche = false;
                    }
                    else if (sumAvgSp > sD->avalancheSpeedThreshold)
                        sD->inAvalanche = true;
                }
                double vsquare = absvalsq(sD->initSpeed);
                energyThisStep += (vsquare + vsquare2) * 0.5 * sD->stepSize * 0.5;
                sD->standardOutputLog << sD->simTime << "\t"
                    << sD->succesfulSteps << "\t"
                    << sD->failedSteps << "\t"
                    << pH->getMaxErrorRatioSqr() << "\t"
                    << sumAvgSp << "\t"
                    << sD->cutOff << "\t";

                if (sD->orderParameterCalculationIsOn)
                    sD->standardOutputLog << orderParameter << "\t";
                else
                    sD->standardOutputLog << "-" << "\t";

                sD->standardOutputLog << sD->externalStressProtocol->getExtStress(StressProtocolStepType::Original) << "\t"
                    << get_wall_time() - lastWriteTimeFinished << "\t";

                if (sD->calculateStrainDuringSimulation)
                    sD->standardOutputLog << sD->totalAccumulatedStrainIncrease << "\t";
                else
                    sD->standardOutputLog << "-" << "\t";


                if (sD->isSpeedThresholdForCutoffChange && sD->speedThresholdForCutoffChange > sumAvgSp)
                {
                    sD->cutOffMultiplier = 1e20;
                    sD->updateCutOff();
                }

                energy += energyThisStep;

                sD->standardOutputLog << vsquare << "\t"
                    << energy << "\t"
                    << get_wall_time() - startTime << std::endl;

                if (sD->isSaveSubConfigs)
                {
                    if ((!sD->inAvalanche && sD->subConfigDelay >= sD->subconfigDistanceCounter) || (sD->inAvalanche && sD->subConfigDelayDuringAvalanche >= sD->subconfigDistanceCounter))
                    {
                        sD->subconfigDistanceCounter = 0;
                        std::stringstream ss;
                        ss << std::setprecision(16);
                        ss << sD->simTime;
                        sD->writeDislocationDataToFile(sD->subConfigPath + ss.str() + ".dconf");
                    }
                    else
                        sD->subconfigDistanceCounter++;
                }

                lastWriteTimeFinished = get_wall_time();
            }
            else
                sD->failedSteps++;

            sD->stepSize = pH->getNewStepSize(sD->stepSize);
            pH->reset();

            if (sD->isMaxStepSizeLimit && sD->maxStepSizeLimit < sD->stepSize)
                sD->stepSize = sD->maxStepSizeLimit;
        }
    }

    sD->writeDislocationDataToFile(sD->endDislocationConfigurationPath);
    std::cout << "Simulation is done (" << get_wall_time() - startTime << " s).\n";
}

void Simulation::calculateXError()
{
    for (unsigned int i = 0; i < sD->dc; i++)
    {
        double tmp = fabs(sD->bigStep_sorted[i].x - sD->secondSmall_sorted[i].x);
        pH->updateError(tmp, i);
    }
}

/**
    @brief integrate:      evolve the dislocation system in time from the actual position t_0 until t_0 + stepsize
    @param stepsize:       how large time step should be made
    @param newDislocation: the suggested new dislocation configuration will be stored here; wont't be in the range of [-0.5:0.5)
    @param old:            the input dislocation configuration, no need to be in the range of [-0.5:0.5)
    @param useSpeed2:      true for the second small step
    @param calcInitSpeed:  if initSpeed can be used as the speeds of the particles; false after a successful step and after the small step
    @param origin:         stress at the beginning of the integration
    @param end:            stress at the end of the integration
*/
void Simulation::integrate(double stepsize, std::vector<DislwoB>& newDisloc, const std::vector<DislwoB>& oldDisloc,
    bool useSpeed2, bool calcInitSpeed, StressProtocolStepType origin, StressProtocolStepType end)
{
    newDisloc = oldDisloc; // initialize newDisloc from oldDisloc
    calculateJacobian(stepsize, oldDisloc);
    calculateSparseFormForJacobian();

    calculateG(stepsize, newDisloc, oldDisloc, useSpeed2, calcInitSpeed, origin, end);
    solveEQSys();
    for (size_t j = 0; j < sD->dc; j++)
        newDisloc[j].x -= sD->x[j];

    calculateG(stepsize, newDisloc, oldDisloc, useSpeed2, false, origin, end);
    solveEQSys();
    for (size_t j = 0; j < sD->dc; j++)
        newDisloc[j].x -= sD->x[j];

    umfpack_di_free_numeric(&sD->Numeric);
}

void Simulation::calculateSpeedsAtStresses(const std::vector<DislwoB>& dis, std::vector<double>& forces_A, std::vector<double>& forces_B, double extStress_A, double extStress_B) const
{
    std::fill(forces_A.begin(), forces_A.end(), 0);

    if (&forces_A != &forces_B)
        std::fill(forces_B.begin(), forces_B.end(), 0);

    for (unsigned int i = 0; i < sD->dc; i++) // typically unefficient
    {
        for (unsigned int j = i + 1; j < sD->dc; j++)
        {
            double dx = dis[i].x - dis[j].x;
            normalize(dx);

            double dy = dis[i].y - dis[j].y;
            normalize(dy);

            double r2 = (dx * dx + dy * dy) * 0.0025; // proportional to distance square
            pH->updateTolerance(r2, i);
            pH->updateTolerance(r2, j);

            double force = sD->tau.xy(dx, dy); // The sign of Burgers vectors can be deduced from i and j

            if ((sD->is_pos_b(i) && sD->is_pos_b(j)) || (!sD->is_pos_b(i) && !sD->is_pos_b(j)))
            {
                forces_A[i] += force;
                forces_A[j] -= force;
            }
            else
            {
                forces_A[i] -= force;
                forces_A[j] += force;
            }

        }

#ifdef USE_POINT_DEFECTS
        // iterate over point defects and calulate their stress contribution
        for (size_t j = 0; j < sD->pc; j++)
        {
            double dx = dis[i].x - sD->points[j].x;
            normalize(dx);

            double dy = dis[i].y - sD->points[j].y;
            normalize(dy);

            double xSqr = X2(dx);
            double ySqr = X2(dy);
            double rSqr = xSqr + ySqr;
            double expXY = exp(-sD->KASQR * rSqr);
            forces_A[i] -= 2 * sD->A * X(dx) * X(dy) * ((1 - expXY) / rSqr - sD->KASQR * expXY) / rSqr * sD->b(i);

            pH->updateTolerance(rSqr, i);
    }
#endif
        if (&forces_A == &forces_B)
            forces_A[i] += extStress_A * sD->b(i);
        else
        {
            forces_B[i] = forces_A[i] + extStress_B * sD->b(i);
            forces_A[i] += extStress_A * sD->b(i);
        }

}
}

void Simulation::calculateSpeedsAtStress(const std::vector<DislwoB>& dis, std::vector<double>& forces, double extStress) const
{
    calculateSpeedsAtStresses(dis, forces, forces, extStress, extStress);
}

void Simulation::calculateSpeeds(const std::vector<DislwoB>& dis, std::vector<double>& forces) const
{
    double extStress = sD->externalStressProtocol->getExtStress(sD->currentStressStateType);
    calculateSpeedsAtStress(dis, forces, extStress);
}

/**
    @brief calculateG:     calculates the g vector
    @param stepsize:       how large time step should be made
    @param newDisloc:      the suggested new dislocation configuration will be stored here; wont't be in the range of [-0.5:0.5)
    @param oldDisloc:      the input dislocation configuration, no need to be in the range of [-0.5:0.5)
    @param useSpeed2:      true for the second small step
    @param calcInitSpeed:  if initSpeed can be used as the speeds of the particles; calcInitSpeed is false after an unsuccessful step, at first small step and at the 2nd NR iteration; true after successful step and at the second step, but both cases only the first NR step
    @param origin:         stress at the beginning of the integration
    @param end:            stress at the end of the integration
*/
void Simulation::calculateG(double stepsize, const std::vector<DislwoB>& newDisloc, const std::vector<DislwoB>& oldDisloc,
    bool useSpeed2, bool calcInitSpeed, StressProtocolStepType origin, StressProtocolStepType end) const
{
    std::vector<double>* isp = &(sD->initSpeed);
    std::vector<double>* csp = &(sD->speed);
    if (useSpeed2)
    {
        isp = &(sD->initSpeed2); // initial speeds
        csp = &(sD->speed2);     // speeds at the end of the step?
    }

    double extStresst_0;
    if (calcInitSpeed)
    {
        // for t_0
        double t_0 = sD->simTime;
        if (origin == StressProtocolStepType::EndOfFirstSmallStep)
            t_0 += sD->stepSize * 0.5;

        sD->externalStressProtocol->calcExtStress(t_0, origin);
        sD->currentStressStateType = origin;
        double extStresst_0 = sD->externalStressProtocol->getExtStress(sD->currentStressStateType);

        // for t_1
        double t_1 = sD->simTime + sD->stepSize; // the time 
        if (end == StressProtocolStepType::EndOfFirstSmallStep)
            t_1 -= sD->stepSize * 0.5;

        sD->externalStressProtocol->calcExtStress(t_1, end);
        sD->currentStressStateType = end;
        double extStresst_1 = sD->externalStressProtocol->getExtStress(sD->currentStressStateType);

        calculateSpeedsAtStresses(oldDisloc, *isp, *csp, extStresst_0, extStresst_1);
    }
    else
    {
        double t_1 = sD->simTime + sD->stepSize; // the time 
        if (end == StressProtocolStepType::EndOfFirstSmallStep)
            t_1 -= sD->stepSize * 0.5;

        sD->externalStressProtocol->calcExtStress(t_1, end);
        sD->currentStressStateType = end;

        calculateSpeeds(newDisloc, *csp);
    }

    for (size_t i = 0; i < sD->dc; i++)
        sD->g[i] = newDisloc[i].x - oldDisloc[i].x - ((1 + sD->dVec[i]) * stepsize * (*csp)[i] + (1 - sD->dVec[i]) * stepsize * (*isp)[i]) / 2;
}

void Simulation::calculateJacobian(double stepsize, const std::vector<DislwoB>& data)
{
    int totalElementCounter = 0;

    for (unsigned int j = 0; j < sD->dc; j++)
    {
        // Previously calculated part
        for (unsigned int i = 0; i < j; i++)
        {
            double v = getElement(j, sD->Ap[i], sD->Ap[i + 1]);
            if (v != 0)
            {
                sD->Ai[totalElementCounter] = i;
                sD->Ax[totalElementCounter++] = v;
            }
        }
        // Add the diagonal element (it will be calculated later and the point defects now)
        sD->Ai[totalElementCounter] = j;

        double dx;
        double dy;
        double tmp = 0;
#ifdef USE_POINT_DEFECTS

        for (size_t l = 0; l < sD->pc; l++)
        {
            dx = data[j].x - sD->points[l].x;
            normalize(dx);
            dy = data[j].y - sD->points[l].y;
            normalize(dy);

            if (pow(sqrt(dx * dx + dy * dy) - sD->cutOff, 2) < 36.8 * sD->cutOffSqr)
            {
                double multiplier = 1;
                if (dx * dx + dy * dy > sD->cutOffSqr)
                {
                    multiplier = exp(-pow(sqrt(dx * dx + dy * dy) - sD->cutOff, 2) * sD->onePerCutOffSqr);
                }
                tmp -= sD->b(j) * (-sD->A * cos(0.2e1 * M_PI * dx) / M_PI * sin(0.2e1 * M_PI * dy) * ((0.1e1 - pow(M_E, -sD->KASQR * ((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) / 0.2e1 +
                    (0.1e1 - cos(0.2e1 * M_PI * dy)) * pow(M_PI, -0.2e1) / 0.2e1))) /
                    ((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) / 0.2e1 + (0.1e1 - cos(0.2e1 * M_PI * dy)) *
                        pow(M_PI, -0.2e1) / 0.2e1) - sD->KASQR * pow(M_E, -sD->KASQR * ((0.1e1 - cos(0.2e1 * M_PI * dx)) *
                            pow(M_PI, -0.2e1) / 0.2e1 +
                            (0.1e1 - cos(0.2e1 * M_PI * dy)) *
                            pow(M_PI, -0.2e1) / 0.2e1))) /
                            ((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) / 0.2e1 + (0.1e1 - cos(0.2e1 * M_PI * dy)) * pow(M_PI, -0.2e1) / 0.2e1)
                    - sD->A * sin(0.2e1 * M_PI * dx) * pow(M_PI, -0.2e1) * sin(0.2e1 * M_PI * dy) * (pow(M_E, -sD->KASQR * ((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) / 0.2e1 +
                    (0.1e1 - cos(0.2e1 * M_PI * dy)) * pow(M_PI, -0.2e1) / 0.2e1)) *
                        sD->KASQR * sin(0.2e1 * M_PI * dx) / M_PI * log(M_E) / ((0.1e1 - cos(0.2e1 * M_PI * dx)) *
                            pow(M_PI, -0.2e1) / 0.2e1 +
                            (0.1e1 - cos(0.2e1 * M_PI * dy)) *
                            pow(M_PI, -0.2e1) / 0.2e1) -
                            (0.1e1 - pow(M_E, -sD->KASQR * ((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) /
                                0.2e1 + (0.1e1 - cos(0.2e1 * M_PI * dy)) *
                                pow(M_PI, -0.2e1) / 0.2e1))) *
                        pow((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) / 0.2e1 +
                        (0.1e1 - cos(0.2e1 * M_PI * dy)) * pow(M_PI, -0.2e1) / 0.2e1, -0.2e1) *
                        sin(0.2e1 * M_PI * dx) / M_PI + sD->KASQR * sD->KASQR * pow(M_E, -sD->KASQR *
                        ((0.1e1 - cos(0.2e1 * M_PI * dx)) *
                            pow(M_PI, -0.2e1) /
                            0.2e1 + (0.1e1 - cos(0.2e1 * M_PI * dy)) *
                            pow(M_PI, -0.2e1) / 0.2e1)) *
                        sin(0.2e1 * M_PI * dx) / M_PI * log(M_E)) / ((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) /
                            0.2e1 + (0.1e1 - cos(0.2e1 * M_PI * dy)) *
                            pow(M_PI, -0.2e1) / 0.2e1) / 0.2e1 + sD->A *
                    pow(sin(0.2e1 * M_PI * dx), 0.2e1) * pow(M_PI, -0.3e1) * sin(0.2e1 * M_PI * dy) * ((0.1e1 - pow(M_E, -sD->KASQR * ((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) / 0.2e1 +
                    (0.1e1 - cos(0.2e1 * M_PI * dy)) * pow(M_PI, -0.2e1) / 0.2e1))) /
                        ((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) / 0.2e1 + (0.1e1 - cos(0.2e1 * M_PI * dy)) *
                            pow(M_PI, -0.2e1) / 0.2e1) - sD->KASQR * pow(M_E, -sD->KASQR * ((0.1e1 - cos(0.2e1 * M_PI * dx)) *
                                pow(M_PI, -0.2e1) / 0.2e1 +
                                (0.1e1 - cos(0.2e1 * M_PI * dy)) *
                                pow(M_PI, -0.2e1) / 0.2e1))) *
                    pow((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) / 0.2e1 + (0.1e1 - cos(0.2e1 * M_PI * dy)) * pow(M_PI, -0.2e1) / 0.2e1, -0.2e1) / 0.2e1) * multiplier;
                // see simulation_tmp_simplified What is this joke? log(M_E)? 0.2e1?
            }
        }

#endif

        sD->Ax[totalElementCounter++] = -tmp * stepsize;
        // Totally new part
        for (unsigned int i = j + 1; i < sD->dc; i++)
        {
            dx = data[i].x - data[j].x;
            normalize(dx);

            dy = data[i].y - data[j].y; // az y értékek konstansok, így ha a diszlokációk rendezve vannak tárolva, akkor csak az ellentétes típusúak között kell max 1 kivonást elvégezni a normalize függvényben
            normalize(dy);

            double distSq = dx * dx + dy * dy;

            if (std::fabs(sqrt(distSq) - sD->cutOff) < 6 * sD->cutOff) // itt a 6-os szorzót simán ki kéne hagyni, meg a következő 3 sort is úgy, ahogy van, TODO
            {
                double multiplier = 1;
                if (distSq > sD->cutOffSqr)
                    multiplier = exp(-pow(sqrt(distSq) - sD->cutOff, 2) * sD->onePerCutOffSqr);

                sD->Ai[totalElementCounter] = i;
                sD->Ax[totalElementCounter++] = stepsize * sD->b(i) * sD->b(j) * sD->tau.xy_diff_x(dx, dy) * multiplier;
            }
        }
        sD->Ap[j + 1] = totalElementCounter;
    }

    for (unsigned int j = 0; j < sD->dc; j++)
    {
        double subSum = 0;
        for (int i = sD->Ap[j]; i < sD->Ap[j + 1]; i++)
        {
            if (sD->Ai[i] == int(j))
                sD->indexes[j] = i;

            subSum -= sD->Ax[i];
        }

        sD->Ax[sD->indexes[j]] = subSum;

        if (subSum > 0) // why
            sD->dVec[j] = std::pow(1 - 1 / (subSum + 1), 2);
        else
            sD->dVec[j] = 0;
    }

    for (unsigned int j = 0; j < sD->dc; j++)
    {
        for (int i = sD->Ap[j]; i < sD->Ap[j + 1]; i++)
            sD->Ax[i] *= (1 + sD->dVec[sD->Ai[i]]) * 0.5;

        sD->Ax[sD->indexes[j]] += 1;
    }
}

double Simulation::calculateStrainIncrement(const std::vector<DislwoB>& old, const std::vector<DislwoB>& newD) const
{
    double ret = 0;
    for (size_t i = 0; i < old.size(); i++)
        ret += sD->b(i) * (newD[i].x - old[i].x); // x is not in the range of [-0.5: 0.5), the difference is the real distance

    return ret;
}
