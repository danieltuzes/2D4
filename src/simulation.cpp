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

#include "simulation.h"

#ifdef BUILD_PYTHON_BINDINGS
#include "simulation_data_wrapper.h"
#endif

#include <umfpack.h>
#define _USE_MATH_DEFINES

#include <iomanip>
#include <numeric>
#include <sstream>
#include <math.h>
#include <ctime>
#include <chrono>

#pragma region utilities
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

inline double E(double x, double y, double K)
{
    return exp(-K * (X2(x) + X2(y)));
}

inline double X_dx(double x)
{
    return cos(2 * M_PI * x);
}

inline double X2_dx(double x)
{
    return sin(2 * M_PI * x) / M_PI;
}

inline double E_dx(double x, double y, double K)
{
    return -E(x, y, K) * K * X2_dx(x);
}


// From:
// https://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows
double get_wall_time()
{
    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_start_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(t_start);
    auto t_start_se = t_start_ms.time_since_epoch();

    double time_in_ms = static_cast<double>(t_start_se.count());

    return time_in_ms / 1000.;
}

#pragma endregion


using namespace sdddstCore;

Simulation::Simulation(std::shared_ptr<SimulationData> sD) :
    succesfulStep(true),
    lastWriteTimeFinished(0),
    initSpeedCalculationIsNeeded(true),
    firstStepRequest(true),
    energy(0),
    sD(sD),
    pH(new PrecisionHandler)
{
    // Format setting
    sD->standardOutputLog << std::scientific << std::setprecision(16);

    pH->setMinPrecisity(sD->prec);
    pH->setSize(sD->dc);
}

Simulation::~Simulation()
{
}

/**
    @brief integrate:	evolve the dislocation system in time
    @param stepsize:	how large time step should be made
    @param newDislocation:	the suggested new dislocation configuration will be stored here; wont't be in the range of [-0.5:0.5)
    @param old:	the input dislocation configuration, no need to be in the range of [-0.5:0.5)
    @param useSpeed2:
    @param calculateInitSpeed:
    @param origin:
    @param end:
*/
void Simulation::integrate(double stepsize, std::vector<Dislocation>& newDislocation, const std::vector<Dislocation>& old,
    bool useSpeed2, bool calculateInitSpeed, StressProtocolStepType origin, StressProtocolStepType end)
{
    calculateJacobian(stepsize, newDislocation);
    calculateSparseFormForJacobian();
    for (size_t i = 0; i < sD->ic; i++)
    {
        if (i > 0)
            calculateG(stepsize, newDislocation, old, useSpeed2, false, false, origin, end);
        else
            calculateG(stepsize, newDislocation, old, useSpeed2, calculateInitSpeed, false, origin, end);

        solveEQSys();
        for (size_t j = 0; j < sD->dc; j++)
            newDislocation[j].x -= sD->x[j];
    }
    umfpack_di_free_numeric(&sD->Numeric);
}

void Simulation::calculateSpeeds(const std::vector<Dislocation>& dis, std::vector<double>& res) const
{
    std::fill(res.begin(), res.end(), 0);

    for (unsigned int i = 0; i < sD->dc; i++) // typically unefficient
    {
        for (unsigned int j = i + 1; j < sD->dc; j++)
        {
            double dx = dis[i].x - dis[j].x;
            normalize(dx);

            double dy = dis[i].y - dis[j].y;
            normalize(dy);

            double tmp = dis[i].b * dis[j].b * sD->tau->xy(dx, dy); // The sign of Burgers vectors can be deduced from i and j

            double r2 = dx * dx + dy * dy;
            pH->updateTolerance(r2, i);
            pH->updateTolerance(r2, j);

            res[i] += tmp;
            res[j] -= tmp;
        }

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
            res[i] -= 2 * sD->A * X(dx) * X(dy) * ((1 - expXY) / rSqr - sD->KASQR * expXY) / rSqr * dis[i].b;

            pH->updateTolerance(rSqr, i);
        }

        res[i] += dis[i].b * sD->externalStressProtocol->getExtStress(sD->currentStressStateType);
    }
}

void Simulation::calculateG(double stepsize, const std::vector<Dislocation>& newDislocation, const std::vector<Dislocation>& old,
    bool useSpeed2, bool calculateInitSpeed, bool useInitSpeedForFirstStep, StressProtocolStepType origin, StressProtocolStepType end) const
{
    std::vector<double>* isp = &(sD->initSpeed);
    std::vector<double>* csp = &(sD->speed);
    if (useSpeed2)
    {
        isp = &(sD->initSpeed2);
        csp = &(sD->speed2);
    }

    if (calculateInitSpeed)
    {
        double t = sD->simTime;
        if (origin == StressProtocolStepType::EndOfFirstSmallStep)
            t += sD->stepSize * 0.5;


        sD->externalStressProtocol->calcExtStress(t, origin);
        sD->currentStressStateType = origin;
        calculateSpeeds(old, *isp);
    }

    if (useInitSpeedForFirstStep)
        csp = isp;
    else
    {
        double t = sD->simTime + sD->stepSize;
        if (end == StressProtocolStepType::EndOfFirstSmallStep)
            t -= sD->stepSize * 0.5;

        sD->externalStressProtocol->calcExtStress(t, end);
        sD->currentStressStateType = end;
        calculateSpeeds(newDislocation, *csp);
    }

    for (size_t i = 0; i < sD->dc; i++)
        sD->g[i] = newDislocation[i].x - (1 + sD->dVec[i]) * 0.5 * stepsize * (*csp)[i] - old[i].x - (1 - sD->dVec[i]) * 0.5 * stepsize * (*isp)[i];
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

void Simulation::calculateJacobian(double stepsize, const std::vector<Dislocation>& data)
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
        double tmp = 0;
        double dx;
        double dy;
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
                tmp -= data[j].b * (-sD->A * cos(0.2e1 * M_PI * dx) / M_PI * sin(0.2e1 * M_PI * dy) * ((0.1e1 - pow(M_E, -sD->KASQR * ((0.1e1 - cos(0.2e1 * M_PI * dx)) * pow(M_PI, -0.2e1) / 0.2e1 +
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
        sD->Ax[totalElementCounter++] = -tmp * stepsize;

        // Totally new part
        for (unsigned int i = j + 1; i < sD->dc; i++)
        {
            dx = data[i].x - data[j].x;
            normalize(dx);

            dy = data[i].y - data[j].y; // az y értékek sosem változnak meg, ezért ha a diszlokációk rendezve vannak tárolva, akkor csak az ellentétes típusúak között kell maximum 1 kivonást elvégezni a normalize függvényben
            normalize(dy);

            if (pow(sqrt(dx * dx + dy * dy) - sD->cutOff, 2) < 36.8 * sD->cutOffSqr) // 36.8?? why
            {
                double multiplier = 1;
                if (dx * dx + dy * dy > sD->cutOffSqr)
                    multiplier = exp(-pow(sqrt(dx * dx + dy * dy) - sD->cutOff, 2) * sD->onePerCutOffSqr);

                sD->Ai[totalElementCounter] = i;
                sD->Ax[totalElementCounter++] = stepsize * data[i].b * data[j].b * sD->tau->xy_diff_x(dx, dy) * multiplier;
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

void Simulation::calculateXError()
{
    for (unsigned int i = 0; i < sD->dc; i++)
    {
        double tmp = fabs(sD->bigStep[i].x - sD->secondSmall[i].x);
        pH->updateError(tmp, i);
    }
}

double Simulation::calculateOrderParameter(const std::vector<double>& speeds) const
{
    double orderParameter = 0;
    for (size_t i = 0; i < sD->dc; i++)
        orderParameter += sD->dislocations[i].b * speeds[i];

    return orderParameter;
}

double Simulation::calculateStrainIncrement(const std::vector<Dislocation>& old, const std::vector<Dislocation>& newD) const
{
    double ret = 0;
    for (size_t i = 0; i < old.size(); i++)
        ret += newD[i].b * (newD[i].x - old[i].x); // x is not in the range of [-0.5: 0.5), the difference is the real distance

    return ret;
}

void Simulation::run()
{
    while (
        ((sD->isTimeLimit && sD->simTime < sD->timeLimit) || !sD->isTimeLimit) &&
        ((sD->isStrainIncreaseLimit && sD->totalAccumulatedStrainIncrease < sD->totalAccumulatedStrainIncreaseLimit) || !sD->isStrainIncreaseLimit) &&
        ((sD->isStepCountLimit && sD->succesfulSteps < sD->stepCountLimit) || !sD->isStepCountLimit) &&
        ((sD->countAvalanches && sD->avalancheCount < sD->avalancheTriggerLimit) || !sD->countAvalanches)
        )
    {
        stepStageI();
        stepStageII();
        stepStageIII();
    }

    sD->writeDislocationDataToFile(sD->endDislocationConfigurationPath);
}

void Simulation::stepStageI()
{
    sD->currentStressStateType = StressProtocolStepType::Original;
    if (firstStepRequest)
    {
        startTime = get_wall_time();
        lastWriteTimeFinished = get_wall_time();
        sD->externalStressProtocol->calcExtStress(sD->simTime, StressProtocolStepType::Original);
        calculateSpeeds(sD->dislocations, sD->initSpeed);
        initSpeedCalculationIsNeeded = false;
        double sumAvgSp = std::accumulate(sD->initSpeed.begin(), sD->initSpeed.end(), 0., [](double a, double b) {return a + fabs(b); }) / double(sD->dc);
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

        firstStepRequest = false;
    }

    // Reset the variables for the integration
    sD->bigStep = sD->dislocations;


    /////////////////////////////////
    /// Integrating procedure begins

    integrate(sD->stepSize, sD->bigStep, sD->dislocations, false, initSpeedCalculationIsNeeded, StressProtocolStepType::Original, StressProtocolStepType::EndOfBigStep);

    // This can not get before the first integration step
    succesfulStep = false;
}

void Simulation::stepStageII()
{
    sD->firstSmall = sD->dislocations;
    integrate(0.5 * sD->stepSize, sD->firstSmall, sD->dislocations, false, false, StressProtocolStepType::Original, StressProtocolStepType::EndOfFirstSmallStep);
}

void Simulation::stepStageIII()
{
    sD->secondSmall = sD->firstSmall;

    integrate(0.5 * sD->stepSize, sD->secondSmall, sD->firstSmall, true, true, StressProtocolStepType::EndOfFirstSmallStep, StressProtocolStepType::EndOfSecondSmallStep);

    double vsquare1 = std::accumulate(sD->initSpeed.begin(), sD->initSpeed.end(), 0., [](double a, double b) {return a + b * b; });
    double vsquare2 = std::accumulate(sD->initSpeed2.begin(), sD->initSpeed2.end(), 0., [](double a, double b) {return a + b * b; });

    double energyThisStep = (vsquare1 + vsquare2) * 0.5 * sD->stepSize * 0.5;

    calculateXError();

    /// Precision related error handling
    if (pH->getMaxErrorRatioSqr() < 1)
    {
        succesfulStep = true;
        initSpeedCalculationIsNeeded = true;

        if (sD->calculateStrainDuringSimulation)
        {
            sD->totalAccumulatedStrainIncrease += calculateStrainIncrement(sD->dislocations, sD->firstSmall);
            sD->totalAccumulatedStrainIncrease += calculateStrainIncrement(sD->firstSmall, sD->secondSmall);
        }

        sD->dislocations.swap(sD->secondSmall);
        for (size_t i = 0; i < sD->dc; i++)
            normalize(sD->dislocations[i].x);

        sD->simTime += sD->stepSize;
        sD->succesfulSteps++;

        double orderParameter = 0;
        if (sD->orderParameterCalculationIsOn)
            orderParameter = calculateOrderParameter(sD->speed);

        sD->currentStressStateType = StressProtocolStepType::Original;
        sD->externalStressProtocol->calcExtStress(sD->simTime, StressProtocolStepType::Original);
        calculateSpeeds(sD->dislocations, sD->initSpeed);
        initSpeedCalculationIsNeeded = false;
        double sumAvgSp = std::accumulate(sD->initSpeed.begin(), sD->initSpeed.end(), 0., [](double a, double b) {return a + fabs(b); }) / sD->dc;
        double vsquare = std::accumulate(sD->initSpeed.begin(), sD->initSpeed.end(), 0., [](double a, double b) {return a + b * b; });

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

        energyThisStep += (vsquare2 + vsquare) * 0.5 * sD->stepSize * 0.5;
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
            {
                sD->subconfigDistanceCounter++;
            }
        }

        lastWriteTimeFinished = get_wall_time();

    }
    else
        sD->failedSteps++;


    sD->stepSize = pH->getNewStepSize(sD->stepSize);
    pH->reset();

    if (sD->isMaxStepSizeLimit && sD->maxStepSizeLimit < sD->stepSize)
    {
        sD->stepSize = sD->maxStepSizeLimit;
    }
}


#ifdef BUILD_PYTHON_BINDINGS
/*const std::vector<Dislocation>& Simulation::getStoredDislocationData()
{
    return sD->dislocations;
}*/

Simulation* Simulation::create(boost::python::object simulationData)
{
    boost::python::extract<PySdddstCore::PySimulationData&> x(simulationData);
    if (x.check())
    {
        PySdddstCore::PySimulationData& tmp = x();
        return new Simulation{ tmp.get() };
    }
    return nullptr;
}
#endif
