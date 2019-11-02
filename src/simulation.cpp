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

// returns wall time in seconds
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

// starts the simulation, all required data are present; modifies the state of the system and creates files
void Simulation::run()
{
    ////////////////////////////////////////////////////////
    // before loop, first run
    ////////////////////////////////////////////////////////

    // cout date and time
    {
        std::time_t start_date_and_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::cout << "Simulation started on " << std::ctime(&start_date_and_time) << std::endl;
    }

    double startTime = get_wall_time();   // start time in seconds
    double lastLogTime = get_wall_time(); // the time of the last write into the logfile
    double energy = 0;                    // 

    // First two line in the log file
    {
        calculateSpeedsAtTime(sD->disl_sorted, sD->initSpeed, sD->simTime); // calculated only for the purpose to be able to write out initial speed values

        double sumAvgSp = std::accumulate(sD->initSpeed.begin(), sD->initSpeed.end(), 0., [](double a, double b) {return a + fabs(b); }) / sD->dc;
        double vsquare = std::accumulate(sD->initSpeed.begin(), sD->initSpeed.end(), 0., [](double a, double b) {return a + b * b; });
        double extStress = sD->externalStressProtocol->extStress(sD->simTime);

        sD->standardOutputLog
            << "# simtime\tsuccessfullsteps\tfailedsteps\tmaxErrorRatioSqr\tsumAvgSp\tcutOff\torder parameter\texternal stress\tstageI - III time\tstrain\tvsquare\tenergy\twall_time_elapsed" << std::endl
            << sD->simTime << "\t"
            << sD->succesfulSteps << "\t"
            << sD->failedSteps << "\t"
            << 0 << "\t"
            << sumAvgSp << "\t"
            << sD->cutOff << "\t"
            << "-" << "\t"
            << extStress << "\t"
            << "-" << "\t"
            << sD->totalAccumulatedStrainIncrease << "\t"
            << vsquare << "\t"
            << energy << "\t"
            << 0 << std::endl;
    }

    while (
        ((sD->isTimeLimit && sD->simTime < sD->timeLimit) || !sD->isTimeLimit) &&
        ((sD->isStrainIncreaseLimit && sD->totalAccumulatedStrainIncrease < sD->totalAccumulatedStrainIncreaseLimit) || !sD->isStrainIncreaseLimit) &&
        ((sD->isStepCountLimit && sD->succesfulSteps < sD->stepCountLimit) || !sD->isStepCountLimit) &&
        ((sD->countAvalanches && sD->avalancheCount < sD->avalancheTriggerLimit) || !sD->countAvalanches)
        )
    {
        //////////////////////////////////////////////////////////////////
        /// step stage I: one large step
        //////////////////////////////////////////////////////////////////
        double t_1 = sD->simTime + sD->stepSize;    // simTime at the end of a large step
        double stepSize = sD->stepSize;             // the size of the time step at the actual stage
        {
            calcJacobianAndSpeedsAtTime(stepSize, sD->disl_sorted, sD->speed2, t_1); //## kiszámolja a jakobit sD->disl_sorted-ból stepSize-ra stepSize/2-re; sD->disl_sorted = config[1]

            for (unsigned int i = 0; i < sD->dc; i++)
                sD->g[i] = -stepSize * ((1 + sD->dVec[i]) * sD->speed2[i] + (1 - sD->dVec[i]) * sD->initSpeed[i]) / 2;
            solveEQSys();
            for (unsigned int i = 0; i < sD->dc; i++)
            {
                sD->bigStep_sorted[i].x = sD->disl_sorted[i].x - sD->x[i]; // bigStep_sorted = config[2]
                sD->bigStep_sorted[i].y = sD->disl_sorted[i].y;
            }

            calculateSpeedsAtTime(sD->bigStep_sorted, sD->speed, t_1);  //## kiszámolja a sebességeket newDisloc-ból, newDisloc = config[2]

            for (unsigned int i = 0; i < sD->dc; i++)
                sD->g[i] = sD->bigStep_sorted[i].x - sD->disl_sorted[i].x - stepSize * ((1 + sD->dVec[i]) * sD->speed[i] + (1 - sD->dVec[i]) * sD->initSpeed[i]) / 2;
            solveEQSys();
            for (unsigned int i = 0; i < sD->dc; i++)
                sD->bigStep_sorted[i].x -= sD->x[i];  // newDisloc = config[3]

            umfpack_di_free_numeric(&sD->Numeric);
        }

        //////////////////////////////////////////////////////////////////
        /// step stageII: first small step
        //////////////////////////////////////////////////////////////////
        stepSize /= 2;
        double t_1p2 = sD->simTime + stepSize; // the time point at the first small step
        {
            calcJacobianFromPrev();  //## kiszámolja a jakobit oldDisloc-ból; oldDisloc = config[1] !!
            for (unsigned int i = 0; i < sD->dc; i++)
                sD->g[i] = -stepSize * ((1 + sD->dVec[i]) * sD->speed2[i] + (1 - sD->dVec[i]) * sD->initSpeed[i]) / 2;
            solveEQSys();
            for (unsigned int i = 0; i < sD->dc; i++)
            {
                sD->firstSmall_sorted[i].x = sD->disl_sorted[i].x - sD->x[i]; // bigStep_sorted = config[2]
                sD->firstSmall_sorted[i].y = sD->disl_sorted[i].y;
            }

            calculateSpeedsAtTime(sD->firstSmall_sorted, sD->speed, t_1p2); //## kiszámolja a sebességeket newDisloc-ból, newDisloc = config[4]

            for (unsigned int i = 0; i < sD->dc; i++)
                sD->g[i] = sD->firstSmall_sorted[i].x - sD->disl_sorted[i].x - stepSize * ((1 + sD->dVec[i]) * sD->speed[i] + (1 - sD->dVec[i]) * sD->initSpeed[i]) / 2;
            solveEQSys();
            for (unsigned int i = 0; i < sD->dc; i++)
                sD->firstSmall_sorted[i].x -= sD->x[i]; // newDisloc = config[5]

            umfpack_di_free_numeric(&sD->Numeric);
        } // firstSmall_sorted = config[5]

        //////////////////////////////////////////////////////////////////
        /// step stage III: second small step
        //////////////////////////////////////////////////////////////////
        {
            calcJacobian(stepSize, sD->firstSmall_sorted); //## kiszámolja a jakobit oldDisloc-ból; oldDisloc = config[5]
            calculateSparseFormForJacobian();

            calculateSpeedsAtTimes(sD->firstSmall_sorted, sD->initSpeed2, sD->speed2, t_1p2, t_1); // ## kiszámolja a sebességeket oldDisloc-ból, oldDisloc = config[5] !!

            for (unsigned int i = 0; i < sD->dc; i++)
                sD->g[i] = -stepSize * ((1 + sD->dVec[i]) * sD->speed2[i] + (1 - sD->dVec[i]) * sD->initSpeed2[i]) / 2;

            solveEQSys();
            for (unsigned int i = 0; i < sD->dc; i++)
            {
                sD->secondSmall_sorted[i].x = sD->firstSmall_sorted[i].x - sD->x[i]; // newDisloc = config[6]
                sD->secondSmall_sorted[i].y = sD->firstSmall_sorted[i].y;
            }

            calculateSpeedsAtTime(sD->secondSmall_sorted, sD->speed2, t_1); // ## kiszámolja a sebességeket oldDisloc-ból, oldDisloc = config[6]

            for (unsigned int i = 0; i < sD->dc; i++)
                sD->g[i] = sD->secondSmall_sorted[i].x - sD->firstSmall_sorted[i].x - stepSize * ((1 + sD->dVec[i]) * sD->speed2[i] + (1 - sD->dVec[i]) * sD->initSpeed2[i]) / 2;

            solveEQSys();
            for (unsigned int i = 0; i < sD->dc; i++)
                sD->secondSmall_sorted[i].x -= sD->x[i];

            umfpack_di_free_numeric(&sD->Numeric);
        }

        double vsquare2 = absvalsq(sD->initSpeed2);
        double energyThisStep = (absvalsq(sD->initSpeed) + vsquare2) * sD->stepSize / 4;

        calculateXError();

        /// Precision related error handling
        if (pH->getMaxErrorRatioSqr() < 1)
        {
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

            calculateSpeedsAtTime(sD->disl_sorted, sD->initSpeed, sD->simTime);

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

            sD->standardOutputLog
                << sD->externalStressProtocol->extStress(sD->simTime) << "\t"
                << get_wall_time() - lastLogTime << "\t";

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

            lastLogTime = get_wall_time();
        }
        else
            sD->failedSteps++;

        sD->stepSize = pH->getNewStepSize(sD->stepSize);
        pH->reset();

        if (sD->isMaxStepSizeLimit && sD->maxStepSizeLimit < sD->stepSize)
            sD->stepSize = sD->maxStepSizeLimit;
    }

    sD->writeDislocationDataToFile(sD->endDislocationConfigurationPath);
    std::cout << "Simulation is done (" << get_wall_time() - startTime << " s).\n";
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

            double force = sD->tau.xy(dx, dy) * sD->b(i) * sD->b(j); // The sign of Burgers vectors can be deduced from i and j

            forces_A[i] += force;
            forces_A[j] -= force;
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

/**
    @brief calcJacobian:    calculates the Jacobian matrix containing the field derivatives multiplied with stepsize; modifies Ai, Ax, Ap, indexes, dVec
    @param stepsize:        how large time step should be made
    @param dislocs:         the actual positions of the dislocations
    @return int:            totalElementCounter, the total number of nonezero elements in the matrix J_{i,j}^k
*/
int Simulation::calcJacobian(double stepsize, const std::vector<DislwoB>& dislocs)
{
    int totalElementCounter = 0;

    for (unsigned int j = 0; j < sD->dc; j++)
    {
        // Previously calculated part
        for (unsigned int i = 0; i < j; i++)
        {
            double v = getElement(j, i);
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
            dx = dislocs[i].x - dislocs[j].x;
            normalize(dx);

            dy = dislocs[i].y - dislocs[j].y; // az y értékek konstansok, így ha a diszlokációk rendezve vannak tárolva, akkor csak az ellentétes típusúak között kell max 1 kivonást elvégezni a normalize függvényben
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

        if (subSum > 0)
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

    return totalElementCounter;
}

// like calcJacobianAndSpeedsAtTime but calculates the forces at two given time point where the forces can be different due to the load protocol
int Simulation::calcJacobianAndSpeedsAtTimes(double stepsize, const std::vector<DislwoB>& dislocs, std::vector<double>& forces_A, std::vector<double>& forces_B, double simTime_A, double simTime_B)
{
    int totalElementCounter = 0;
    std::fill(forces_A.begin(), forces_A.end(), 0);
    double extStress_A = sD->externalStressProtocol->extStress(simTime_A);
    double extStress_B = sD->externalStressProtocol->extStress(simTime_B);

    if (&forces_A != &forces_B)
        std::fill(forces_B.begin(), forces_B.end(), 0);

    for (unsigned int i = 0; i < sD->dc; i++)
    {
        // Previously calculated part
        for (unsigned int j = 0; j < i; j++)
        {
            double v = getElement(i, j);
            if (v != 0)
            {
                sD->Ai[totalElementCounter] = j;
                sD->Ax[totalElementCounter] = v;
                totalElementCounter++;
            }
        }
        // Add the diagonal element (it will be calculated later and the point defects now)
        sD->Ai[totalElementCounter] = i;

        double tmp = 0;
#ifdef USE_POINT_DEFECTS

        for (size_t l = 0; l < sD->pc; l++)
        {
            double dx;
            double dy;

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

        sD->Ax[totalElementCounter] = -tmp * stepsize;
        totalElementCounter++;
        // Totally new part
        for (unsigned int j = i + 1; j < sD->dc; j++)
        {
            double dx = dislocs[i].x - dislocs[j].x;
            normalize(dx);

            double dy = dislocs[i].y - dislocs[j].y;
            normalize(dy);

            double distSq = dx * dx + dy * dy;
            pH->updateTolerance(distSq * 0.0025, i);
            pH->updateTolerance(distSq * 0.0025, j);

            double exp2pix = exp(2 * M_PI * dx);
            double cos2piy = cos(2 * M_PI * dy);

            if (std::fabs(sqrt(distSq) - sD->cutOff) < 6 * sD->cutOff) // itt a 6-os szorzót simán ki kéne hagyni, meg a következő 3 sort is úgy, ahogy van, TODO
            {
                double multiplier = 1;
                if (distSq > sD->cutOffSqr)
                    multiplier = exp(-pow(sqrt(distSq) - sD->cutOff, 2) * sD->onePerCutOffSqr);

                sD->Ai[totalElementCounter] = j;
                sD->Ax[totalElementCounter] = stepsize * sD->b(j) * sD->b(i) * sD->tau.xy_diff_x(dx, dy, exp2pix, cos2piy) * multiplier;
                totalElementCounter++;
            }

            double force = sD->tau.xy(dx, dy, exp2pix, cos2piy) * sD->b(i) * sD->b(j); // The sign of Burgers vectors can be deduced from i and j
            forces_A[i] += force;
            forces_A[j] -= force;
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
        sD->Ap[i + 1] = totalElementCounter;
        if (&forces_A == &forces_B)
            forces_A[i] += extStress_A * sD->b(i);
        else
        {
            forces_B[i] = forces_A[i] + extStress_B * sD->b(i);
            forces_A[i] += extStress_A * sD->b(i);
        }
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

        if (subSum > 0)
            sD->dVec[j] = std::pow(1 - 1 / (subSum + 1), 2);
        else
            sD->dVec[j] = 0;
    }
       
    for (unsigned int j = 0; j < sD->dc; j++)
    {
        for (int i = sD->Ap[j]; i < sD->Ap[j + 1]; i++)
            sD->Ax[i] *= (1 + sD->dVec[sD->Ai[i]]) / 2;

        sD->Ax[sD->indexes[j]] += 1;
    }

    calculateSparseFormForJacobian();

    return totalElementCounter;
}

// like calcJacobian but also calculates the forces at the given time
int Simulation::calcJacobianAndSpeedsAtTime(double stepsize, const std::vector<DislwoB>& dislocs, std::vector<double>& forces, double simTime)
{
    return calcJacobianAndSpeedsAtTimes(stepsize, dislocs, forces, forces, simTime, simTime);
}

// Calculates the new Jacobian J_{i,j}^k from the previous one by halfing the non-diagonal elements and also the weights
void Simulation::calcJacobianFromPrev()
{
    std::vector<double> new_weig(sD->dc);                                // the new weight values
    std::vector<double> weig_fac(sD->dc);                                // for the new J_{i,j} elements

    // calculates the new_w values from the old dVec values
    std::transform(sD->dVec.begin(), sD->dVec.end(), new_weig.begin(), [](double old_weight) {return old_weight == 0 ? 0 : old_weight / pow(2 - sqrt(old_weight), 2); });

    // calcualtes the weight factors needed to calculate J_{i,j}^k
    std::transform(
        sD->dVec.begin(), sD->dVec.end(),
        new_weig.begin(),
        weig_fac.begin(),
        [](double old_val, double new_val) {return (1 + new_val) / (1 + old_val); }
    );

    for (unsigned int j = 0; j < sD->dc; j++)                   // for all d-d interaction
    {
        for (int i = sD->Ap[j]; i < sD->Ap[j + 1]; i++)         // go through the row indicies Ai for column j
        {
            double factor = weig_fac[sD->Ai[i]];
            if (sD->Ai[i] == j)                                 // diagonal element
                sD->Ax[i] = 1 + (sD->Ax[i] - 1) / 2 * factor;
            else                                                // offdiagonal element
                sD->Ax[i] = sD->Ax[i] / 2 * factor;
        }
    }
    sD->dVec = std::move(new_weig);

    calculateSparseFormForJacobian();
}

// calculates the Burgers' vector weighted sum of the displacements
double Simulation::calculateStrainIncrement(const std::vector<DislwoB>& old, const std::vector<DislwoB>& newD) const
{
    double ret = 0;
    for (size_t i = 0; i < old.size(); i++)
        ret += sD->b(i) * (newD[i].x - old[i].x); // x is not in the range of [-0.5: 0.5), the difference is the real distance

    return ret;
}

// calculates the forces (therefore, the speed too) between all d-d and d-p (d: dislocation, p: fixed point defect) if dislocations are at dis and the force will be calculated from simTime
void Simulation::calculateSpeedsAtTimes(const std::vector<DislwoB>& dis, std::vector<double>& forces_A, std::vector<double>& forces_B, double simTime_A, double simTime_B) const
{
    double extStress_A = sD->externalStressProtocol->extStress(simTime_A);
    double extStress_B = sD->externalStressProtocol->extStress(simTime_B);
    calculateSpeedsAtStresses(dis, forces_A, forces_B, extStress_A, extStress_B);
}

// calculates the forces (therefore, the speed too) between all d-d and d-p (d: dislocation, p: fixed point defect) if dislocations are at dis and the force will be calculated from simTime
void Simulation::calculateSpeedsAtTime(const std::vector<DislwoB>& dis, std::vector<double>& forces, double simTime) const
{
    calculateSpeedsAtTimes(dis, forces, forces, simTime, simTime);
}

// calculates the difference of the large and two small steps and store it for all dislocs and saves the (largest relative to toleranceAndError[ID].first)
void Simulation::calculateXError()
{
    for (unsigned int i = 0; i < sD->dc; i++)
    {
        double tmp = fabs(sD->bigStep_sorted[i].x - sD->secondSmall_sorted[i].x);
        pH->updateError(tmp, i);
    }
}

// constructor: saves the sD pointer, creates pH object, sets output format
Simulation::Simulation(std::shared_ptr<SimulationData> sD) :
    sD(sD),
    pH(new PrecisionHandler)
{
    // Format setting
    sD->standardOutputLog << std::scientific << std::setprecision(16);

    pH->setMinPrecisity(sD->prec);
    pH->setSize(sD->dc);
}

// returns the i,j element of the dense matrix A
double Simulation::getElement(int i, int j) const
{
    return getElement(i, sD->Ap[j], sD->Ap[j + 1]);
}

double Simulation::getElement(int j, int si, int ei) const
{
    int len = ei - si;
    if (len > 1)
    {
        int tmp = len / 2;

        if (sD->Ai[si + tmp] > j)
            return getElement(j, si, si + tmp); // this is not inefficient https://softwareengineering.stackexchange.com/a/182334/350074
        else
            return getElement(j, si + tmp, ei);

    }
    else if (sD->Ai[si] == j)
        return sD->Ax[si];

    return 0;
}

// modifies only Symbolic, constant for Ap, Ai, Ax; doesn't read dVec, indexes, speeds or positions
void Simulation::calculateSparseFormForJacobian()
{
    (void)umfpack_di_symbolic(sD->dc, sD->dc, sD->Ap, sD->Ai, sD->Ax, &(sD->Symbolic), sD->null, sD->null);
    (void)umfpack_di_numeric(sD->Ap, sD->Ai, sD->Ax, sD->Symbolic, &(sD->Numeric), sD->null, sD->null);
    umfpack_di_free_symbolic(&(sD->Symbolic));
}

// solves A * Δx = g for Δx
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
