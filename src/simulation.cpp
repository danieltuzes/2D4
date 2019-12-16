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

using namespace sdddstCore;

// constructor: saves the sD pointer, creates pH object, sets output format
Simulation::Simulation(std::shared_ptr<SimulationData> sD) :
    sD(sD),
    pH(new PrecisionHandler)
{
    // Format setting
    sD->standardOutputLog << std::scientific << std::setprecision(16);

    pH->setMinPrecisity(sD->prec);
    pH->setDipolePrecisity(sD->dipole_prec);
    pH->setSize(sD->dc);
}

// starts the simulation, all required data are present; modifies the state of the system and creates files
void Simulation::run()
{
#pragma region stage 0: before loop, first run

    // cout date and time
    {
        std::time_t start_date_and_time = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        std::cout << "Simulation started on " << std::ctime(&start_date_and_time) << std::endl;
    }

    double startTime = get_wall_time();     // start time in seconds
    double lastLogTime = get_wall_time();   // the time of the last write into the logfile
    double energy = 0;                      // 

    // First two lines in the log file
    {
        calculateSpeedsAtTime(sD->disl_sorted, sD->initSpeed, sD->simTime); // calculated only for the purpose to be able to write out initial speed values

        double sumAvgSp = std::accumulate(sD->initSpeed.begin(), sD->initSpeed.end(), 0., [](double a, double b) {return a + fabs(b); }) / sD->dc;
        double vsquare = std::accumulate(sD->initSpeed.begin(), sD->initSpeed.end(), 0., [](double a, double b) {return a + b * b; });
        double extStress = sD->externalStressProtocol->extStress(sD->simTime);

        sD->standardOutputLog <<
            "# simTime(" << 1 << ")\t"
            "successfullSteps(" << 2 << ")\t"
            "failedSteps(" << 3 << ")\t"
            "maxErrorRatioSqr(" << 4 << ")\t"
            "maxErrorRatioID(" << 5 << ")\t"
            "sumAvgSp(" << 6 << ")\t"
            "cutOff(" << 7 << ")\t"
            "orderParameter(" << 8 << ")\t"
            "externalStress(" << 9 << ")\t"
            "timeStageI - III(" << 10 << ")\t"
            "strain(" << 11 << ")\t"
            "vSquare(" << 12 << ")\t"
            "energy(" << 13 << ")\t"
            "wallTimeElapsed(" << 14 << ")\t"
            "stepSize(" << 15 << ")\t" << std::endl;
        sD->standardOutputLog <<
            sD->simTime << "\t" <<
            sD->succesfulSteps << "\t" <<
            sD->failedSteps << "\t" <<
            0 << "\t" <<
            0 << "\t" <<
            sumAvgSp << "\t" <<
            sD->cutOff << "\t" <<
            "-" << "\t" <<
            extStress << "\t" <<
            "-" << "\t" <<
            sD->totalAccumulatedStrainIncrease << "\t" <<
            vsquare << "\t" <<
            energy << "\t" <<
            0 << "\t" <<
            sD->stepSize << std::endl;
    }
#pragma endregion

    while (!interrupt())
    {
#pragma region step stage I: one large step
        double t_0__ = sD->simTime;                     // simTime at the beginning of a large step
        double t_1__ = sD->simTime + sD->stepSize;      // simTime at the end of a large step
        double t_1p2 = sD->simTime + sD->stepSize / 2;  // the time point at the first small step
        // int nz;                                         // the number of non0 elements in the Jacobian
        {
            //## calculates the Jacobian from disl_sorted with stepSize (therefore at stage II Jacobian can be deduced for stepSize/2); disl_sorted = config[1]
            // speed is also calculated from config[1] at time t_1 (therefore at stage II speed can be deduced at time t_1p2), and in case of failure step, for the next steps stage II can be reused again !!
            calcJacobianAndSpeedsAtTimes(sD->stepSize, sD->disl_sorted, sD->initSpeed, sD->speed, t_0__, t_1__);

            for (unsigned int i = 0; i < sD->dc; i++)
                sD->g[i] = -sD->stepSize * ((1 + sD->dVec[i]) * sD->speed[i] + (1 - sD->dVec[i]) * sD->initSpeed[i]) / 2;  // initSpeed has been previously already calculated
            solveEQSys("stage I, 1.");

            for (unsigned int i = 0; i < sD->dc; i++)
            {
                sD->bigStep_sorted[i].x = sD->disl_sorted[i].x - sD->x[i]; // bigStep_sorted = config[2]
                sD->bigStep_sorted[i].y = sD->disl_sorted[i].y;
            }

            calculateSpeedsAtTime(sD->bigStep_sorted, sD->speed2, t_1__);  //## calculates speeds from sD->bigStep_sorted = config[2]

            calcGSolveAndUpdate(sD->bigStep_sorted, sD->disl_sorted, sD->stepSize, sD->speed2, sD->initSpeed, "stage I, 2.");
        }
#pragma endregion

#pragma region step stage II: first small step
        {
            calcJacobianAndSpeedsFromPrev(t_0__, t_1p2, t_1__);  //## calculates Jacobian and speed from disl_sorted = config[1], speed is at time t_1p2 = sD->simTime + sD->stepSize / 2 !!

            for (unsigned int i = 0; i < sD->dc; i++)
                sD->g[i] = -sD->stepSize / 2 * ((1 + sD->dVec[i]) * sD->speed[i] + (1 - sD->dVec[i]) * sD->initSpeed[i]) / 2;
            solveEQSys("stage II, 3.");
            for (unsigned int i = 0; i < sD->dc; i++)
            {
                sD->firstSmall_sorted[i].x = sD->disl_sorted[i].x - sD->x[i]; // firstSmall_sorted = config[4]
                sD->firstSmall_sorted[i].y = sD->disl_sorted[i].y;
            }

            calculateSpeedsAtTime(sD->firstSmall_sorted, sD->speed2, t_1p2); //## calculates speed from firstSmall_sorted = config[4]

            calcGSolveAndUpdate(sD->firstSmall_sorted, sD->disl_sorted, sD->stepSize / 2, sD->speed2, sD->initSpeed, "stage II, 4.");
        }
#pragma endregion

#pragma region step stage III: second small step
        {
            calcJacobianAndSpeedsAtTimes(sD->stepSize / 2, sD->firstSmall_sorted, sD->initSpeed2, sD->speed2, t_1p2, t_1__); //## calculates the Jacobian from firstSmall_sorted = config[5]

            for (unsigned int i = 0; i < sD->dc; i++)
                sD->g[i] = -sD->stepSize / 2 * ((1 + sD->dVec[i]) * sD->speed2[i] + (1 - sD->dVec[i]) * sD->initSpeed2[i]) / 2;

            solveEQSys("stage III, 5.");

            for (unsigned int i = 0; i < sD->dc; i++)
            {
                sD->secondSmall_sorted[i].x = sD->firstSmall_sorted[i].x - sD->x[i]; // secondSmall_sorted = config[6]
                sD->secondSmall_sorted[i].y = sD->firstSmall_sorted[i].y;
            }

            calculateSpeedsAtTime(sD->secondSmall_sorted, sD->speed2, t_1__); // ## calculates speeds from secondSmall_sorted = config[6]

            calcGSolveAndUpdate(sD->secondSmall_sorted, sD->firstSmall_sorted, sD->stepSize / 2, sD->speed2, sD->initSpeed2, "stage III, 6.");
        }
#pragma endregion

#pragma region step stage IV: accept or retry disl_sorted 
        calculateXError();

        if (pH->getMaxErrorRatioSqr() < 1) // accept step
        {
            sD->simTime += sD->stepSize;
            sD->succesfulSteps++;

            // for simulation criterium and output
            if (sD->calculateStrainDuringSimulation)
            {
                sD->totalAccumulatedStrainIncrease += calculateStrainIncrement(sD->disl_sorted, sD->firstSmall_sorted);
                sD->totalAccumulatedStrainIncrease += calculateStrainIncrement(sD->firstSmall_sorted, sD->secondSmall_sorted);
            }

            sD->disl_sorted.swap(sD->secondSmall_sorted);
            for (auto& disl : sD->disl_sorted)
                normalize(disl.x);

            double sumAvgSp = std::accumulate(sD->initSpeed2.begin(), sD->initSpeed2.end(), 0., [](double a, double b) {return a + fabs(b); }) / sD->dc; // for logfile and cutoff modification

            // logfile write
            {
                double vsquare_end = absvalsq(sD->speed2);  // best approximation for speeds at the end of the 2nd small step
                energy += (absvalsq(sD->initSpeed) + 2 * absvalsq(sD->initSpeed2) + vsquare_end) * sD->stepSize / 4;

                sD->standardOutputLog << sD->simTime << "\t"
                    << sD->succesfulSteps << "\t"
                    << sD->failedSteps << "\t"
                    << pH->getMaxErrorRatioSqr() << "\t"
                    << pH->maxErrorRatioID() << "\t"
                    << sumAvgSp << "\t"
                    << sD->cutOff << "\t";

                if (sD->orderParameterCalculationIsOn)
                    sD->standardOutputLog << calculateOrderParameter(sD->speed) << "\t";
                else
                    sD->standardOutputLog << "-" << "\t";

                sD->standardOutputLog
                    << sD->externalStressProtocol->extStress(sD->simTime) << "\t"
                    << get_wall_time() - lastLogTime << "\t";

                if (sD->calculateStrainDuringSimulation)
                    sD->standardOutputLog << sD->totalAccumulatedStrainIncrease << "\t";
                else
                    sD->standardOutputLog << "-" << "\t";

                sD->standardOutputLog << vsquare_end << "\t"
                    << energy << "\t"
                    << get_wall_time() - startTime << "\t"
                    << sD->stepSize << std::endl;

                lastLogTime = get_wall_time();
            }

            if (sD->isSpeedThresholdForCutoffChange && sD->speedThresholdForCutoffChange > sumAvgSp) // must be after printout to logfile
            {
                sD->cutOffMultiplier = 1e20;
                sD->updateCutOff();
            }

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

            if (sD->isSaveSubConfigs)
            {
                if ((!sD->inAvalanche && sD->subConfigDelay <= sD->subconfigDistanceCounter) ||
                    (sD->inAvalanche && sD->subConfigDelayDuringAvalanche <= sD->subconfigDistanceCounter))
                {
                    sD->subconfigDistanceCounter = 0;
                    std::stringstream ss;
                    ss << std::setprecision(9);
                    ss << sD->simTime;
                    sD->writeDislocationDataToFile(sD->subConfigPath + ss.str() + ".dconf");
                }
                else
                    sD->subconfigDistanceCounter++;
            }

        }
        else
            sD->failedSteps++;

        sD->stepSize = pH->getNewStepSize(sD->stepSize);

        // if stepSize was limited last time, the proposed stepSize has been saved to stepSizeBeforeWriteout, so there is no need to start with a very small stepSize
        if (sD->stepSizeBeforeWriteout != 0)
        {
            if (pH->getMaxErrorRatioSqr() < 1)  // use the increased step size in case it was a successful step
                sD->stepSize = std::max(sD->stepSize, sD->stepSizeBeforeWriteout);

            sD->stepSizeBeforeWriteout = 0;     // otherwise, forget about the suggestion
        }

        if (sD->isSaveSubConfigs && sD->subConfigTimes != 0)
        {
            double remainder = 0; // the upcoming step size
            if (sD->subConfigTimesType == 'a')
                remainder = sD->subConfigTimes - fmod(sD->simTime, sD->subConfigTimes);
            else if (sD->subConfigTimesType == 'b')
            {
                double exponent = nextafter(log(sD->simTime / sD->initStepSize) / log(sD->subConfigTimes), INFINITY);
                int nextExp = std::max(int(exponent) + 1, 1);
                double nextTime = sD->initStepSize * pow(sD->subConfigTimes, nextExp);
                remainder = nextTime - sD->simTime;
            }

            if (remainder < sD->stepSize)
            {
                sD->stepSizeBeforeWriteout = sD->stepSize;      // after the write out, the next time step will be at least stepSizeBeforeWriteout, because stepSize can be decreased a lot in the next line
                sD->stepSize = nextafter(remainder, INFINITY);  // the new stepSize can be unnecessarily small for the precision, but is a requirement for further analysis of the output 
                sD->subconfigDistanceCounter = sD->subConfigDelay + sD->subConfigDelayDuringAvalanche; // writeDislocationDataToFile will be triggered next time
            }
        }

        pH->reset();

        if (sD->isMaxStepSizeLimit && sD->maxStepSizeLimit < sD->stepSize)
            sD->stepSize = sD->maxStepSizeLimit;

#pragma endregion
    }

    sD->writeDislocationDataToFile(sD->endDislocationConfigurationPath);
    std::cout << "Simulation is done (" << get_wall_time() - startTime << " s and " << sD->succesfulSteps + sD->failedSteps << " steps).\n";
}

#pragma region Functions directly called from run

// calculates the forces (therefore, the speeds too) between all d-d and d-p (d: dislocation, p: fixed point defect)
// dislocations are at dis and the force will be calculated from simTime
void Simulation::calculateSpeedsAtTime(const std::vector<DislwoB>& dis, std::vector<double>& forces, double simTime) const
{
    std::fill(forces.begin(), forces.end(), 0);
    double extStress = sD->externalStressProtocol->extStress(simTime);

    for (unsigned int i = 0; i < sD->dc; i++) // typically unefficient
    {
        for (unsigned int j = i + 1; j < sD->dc; j++)
        {
            double dx = dis[i].x - dis[j].x;
            normalize(dx);

            double dy = dis[i].y - dis[j].y;
            normalize(dy);

            double distSq = dx * dx + dy * dy; // proportional to distance square
            pH->updateTolerance(distSq, i);
            pH->updateTolerance(distSq, j);

            double force = sD->tau.xy(dx, dy) * sD->b(i) * sD->b(j); // The sign of Burgers vectors can be deduced from i and j

            forces[i] += force;
            forces[j] -= force;
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
        forces[i] += extStress * sD->b(i);
    }
}

/**
@brief calcJacobianAndSpeedsAtTimes:    calculates the Jacobian matrix containing the field derivatives multiplied with stepsize; modifies Ai, Ax, Ap, indexes, dVec; also calculates the forces at two different time
@param stepsize:                        how large time step should be made
@param dislocs:                         the actual positions of the dislocations
@param forces_A:                        the estimated speeds of the particles at simTime_A
@param forces_B:                        the other estimated speeds of the particles at simTime_B
@param simTime_A:                       the which time should the external stress be evaluated to estimate the forces
@param simTime_B:                       the which other time should the external stress be evaluated to estimate the forces
@return int:                            totalElementCounter, the total number of nonezero elements in the matrix J_{i,j}^k
*/
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
        if (sD->currentStorageSize - totalElementCounter < sD->dc) // if free memory is less than sD->dc
            sD->increaseCurrentStorageSize(totalElementCounter);

        double subSum = 0;
        // Previously calculated part
        for (unsigned int j = 0; j < i; j++)
        {
            double v = getElement(i, j);
            if (v != 0)
            {
                sD->Ai[totalElementCounter] = j;
                sD->Ax[totalElementCounter] = v;

                subSum -= v;

                totalElementCounter++;
            }
        }
        // Add the diagonal element (it will be calculated later and the point defects now)
        sD->Ai[totalElementCounter] = i;
        sD->indexes[i] = totalElementCounter;

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
            pH->updateTolerance(distSq, i);
            pH->updateTolerance(distSq, j);

            double exp2pix = exp2(PI2LN2INV * dx);
            double cos2piy = cos(2 * M_PI * dy);

            double force = sD->tau.xy(dx, dy, exp2pix, cos2piy) * sD->b(i) * sD->b(j); // The sign of Burgers vectors can be deduced from i and j

            double cutOffRangeSq = 49; // the weight in the Jacobian decreases after the distance reaches cutOff, and at cutOffRange * cutOff, wight is 0
            if (sD->heavisideCutoff)
                cutOffRangeSq = 1;

            if (distSq < cutOffRangeSq * sD->cutOffSqr)
            {
                double v = stepsize * sD->b(j) * sD->b(i) * sD->tau.xy_diff_x(dx, dy, exp2pix, cos2piy);

                if (!sD->heavisideCutoff && distSq > sD->cutOffSqr)
                    v *= exp(-pow(sqrt(distSq) - sD->cutOff, 2) * sD->onePerCutOffSqr);

                sD->Ai[totalElementCounter] = j;
                sD->Ax[totalElementCounter] = v;

                subSum -= v;

                totalElementCounter++;
            }

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

        sD->Ax[sD->indexes[i]] = subSum;
        sD->dVec[i] = weight(subSum, sD->weightFunc);
    }

    for (unsigned int i = 0; i < sD->dc; i++)
    {
        for (int j = sD->Ap[i]; j < sD->Ap[i + 1]; j++)
            sD->Ax[j] *= (1 + sD->dVec[sD->Ai[j]]) / 2;

        sD->Ax[sD->indexes[i]] += 1;
    }

    calculateSparseFormForJacobian();

    return totalElementCounter;
}

/**
@brief calcJacobian:    like calcJacobianAndSpeedsAtTimes: calculates the Jacobian matrix containing the field derivatives multiplied with stepsize; modifies Ai, Ax, Ap, indexes, dVec; also calculates the force but at only 1 time point
@param stepsize:        how large time step should be made
@param dislocs:         the actual positions of the dislocations
@param forces:          the speeds of the particles
@param simTime:         the value of the external force will be calcualted from tim
@return int:            totalElementCounter, the total number of nonezero elements in the matrix J_{i,j}^k
*/
int Simulation::calcJacobianAndSpeedsAtTime(double stepsize, const std::vector<DislwoB>& dislocs, std::vector<double>& forces, double simTime)
{
    return calcJacobianAndSpeedsAtTimes(stepsize, dislocs, forces, forces, simTime, simTime);
}

// Calculates the new Jacobian J_{i,j}^k from the one in the memory by halfing the non-diagonal elements and also the weights; recalculates speed too bc of the different external stress values
// baseTime: the the time the Jacobian was in the correct state
// newTime:  the new time point, where the new configuration should be predicted
// oldTime:  the old time point, where the last prediction was made with the Jacobian
void Simulation::calcJacobianAndSpeedsFromPrev(double baseTime, double newTime, double oldTime)
{
    //////////////////////////////////////////////////////////////////
    // Calculating the Jacobian
    //////////////////////////////////////////////////////////////////

    std::vector<double> new_weig(sD->dc);                               // the new weight values
    std::vector<double> weig_fac(sD->dc);                               // the multiplying factor for the new J_{i,j} elements

    // calculates the new_w values from the old dVec values

    double factor = (newTime - baseTime) / (oldTime - baseTime);        // usually, it is a factor of 1/2, if 1 large step = 2 small steps
    char T = sD->weightFunc;
    std::transform(
        sD->dVec.begin(), sD->dVec.end(),
        new_weig.begin(),
        [factor, T](double old_weight)
        {
            if (old_weight == 0)
                return 0.;
            else
            {
                double prevSum = weightInv(old_weight, T);
                return weight(prevSum * factor, T);
            }
        });

    // calcualtes the weight factors needed to calculate J_{i,j}^k
    std::transform(
        sD->dVec.begin(), sD->dVec.end(),
        new_weig.begin(),
        weig_fac.begin(),
        [](double old_val, double new_val)
        {
            return (1 + new_val) / (1 + old_val);
        });

    for (unsigned int j = 0; j < sD->dc; j++)                   // for all d-d interaction
    {
        for (int i = sD->Ap[j]; i < sD->Ap[j + 1]; i++)         // go through the row indicies Ai for column j
        {
            double factor = weig_fac[sD->Ai[i]];
            if (sD->Ai[i] == int(j))                            // diagonal element
                sD->Ax[i] = 1 + (sD->Ax[i] - 1) / 2 * factor;
            else                                                // offdiagonal element
                sD->Ax[i] = sD->Ax[i] / 2 * factor;
        }
    }
    sD->dVec = std::move(new_weig);

    calculateSparseFormForJacobian();

    //////////////////////////////////////////////////////////////////
    // Calculating the new speeds
    //////////////////////////////////////////////////////////////////
    double extStress_old = sD->externalStressProtocol->extStress(oldTime);
    double extStress_new = sD->externalStressProtocol->extStress(newTime);
    double extStress_dif = extStress_new - extStress_old;

    for (unsigned int i = 0; i < sD->dc; ++i)
        sD->speed[i] += extStress_dif * sD->b(i);
}

/**
@brief calcGSolveAndUpdate: calculates the new g vector, solves the linear equations and updates the dislocation positions
@param new_disloc:          the container of the target dislocation arrangement
@param old_config:          the actual dislocation configuration
@param stepSize:            the time size of the step
@param endSpeed:            the estimated speeds at the end of the step
@param initSpeed:           the speeds at the beginning of the step
*/
void Simulation::calcGSolveAndUpdate(std::vector<DislwoB>& new_disloc, const std::vector<DislwoB>& old_config, double stepSize, const std::vector<double>& endSpeed, const std::vector<double>& initSpeed, std::string label)
{
    for (unsigned int i = 0; i < sD->dc; i++)
        sD->g[i] = new_disloc[i].x - old_config[i].x - stepSize * ((1 + sD->dVec[i]) * endSpeed[i] + (1 - sD->dVec[i]) * initSpeed[i]) / 2;
    solveEQSys(label);
    for (unsigned int i = 0; i < sD->dc; i++)
        new_disloc[i].x -= sD->x[i];

    umfpack_di_free_numeric(&sD->Numeric);
}

#pragma endregion

#pragma region Functions related to the sparse matrix handling and solving


// modifies only Symbolic, constant for Ap, Ai, Ax; doesn't read dVec, indexes, speeds or positions
void Simulation::calculateSparseFormForJacobian()
{
    (void)umfpack_di_symbolic(sD->dc, sD->dc, sD->Ap, sD->Ai, sD->Ax, &(sD->Symbolic), sD->null, sD->null);
    (void)umfpack_di_numeric(sD->Ap, sD->Ai, sD->Ax, sD->Symbolic, &(sD->Numeric), sD->null, sD->null); // umfpack_di_free_numeric is in calcGSolveAndUpdate
    umfpack_di_free_symbolic(&(sD->Symbolic));
}

// solves A * Δx = g for Δx
void Simulation::solveEQSys(std::string label)
{
    int status = umfpack_di_solve(UMFPACK_A, sD->Ap, sD->Ai, sD->Ax, sD->x, sD->g.data(), sD->Numeric, sD->null, sD->null);
    if (status == UMFPACK_WARNING_singular_matrix)
    {
        std::cerr
            << "Warning: singular matrix at time " << sD->simTime
            << " using stepSize " << sD->stepSize
            << " after " << sD->succesfulSteps << " succesful and " << sD->failedSteps << " failed steps at position labeled as " << label << ". All nonfinite result values should be removed." << std::endl;
        for (unsigned int i = 0; i < sD->dc; ++i)
            if (!((-1 < sD->x[i]) && (sD->x[i] < 1)))
            {
                std::cerr << "x[" << i << "] = " << sD->x[i] << " is problematic, a new value 0 is assigned to it to avoid further errors\n";
                sD->x[i] = 0;
            }
    }
}

// return the jth line element that is in between si and ei indices
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

// returns the i,j element of the dense matrix A
double Simulation::getElement(int i, int j) const
{
    return getElement(i, sD->Ap[j], sD->Ap[j + 1]);
}

#pragma endregion

#pragma region Functions required for analysis and outputs

// calculates the difference of the large and two small steps and store it for all dislocs and saves the (largest relative to toleranceAndError[ID].first)
void Simulation::calculateXError()
{
    for (unsigned int i = 0; i < sD->dc; i++)
    {
        double tmp = fabs(sD->bigStep_sorted[i].x - sD->secondSmall_sorted[i].x);
        pH->updateMaxErrorToleranceRatioSq(tmp, i);
        //pH->updateError(tmp, i);
    }
}

// calcualtes the deformation speed: Burgers' vector signed sum of the speeds
double Simulation::calculateOrderParameter(const std::vector<double>& speeds) const
{
    double orderParameter = 0;
    for (size_t i = 0; i < sD->dc; i++)
        orderParameter += sD->b(i) * speeds[i];

    return orderParameter;
}

// calculates the Burgers' vector weighted sum of the displacements
double Simulation::calculateStrainIncrement(const std::vector<DislwoB>& old, const std::vector<DislwoB>& newD) const
{
    double ret = 0;
    for (size_t i = 0; i < old.size(); i++)
        ret += sD->b(i) * (newD[i].x - old[i].x); // x is not in the range of [-0.5: 0.5), the difference is the real distance

    return ret;
}

bool Simulation::interrupt() const
{
    if (sD->isTimeLimit && sD->simTime > sD->timeLimit)
    {
        std::cout << "time-limit (" << sD->timeLimit << ") has been reached, "
            "simTime = " << sD->simTime << "\n";
        return true;
    }
    if (sD->isStrainIncreaseLimit && sD->totalAccumulatedStrainIncrease > sD->totalAccumulatedStrainIncreaseLimit)
    {
        std::cout << "strain-increase-limit (" << sD->totalAccumulatedStrainIncreaseLimit << ") has been reached, "
            "strain = " << sD->totalAccumulatedStrainIncrease << "\n";
        return true;
    }
    if (sD->isStepCountLimit && sD->succesfulSteps > sD->stepCountLimit)
    {
        std::cout << "step-count-limit (" << sD->stepCountLimit << ") has been reached, "
            "succesfulSteps = " << sD->succesfulSteps << "\n";
        return true;
    }
    if (sD->countAvalanches && sD->avalancheCount > sD->avalancheTriggerLimit)
    {
        std::cout << "avalanche-detection-limit (" << sD->avalancheTriggerLimit << ") has been reached, "
            "avalancheCount = " << sD->avalancheCount << "\n";
        return true;
    }

    return false;
}

#pragma endregion

#ifdef DEBUG_VERSION

// prints out the selected container's x values to fname 
template <typename T>
void Simulation::printOut(std::string fname, const std::vector<T>& m_vector) const
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

// prints out size number of elements from the selected array's values to fname
template <typename T>
void Simulation::printOut(std::string fname, T* array, int size) const
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
void Simulation::printAll(std::string fname, unsigned int nz) const
{
    printOut("speed" + fname, sD->speed);
    printOut("speed2" + fname, sD->speed2);
    printOut("initSpeed" + fname, sD->initSpeed);
    printOut("initSpeed2" + fname, sD->initSpeed2);
    printOut("disl_sorted" + fname, sD->disl_sorted);
    printOut("firstSmall_sorted" + fname, sD->firstSmall_sorted);
    printOut("secondSmall_sorted" + fname, sD->secondSmall_sorted);
    printOut("bigStep_sorted" + fname, sD->bigStep_sorted);
    printOut("dVec" + fname, sD->dVec);
    printOut("indexes" + fname, sD->indexes);
    printOut("g" + fname, sD->g);
    printOut("Ax" + fname, sD->Ax, nz);
    printOut("Ai" + fname, sD->Ai, nz);
    printOut("Ap" + fname, sD->Ap, sD->dc + 1);
    printOut("x" + fname, sD->x, sD->dc);
    Ax_nonSparse_to_file(fname);

}

// checks if all values in the container are finite
bool Simulation::isFinite(std::vector<double> m_vector, double lb, double ub) const
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
bool Simulation::isFinite(std::vector<DislwoB> m_vector, double lb, double ub) const
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
bool Simulation::isFinite(double* m_array, size_t nz, double lb, double ub) const
{
    std::vector<double> m_vector(m_array, m_array + nz);
    return isFinite(m_vector, lb, ub);
}

// checks if all the containers and arrays up to nz number of elements contain only finite values and
// print out results to labeled filenames, label should match ^[\w,\s-]+
bool Simulation::isAllFinite(size_t nz, std::string label) const
{
    if (sD->simTime < 8.433e-04)
        return true;

    bool allfinite = true;
    if (!isFinite(sD->speed, -1E18, 1E18))
    {
        std::cerr << "speed can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(sD->speed2, -1E18, 1E18))
    {
        std::cerr << "speed2 can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(sD->initSpeed, -1E18, 1E18))
    {
        std::cerr << "initSpeed can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(sD->initSpeed2, -1E18, 1E18))
    {
        std::cerr << "initSpeed2 can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(sD->disl_sorted, -1, 1))
    {
        std::cerr << "disl_sorted can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(sD->firstSmall_sorted, -1, 1))
    {
        std::cerr << "firstSmall_sorted can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(sD->secondSmall_sorted, -1, 1))
    {
        std::cerr << "secondSmall_sorted can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(sD->bigStep_sorted, -1, 1))
    {
        std::cerr << "bigStep_sorted can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(sD->dVec, 0, 1))
    {
        std::cerr << "dVec can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(sD->g, -1E18, 1E18))
    {
        std::cerr << "g can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(sD->Ax, nz, -(int)(sD->dc), sD->dc))
    {
        std::cerr << "Ax can lead to nonfinite result\n";
        allfinite = false;
    }
    if (!isFinite(sD->x, sD->dc, -1, 1))
    {
        std::cerr << "x can lead to nonfinite result\n";
        allfinite = false;
    }

    if (!allfinite)
    {
        std::cerr
            << "The error was found at label "
            << label << "\t"
            << sD->succesfulSteps << "\t"
            << sD->failedSteps << "\t"
            << sD->stepSize << "\t"
            << sD->simTime << "\n"
            << "Debugging files will be creatred."
            << std::endl;

        std::stringstream ss;
        ss << label << "_" << sD->succesfulSteps << "_" << sD->failedSteps << ".txt";
        printAll(ss.str(), nz);
        std::cerr << std::endl;
    }
    return allfinite;
}

void Simulation::Ax_nonSparse_to_file(std::string fname) const
{
    std::ofstream of("Ax_nonSparse" + fname, std::ios_base::app);
    if (!of)
    {
        std::cerr << "Cannot create or append to file " << fname << std::endl;
        return;
    }

    std::vector<std::vector<double>> ret(sD->dc, std::vector<double>(sD->dc, 0));
    for (size_t i = 0; i < sD->dc; ++i)
    {
        for (size_t j = 0; j < sD->dc; ++j)
            of << getElement(i, j) << "\t";
        of << "\n";
    }
}


#endif