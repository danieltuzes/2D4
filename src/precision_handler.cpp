// 
// precision_handler.cpp : contains the function defintions for precision_handler.h, also included in simulation.h

#include "constants.h"
#include "precision_handler.h"

using namespace sdddstCore;

PrecisionHandler::PrecisionHandler() :
    minPrecisity(0),
    minPrecisitySq(0),
    maxErrorToleranceRatioSq(0),
    selectedID(0),
    use_dipoleCrit(0),
    dipolePrecisity(0),
    dipolePrecisitySq(0) {}

// if use_dipoleCrit, sets the tolerance container to its proper size
void PrecisionHandler::setSize(unsigned int size)
{
    if (use_dipoleCrit)
        tolerance.resize(size, minPrecisitySq);
}

// if use_dipoleCrit, resets the tolerance to minPrecisity, i.e. it forgets about the closest dislocations
void PrecisionHandler::reset()
{
    maxErrorToleranceRatioSq = 0;
    if (use_dipoleCrit)
        tolerance.assign(tolerance.size(), minPrecisitySq);
}

// if use_dipoleCrit is not 0 then if distanceSq is smaller than the prescribed minPrecisitySq, the function overwrites it for that particle with distanceSq
void PrecisionHandler::updateTolerance(double distSq, unsigned int ID)
{
    if (use_dipoleCrit && distSq < minPrecisitySq && distSq < tolerance[ID])
        tolerance[ID] = distSq;
}

// calculates if error compared to the required precision (absolute or dipolePrecisity, depending on use_dipoleCrit) and stores the renitent's ID to selectedID
void PrecisionHandler::updateMaxErrorToleranceRatioSq(double error, unsigned int ID)
{
    double tmp;
    if (use_dipoleCrit)
        tmp = error * error / tolerance[ID];
    else
        tmp = error * error / minPrecisitySq;

    if (tmp > maxErrorToleranceRatioSq)
    {
        maxErrorToleranceRatioSq = tmp;
        selectedID = ID;
    }
}

// heuristic guess for the new stepsize
double PrecisionHandler::getNewStepSize(double oldStepSize) const
{
    double factor = std::min(2., pow(getMaxErrorRatioSqr(), -1. / 6)); // heuristic multiplier

    return 0.9 * oldStepSize * factor;
}

// sets the new absolute precisity minPrecisity
void PrecisionHandler::setMinPrecisity(double value)
{
    minPrecisity = value;
    minPrecisitySq = value * value;
}

// sets the new dipole precisity dipolePrecisity
void PrecisionHandler::setDipolePrecisity(double value)
{
    if (value == 0)
    {
        use_dipoleCrit = false;
        dipolePrecisity = 1;
        dipolePrecisitySq = 1;
    }
    else
    {
        use_dipoleCrit = true;
        dipolePrecisity = value;
        minPrecisity /= value;

        dipolePrecisitySq = value * value;
        minPrecisitySq /= value * value;
    }
}

// returns the prevoiusly calculated maxErrorToleranceRatioSq, and due to the technical multiplication by dipolePrecisitySq, devides by dipolePrecisitySq
double PrecisionHandler::getMaxErrorRatioSqr() const
{
    return maxErrorToleranceRatioSq / dipolePrecisitySq;
}

// returns for which ID was the error the largest
unsigned int PrecisionHandler::maxErrorRatioID() const
{
    return selectedID;
}
