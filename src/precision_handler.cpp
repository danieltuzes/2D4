// 
// precision_handler.cpp : contains the function defintions for precision_handler.h, also included in simulation.h

#include "constants.h"
#include "precision_handler.h"

using namespace sdddstCore;

PrecisionHandler::PrecisionHandler() :
    minPrecisity(1e-6),
    minPrecisitySqr(1e-12),
    maxErrorRatioSqr(0),
    selectedID(0)
{
    //Nothing to do
}

void PrecisionHandler::setSize(unsigned int size)
{
    toleranceAndError.resize(size, std::pair<double, double>(minPrecisitySqr, 0));
}

void PrecisionHandler::reset()
{
    maxErrorRatioSqr = 0;
    toleranceAndError.assign(
        toleranceAndError.size(),
        std::pair<double, double>(minPrecisitySqr, 0)
    );
}

// if distanceSQ_fraction is smaller than the prescribed minPrecisitySqr, the function overwrites it for that particle with distanceSQ_fraction
void PrecisionHandler::updateTolerance(double distanceSQ_fraction, unsigned int ID)
{
    if (distanceSQ_fraction < minPrecisitySqr &&
        distanceSQ_fraction < toleranceAndError[ID].first)
        toleranceAndError[ID].first = distanceSQ_fraction;
}

void PrecisionHandler::updateError(double error, unsigned int ID)
{
    if (toleranceAndError[ID].second < error)
    {
        toleranceAndError[ID].second = error;
        double tmp = error * error / toleranceAndError[ID].first;
        if (tmp > maxErrorRatioSqr)
        {
            maxErrorRatioSqr = tmp;
            selectedID = ID;
        }
    }
}

double PrecisionHandler::getNewStepSize(double oldStepSize) const
{
    if (0 == maxErrorRatioSqr)
        return oldStepSize * 2;

    double factor = std::min(2., pow(maxErrorRatioSqr, -1. / 6)); // heuristic multiplier (why?)

    return 0.9 * oldStepSize * factor;
}

double PrecisionHandler::getMinPrecisity() const
{
    return minPrecisity;
}

void PrecisionHandler::setMinPrecisity(double value)
{
    minPrecisity = value;
    minPrecisitySqr = value * value;
}

double PrecisionHandler::getMaxErrorRatioSqr() const
{
    return maxErrorRatioSqr;
}