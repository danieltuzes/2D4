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

PrecisionHandler::~PrecisionHandler()
{
    // Nothing to do
}

void PrecisionHandler::setSize(unsigned int size)
{
    if (toleranceAndError.size() < size)
        while (toleranceAndError.size() != size)
            toleranceAndError.emplace_back(minPrecisitySqr, 0);
    else
        toleranceAndError.resize(size);
}

unsigned long PrecisionHandler::getSize() const
{
    return toleranceAndError.size();
}

void PrecisionHandler::reset()
{
    maxErrorRatioSqr = 0;
    for (auto& i : toleranceAndError)
    {
        i.first = minPrecisitySqr;
        i.second = 0;
    }
}

void PrecisionHandler::updateTolerance(double distanceSqr, unsigned int ID) // const int& is more expensive than unsigned int
{
    double tmp = distanceSqr * 0.25 * 1e-2;
    if (tmp < minPrecisitySqr && tmp < toleranceAndError[ID].first)
    {
        if (tmp == 0)
        {
            toleranceAndError[ID].first = EPS;
            std::cout << "Two dislocations are in the same place!\n";
        }
        else
            toleranceAndError[ID].first = tmp;

    }
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

    double factor = std::min(2., 1. / pow(maxErrorRatioSqr, -1. / 6)); // heuristic multiplier (why?)

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

#ifdef BUILD_PYTHON_BINDINGS

std::string PrecisionHandler::__str__() const
{
    return "N: " + std::to_string(toleranceAndError.size()) +
        " Prec: " + std::to_string(minPrecisity) +
        " New stepsize factor: " + std::to_string(getNewStepSize(1.0));
}

std::string PrecisionHandler::__repr__() const
{
    return __str__();
}

#endif
