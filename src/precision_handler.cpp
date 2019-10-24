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