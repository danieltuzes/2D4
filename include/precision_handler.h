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

#ifndef SDDDST_CORE_PRECISION_HANDLER_H
#define SDDDST_CORE_PRECISION_HANDLER_H

#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
#include <iostream>

namespace sdddstCore {

class PrecisionHandler
{
public:
    PrecisionHandler();
    ~PrecisionHandler();

    void setSize(unsigned int size);
    unsigned long getSize() const;

    void reset();
    void updateTolerance(double distanceSqr, unsigned int ID);

    void updateError(double error, unsigned int ID);

    double getNewStepSize(double oldStepSize) const;

    double getMinPrecisity() const;
    void setMinPrecisity(double value);

    double getMaxErrorRatioSqr() const;

#ifdef BUILD_PYTHON_BINDINGS
    std::string __str__() const;
    std::string __repr__() const;
#endif

private:
    std::vector<std::pair<double, double> > toleranceAndError;
    double minPrecisity; // set from position-precision
    double minPrecisitySqr; // square of minPrecisity (stored to save calculation time?)
    double maxErrorRatioSqr;
    unsigned int selectedID;
};

}

#endif
