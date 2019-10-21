
// precision_handler.h : contains the function declarations for precision_handler.cpp, also included in simulation.h

#ifndef SDDDST_CORE_PRECISION_HANDLER_H
#define SDDDST_CORE_PRECISION_HANDLER_H

#define VERSION_precision_handler 0.1

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

    private:
        std::vector<std::pair<double, double> > toleranceAndError;
        double minPrecisity; // set from position-precision
        double minPrecisitySqr; // square of minPrecisity (stored to save calculation time?)
        double maxErrorRatioSqr;
        unsigned int selectedID;
    };

}

#endif
