// 
// precision_handler.h : contains the function declarations for precision_handler.cpp, also included in simulation.h

/*
# 0.3
* use_dipoleCrit is added to tell if dipole precisity criterium should be used or not
* dipolePrecisity is to measure how close the two predictions should be for a pivot dislocation with respect to its nearest other dislocation distance multiplied with this number
* dipolePrecisitySq: dipolePrecisity * dipolePrecisity
* setDipolePrecisity modifies minPrecisity (devides with dipolePrecisity) and its square so that in the interaction calculation there will be no need for multiplications with the distance. But then when the error (square) is devided by the target precision (square), it must be compared with dipolePrecisity (square).
* to not to confuse names, simulation's calculateXError calls updateMaxErrorToleranceRatioSq, and getMaxErrorRatioSqr returns maxErrorToleranceRatioSq / dipolePrecisitySq

# 0.2
selectedID is returned maxErrorRatioID();

# 0.1: first version tracked file
*/

#ifndef SDDDST_CORE_PRECISION_HANDLER_H
#define SDDDST_CORE_PRECISION_HANDLER_H

#define VERSION_precision_handler 0.3

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

        void setSize(unsigned int size); // if use_dipoleCrit, sets the tolerance container to its proper size

        void reset(); // if use_dipoleCrit, resets the tolerance to minPrecisity, i.e. it forgets about the closest dislocations
        
        void updateTolerance(double distanceSq, unsigned int ID); // if use_dipoleCrit is not 0 then if distanceSq is smaller than the prescribed minPrecisitySq, the function overwrites it for that particle with distanceSq

        void updateMaxErrorToleranceRatioSq(double error, unsigned int ID); // calculates if error compared to the required precision (absolute or dipolePrecisity, depending on use_dipoleCrit) and stores the renitent's ID to selectedID

        double getNewStepSize(double oldStepSize) const; // heuristic guess for the new stepsize

        void setMinPrecisity(double value);     // sets the new absolute precisity minPrecisity
        void setDipolePrecisity(double value);  // sets the new dipole precisity dipolePrecisity
        
        // returns the prevoiusly calculated maxErrorToleranceRatioSq, and due to the technical multiplication by dipolePrecisitySq, devides by dipolePrecisitySq
        double getMaxErrorRatioSqr() const;

        // returns for which ID was the error the largest
        unsigned int maxErrorRatioID() const;


    private:
        // min(distance square, minPrecisitySq); the former toleranceAndError's second value was unnecessary, this is the remaining first value
        std::vector<double> tolerance;

        double minPrecisity;                // set from position-precision; to avoid numerous multiplications; technically, it is multiplied with dipolePrecisity but its ratio with dipolePrecisity is checked at the beginning of stage IV
        double minPrecisitySq;              // square of minPrecisity (stored to save calculation time?)
        double maxErrorToleranceRatioSq;    // the largest "error square divided by tolerance" for all particles
        unsigned int selectedID;            // the ID of the dislocation with the highest errorRatioSqr
        
        // dipoleCrit enables a feature to measure how close the two predictions should be for a pivot dislocation with respect to its nearest other dislocation distance multiplied with this number
        // it is used at stage IV whether or not to accept a step

        bool use_dipoleCrit;        // if the program should compare error with the closest dislocation too, or with minPrecisity only
        double dipolePrecisity;     // technically, this value multiplies minPrecisity, but its value is compared against dipolePrecisity, or to be more precise, their square
        double dipolePrecisitySq;   // dipolePrecisity * dipolePrecisity
    };

}

#endif
