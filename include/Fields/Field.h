// 
// Field.h : contains function declaration for Field.cpp, also included in AnalyticField.h

/* 
# 0.7
functions xy and xy_diff_x takes additional arguments for precalculated hyperbolic and trigonometric function values, helps only if !USE_IEEE_HYPERBOLIC

# 0.6
* AnalyticField is merged into Field to simplify code by removing unused code
* AnalyticField and PeriodicShearStressElte are removed

# 0.1
First version tracked version of the file
*/

#ifndef SDDDST_CORE_FIELD_H
#define SDDDST_CORE_FIELD_H

#define VERSION_field 0.7

#include "constants.h"
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>

#include <string>
#include <memory>

namespace sdddstCore {

    class Field
    {
    public:
        Field();

        // returns the field at point dx, dy without knowing their hyperbolic and trigonometric function values
        double xy(double dx, double dy) const;

        // returns the derivative of the field in x direction at point dx, dy without knowing their hyperbolic and trigonometric function values
        double xy_diff_x(double dx, double dy) const;
        
        // returns the field at point dx, dy knowing their hyperbolic and trigonometric function values
        double xy(double dx, double dy, double exp2pix, double cos2piy) const;
        
        // returns the derivative of the field in x direction at point dx, dy knowing their hyperbolic and trigonometric function values
        double xy_diff_x(double dx, double dy, double exp2pix, double cos2piy) const;
    };

}

#endif
