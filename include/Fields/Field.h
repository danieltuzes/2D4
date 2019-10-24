// 
// Field.h : contains function declaration for Field.cpp, also included in AnalyticField.h

/* 
# 0.6
* AnalyticField is merged into Field to simplify code by removing unused code
* AnalyticField and PeriodicShearStressElte are removed

# 0.1
First version tracked version of the file
*/

#ifndef SDDDST_CORE_FIELD_H
#define SDDDST_CORE_FIELD_H

#define VERSION_field 0.6

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

        double xy(double dx, double dy) const;
        double xy_diff_x(double dx, double dy) const;
    };

}

#endif
