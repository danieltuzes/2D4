// 
// AnalyticField.h : contains function declarations for AnalyticField.cpp, also included in project_parser.cpp

/*changelog
# 0.3
xy and xy_diff_x are consts


# 0.2
implementation in case USE_IEEE_HYPERBOLIC is false defined

# 0.1
first version with #define VERSION_analytic_field
*/

#ifndef SDDDST_CORE_ANALYTIC_FIELD_H
#define SDDDST_CORE_ANALYTIC_FIELD_H

#define VERSION_analytic_field 0.3

#include "Fields/Field.h"
#include "constants.h"

#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>

namespace sdddstCore {

    class AnalyticField : public Field
    {
    public:
        AnalyticField();

        virtual double xy(double dx, double dy) const;
        virtual double xy_diff_x(double dx, double dy) const;
    };

}

#endif
