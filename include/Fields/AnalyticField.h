
// AnalyticField.h : contains function declarations for AnalyticField.cpp, also included in project_parser.cpp

#ifndef SDDDST_CORE_ANALYTIC_FIELD_H
#define SDDDST_CORE_ANALYTIC_FIELD_H

#define VERSION_analytic_field 0.1

#include "Fields/Field.h"
#include "constants.h"

#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>

namespace sdddstCore {

class AnalyticField: public Field
{
public:
    AnalyticField();

    virtual double xy(double dx, double dy);
    virtual double xy_diff_x(double dx, double dy);
};

}

#endif
