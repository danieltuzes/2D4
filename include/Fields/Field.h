
// Field.h : contains function declaration for Field.cpp, also included in AnalyticField.h

#ifndef SDDDST_CORE_FIELD_H
#define SDDDST_CORE_FIELD_H

#define VERSION_field 0.1

#include <string>
#include <memory>

namespace sdddstCore {

    class Field
    {
    public:
        Field();
        virtual ~Field();

        virtual double xy(double dx, double dy) const;
        virtual double xy_diff_x(double dx, double dy) const;
    };

}

#endif
