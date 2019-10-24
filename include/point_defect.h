
// point_defect.h : contains the function declarataion for simulation_data.cpp via simulation_data.h; also included in project_parser.h

#ifndef SDDDST_CORE_POINT_DEFECT_H
#define SDDDST_CORE_POINT_DEFECT_H

#include <string>

#define VERSION_point_defect 0.1

namespace sdddstCore {

    struct PointDefect
    {
        double x;
        double y;

        bool operator==(PointDefect a)
        {
            if (a.x == x && a.y == y)
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    };

}

#endif
