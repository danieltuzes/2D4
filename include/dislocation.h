//
// dislocation.h : contains function definitions, included in stress_protocol.h, simulation_data.h, simulation.h

/*
# 0.2
std::ofstream& sdddstCore::operator<<(std::ofstream& o, const DislwoB& disl) is added

# 0.1
First version tracked file
*/

#ifndef SDDDST_CORE_DISLOCATION_H
#define SDDDST_CORE_DISLOCATION_H

#define VERSION_dislocations 0.2

#include <string>
#include <fstream>

namespace sdddstCore
{
    // dislocation without Burger's vector, that can be deduced from the ID of the element in the vector
    class DislwoB
    {
    public:
        DislwoB() : x(0), y(0) {};
        DislwoB(double x, double y) : x(x), y(y) {};

        double x;
        double y;
        // Burgers vector will be deduced from the index of the dislocation

        friend std::ofstream& operator<<(std::ofstream& o, const DislwoB& disl);
    };

    // dislocation with ID, ID is in increasing from 0
    class DislwId
    {
    public:
        DislwId() : x(0), y(0), b(0), id(0) {};
        DislwId(double x, double y, double b, size_t id) : x(x), y(y), b(b), id(id) {};

        double x;
        double y;
        double b;
        size_t id;
    };

    // the original class defined by PGabor
    class Dislocation
    {
    public:
        Dislocation() : x(0), y(0), b(0) {};
        Dislocation(double x, double y, double b) : x(x), y(y), b(b) {};

        double x;
        double y;
        double b;
    };
}

#endif
