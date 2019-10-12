
// dislocation.h : contains function definitions, included in stress_protocol.h, simulation_data.h, simulation.h

#ifndef SDDDST_CORE_DISLOCATION_H
#define SDDDST_CORE_DISLOCATION_H

#define VERSION_dislocations 0.1

#include <string>

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
#ifdef BUILD_PYTHON_BINDINGS
        std::string __str__() const
        {
            return "x: " + std::to_string(this->x) +
                " y: " + std::to_string(this->y) +
                " b: " + std::to_string(this->b);
        }
        std::string __repr__() const
        {
            return __str__();
        }
#endif
    };
}

#endif
