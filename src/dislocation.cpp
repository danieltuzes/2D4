#include "dislocation.h"

namespace sdddstCore
{
    std::ofstream& operator<<(std::ofstream& o, const DislwoB& disl)
    {
        o << disl.x << "\t" << disl.y;
        return o;
    }
}