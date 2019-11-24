#include "dislocation.h"

std::ofstream& sdddstCore::operator<<(std::ofstream& o, const DislwoB& disl)
{
    o << disl.x << "\t" << disl.y;
    return o;
}
