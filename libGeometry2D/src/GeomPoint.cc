#include "GeomPoint.hh"

#include <sstream>

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

GEOM_HOST
std::ostream& operator<<(std::ostream& stream, const geom_p2d& p)
{
    return stream << "(" << p.x << ", " << p.y << ")";
}
GEOM_HOST
std::istream& operator>>(std::istream& stream, geom_p2d& p)
{
    abort();
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

GEOM_HOST
std::string geom_p2d::str(const geom_p2d& p1)
{
    std::stringstream stream;
    stream << p1;
    return stream.str();
}
