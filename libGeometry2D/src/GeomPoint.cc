#include "GeomPoint.hh"

#include <sstream>

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

MOTH_HOST
std::string moth_p2d::str(const moth_p2d& p1)
{
    std::stringstream stream;
    stream << p1;
    return stream.str();
}

MOTH_HOST MOTH_CORE
std::ostream& operator<<(std::ostream& stream, const moth_p2d& p)
{
    return stream << "(" << p.x << ", " << p.y << ")";
}
MOTH_HOST MOTH_CORE
std::istream& operator>>(std::istream& stream, moth_p2d& p)
{
    abort();
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

MOTH_HOST
std::string moth_p3d::str(const moth_p3d& p1)
{
    std::stringstream stream;
    stream << p1;
    return stream.str();
}

MOTH_HOST MOTH_CORE
std::ostream& operator<<(std::ostream& stream, const moth_p3d& p)
{
    return stream << "(" << p.x << ", " << p.y << ", " << p.z << ")";
}
MOTH_HOST MOTH_CORE
std::istream& operator>>(std::istream& stream, moth_p3d& p)
{
    abort();
}
