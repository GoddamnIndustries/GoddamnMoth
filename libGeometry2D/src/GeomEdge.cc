#include "libGeometry2D/src/GeomEdge.hh"
#include "libCommon/src/CommTest.hh"

#include <sstream>

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

MOTH_HOST
std::string moth_e2d::str(const moth_e2d& e)
{
    std::stringstream stream;
    stream << e;
    return stream.str();
}

MOTH_HOST MOTH_CORE
std::ostream& operator<<(std::ostream& stream, const moth_e2d& e)
{
    stream << "(" << e.p1 << ", " << e.p2 << ")";
    return stream;
}

MOTH_HOST MOTH_CORE
std::istream& operator>>(std::istream& stream, moth_e2d& e)
{
    std::abort();
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

MOTH_HOST
std::string moth_e3d::str(const moth_e3d& e)
{
    std::stringstream stream;
    stream << e;
    return stream.str();
}

MOTH_HOST MOTH_CORE
std::ostream& operator<<(std::ostream& stream, const moth_e3d& e)
{
    stream << "(" << e.p1 << ", " << e.p2 << ")";
    return stream;
}

MOTH_HOST MOTH_CORE
std::istream& operator>>(std::istream& stream, moth_e3d& e)
{
    std::abort();
}

// >>>>----------------------------------------------------------------------------<<<< //
// >>>>----------------------------------------------------------------------------<<<< //

COMM_UNIT_TEST()
{
    moth_e2d e;

    // Intersection of the non-collinear edges.
    COMM_UNIT_VERIFY_F(moth_e2d::intersect({{0.0,5.0}, {2.0,5.0}}, {{2.0, 0.0}, {2.0, 2.0}}, e));
    COMM_UNIT_VERIFY_T(moth_e2d::intersect({{0.0,0.0}, {2.0,2.0}}, {{0.0, 2.0}, {2.0, 0.0}}, e) &&
                       e.p1 == moth_p2d{1.0, 1.0} && e.p2 == e.p1);

    // Intersection of the non-collinear on the ends.
    // ( Note that edge is e is (s, t], not [s, t]. )
    COMM_UNIT_VERIFY_F(moth_e2d::intersect({{0.0,0.0}, {2.0,0.0}}, {{2.0, 0.0}, {2.0, 2.0}}, e));
    COMM_UNIT_VERIFY_T(moth_e2d::intersect({{0.0,0.0}, {2.0,0.0}}, {{2.0, 2.0}, {2.0, 0.0}}, e) &&
                       e.p1 == moth_p2d{2.0, 0.0} && e.p2 == e.p1);

    // Intersection of the collinear edges not on the single line.
    COMM_UNIT_VERIFY_F(moth_e2d::intersect({{0.0,2.0}, {2.0,0.0}}, {{6.0, 2.0}, {4.0, 4.0}}, e));

    // Intersection of the collinear edges on single line.
    COMM_UNIT_VERIFY_F(moth_e2d::intersect({{1.0,0.0}, {4.0,0.0}}, {{6.0, 0.0}, {5.0, 0.0}}, e));
    COMM_UNIT_VERIFY_T(moth_e2d::intersect({{1.0,0.0}, {4.0,0.0}}, {{5.0, 0.0}, {3.0, 0.0}}, e) &&
                       e == moth_e2d{{3.0,0.0}, {4.0,0.0}});
    COMM_UNIT_VERIFY_T(moth_e2d::intersect({{1.0,0.0}, {4.0,0.0}}, {{3.0, 0.0}, {2.0, 0.0}}, e) &&
                       e == moth_e2d{{2.0,0.0}, {3.0,0.0}});
    COMM_UNIT_VERIFY_T(moth_e2d::intersect({{1.0,0.0}, {4.0,0.0}}, {{0.0, 0.0}, {2.0, 0.0}}, e) &&
                       e == moth_e2d{{1.0,0.0}, {2.0,0.0}});
};

#if 0
COMM_UNIT_TEST()
{
    moth_e2d e;

    // Reflection with the CCW normal.
    // ( e2 is a diode wall: in first case particle simply files through it, reflects in second. )
    COMM_UNIT_VERIFY_F(moth_e2d::reflect({{0.0,2.0}, {4.0,2.0}}, {{0.0,4.0}, {1.0,3.0}}, e));
    COMM_UNIT_VERIFY_T(moth_e2d::reflect({{0.0,2.0}, {4.0,2.0}}, {{0.0,4.0}, {4.0,0.0}}, e) &&
                       e == moth_e2d{{0.0,4.0}, {4.0, 0.0}});
    COMM_UNIT_VERIFY_T(moth_e2d::reflect({{0.0,2.0}, {4.0,2.0}}, {{4.0,0.0}, {1.0,3.0}}, e) &&
                       e == moth_e2d{{4.0,0.0}, {1.0, 1.0}});

    COMM_UNIT_VERIFY_F(moth_e2d::reflect({{2.0,0.0}, {0.0,0.0}}, {{1.0,0.0}, {1.0,1.0}}, e));
    COMM_UNIT_VERIFY_T(moth_e2d::reflect({{2.0,0.0}, {0.0,0.0}}, {{1.0,1.0}, {1.0,0.0}}, e) &&
                       e == moth_e2d{{1.0,1.0}, {1.0,0.0}});
};
#endif
