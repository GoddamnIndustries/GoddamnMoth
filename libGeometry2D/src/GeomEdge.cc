#include "libGeometry2D/src/GeomEdge.hh"
#include "libCommon/src/CommTest.hh"

#include <sstream>

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

std::ostream& operator<<(std::ostream& stream, const geom_e2d& e)
{
    stream << "(" << e.s << ", " << e.t << ")";
    return stream;
}

std::istream& operator>>(std::istream& stream, geom_e2d& e)
{
    std::abort();
}

std::string geom_e2d::str(const geom_e2d& e)
{
    std::stringstream stream;
    stream << e;
    return stream.str();
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

std::ostream& operator<<(std::ostream& stream, const geom_e2d_list* poly)
{
    stream << "(";
    const geom_e2d_list* head = poly;
    do {
        stream << poly->point;
        if (poly->next != head) {
            stream << ", ";
        }
    } while (geom_e2d_list::move(poly, head));
    stream << ")";
    return stream;
}

std::istream& operator>>(std::istream& stream, geom_e2d_list* poly)
{
    abort();
}

std::string geom_e2d_list::str(const geom_e2d_list* poly)
{
    std::stringstream stream;
    stream << poly;
    return stream.str();
}

// >>>>----------------------------------------------------------------------------<<<< //
// >>>>----------------------------------------------------------------------------<<<< //

COMM_UNIT_TEST()
{
    geom_e2d e;

    // Intersection of the non-collinear edges.
    COMM_UNIT_VERIFY_F(geom_e2d::intersect({{0.0,5.0}, {2.0,5.0}}, {{2.0, 0.0}, {2.0, 2.0}}, e));
    COMM_UNIT_VERIFY_T(geom_e2d::intersect({{0.0,0.0}, {2.0,2.0}}, {{0.0, 2.0}, {2.0, 0.0}}, e) &&
                       e.s == geom_p2d{1.0, 1.0} && e.t == e.s);

    // Intersection of the non-collinear on the ends.
    // ( Note that edge is e is (s, t], not [s, t]. )
    COMM_UNIT_VERIFY_F(geom_e2d::intersect({{0.0,0.0}, {2.0,0.0}}, {{2.0, 0.0}, {2.0, 2.0}}, e));
    COMM_UNIT_VERIFY_T(geom_e2d::intersect({{0.0,0.0}, {2.0,0.0}}, {{2.0, 2.0}, {2.0, 0.0}}, e) &&
                       e.s == geom_p2d{2.0, 0.0} && e.t == e.s);

    // Intersection of the collinear edges not on the single line.
    COMM_UNIT_VERIFY_F(geom_e2d::intersect({{0.0,2.0}, {2.0,0.0}}, {{6.0, 2.0}, {4.0, 4.0}}, e));

    // Intersection of the collinear edges on single line.
    COMM_UNIT_VERIFY_F(geom_e2d::intersect({{1.0,0.0}, {4.0,0.0}}, {{6.0, 0.0}, {5.0, 0.0}}, e));
    COMM_UNIT_VERIFY_T(geom_e2d::intersect({{1.0,0.0}, {4.0,0.0}}, {{5.0, 0.0}, {3.0, 0.0}}, e) &&
                       e == geom_e2d{{3.0,0.0}, {4.0,0.0}});
    COMM_UNIT_VERIFY_T(geom_e2d::intersect({{1.0,0.0}, {4.0,0.0}}, {{3.0, 0.0}, {2.0, 0.0}}, e) &&
                       e == geom_e2d{{2.0,0.0}, {3.0,0.0}});
    COMM_UNIT_VERIFY_T(geom_e2d::intersect({{1.0,0.0}, {4.0,0.0}}, {{0.0, 0.0}, {2.0, 0.0}}, e) &&
                       e == geom_e2d{{1.0,0.0}, {2.0,0.0}});
};

COMM_UNIT_TEST()
{
    geom_e2d e;

    // Reflection with the CCW normal.
    // ( e2 is a diode wall: in first case particle simply files through it, reflects in second. )
    COMM_UNIT_VERIFY_F(geom_e2d::reflect({{0.0,2.0}, {4.0,2.0}}, {{0.0,4.0}, {1.0,3.0}}, e));
    COMM_UNIT_VERIFY_T(geom_e2d::reflect({{0.0,2.0}, {4.0,2.0}}, {{0.0,4.0}, {4.0,0.0}}, e) &&
                       e == geom_e2d{{0.0,4.0}, {4.0, 0.0}});
    COMM_UNIT_VERIFY_T(geom_e2d::reflect({{0.0,2.0}, {4.0,2.0}}, {{4.0,0.0}, {1.0,3.0}}, e) &&
                       e == geom_e2d{{4.0,0.0}, {1.0, 1.0}});

    COMM_UNIT_VERIFY_F(geom_e2d::reflect({{2.0,0.0}, {0.0,0.0}}, {{1.0,0.0}, {1.0,1.0}}, e));
    COMM_UNIT_VERIFY_T(geom_e2d::reflect({{2.0,0.0}, {0.0,0.0}}, {{1.0,1.0}, {1.0,0.0}}, e) &&
                       e == geom_e2d{{1.0,1.0}, {1.0,0.0}});
};
