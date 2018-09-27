#include "libGeometry2D/src/GeomEdge.hh"
#include "libCommon/src/CommTest.hh"

#include <sstream>

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

bool geom_e2d::intersect(const geom_e2d& e1, geom_e2d e2, geom_e2d& e)
{
    geom_p2d v1 = e1.t - e1.s;
    geom_p2d v2 = e2.t - e2.s;
    geom_real_t det = geom_p2d::det(v1, v2);
    if(det != 0.0) {
        // Edges are not collinear.
        geom_p2d v = e2.s - e1.s;
        geom_real_t alpha = geom_p2d::det(v, v2) / det;
        geom_real_t gamma = geom_p2d::det(v, v1) / det;
        if ((0.0 < alpha && alpha <= 1.0) && (0.0 < gamma && gamma <= 1.0)) {
            // Edges intersect on point.
            e.s = e1.s + alpha * v1;
            e.t = e.s;
            return true;
        }
    } else {
        // Edges are collinear.
        if (geom_p2d::dot(v1, v2) < 0.0) {
            std::swap(e2.s, e2.t);
        }

        geom_p2d q1 = e2.s - e1.t;
        geom_p2d q2 = e2.t - e1.s;
        if (geom_p2d::det(q1, q2) == 0.0) {
            // Edges are on one line.
            if (geom_p2d::len(v1) + geom_p2d::len(v2) != fabs(geom_p2d::len(q2) - geom_p2d::len(q1))) {
                // Edges intersect on line.
                e.s = geom_p2d::max(geom_p2d::min(e2.s, e2.t), geom_p2d::min(e1.s, e1.t));
                e.t = geom_p2d::min(geom_p2d::max(e2.s, e2.t), geom_p2d::max(e1.s, e1.t));
                return true;
            }
        }
    }
    return false;
}

bool geom_e2d::reflect(const geom_e2d& e1, const geom_e2d& e2, geom_e2d& e)
{
    if (intersect(e1, e2, e) && e.s == e.t) {
        geom_p2d p = e.s;
        geom_real_t l1 = geom_p2d::len(e2.s - p);
        geom_real_t l2 = geom_p2d::len(e2.t - p);
        geom_p2d v = (e2.s - p) / l1;
        geom_p2d n = geom_p2d::normal(e1.s - p);

        geom_real_t cos = geom_p2d::dot(v, n);
        geom_real_t sin = geom_p2d::det(v, n);
        if (cos > 0.0) {
            geom_p2d q{n.x * cos - n.y * sin, n.x * sin + n.y * cos};
            e.s = e2.s;
            e.t = p + l2 * q;
        } else {
            e = e2;
        }
        return true;
    }
    e = e2;
    return false;
}

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
};
