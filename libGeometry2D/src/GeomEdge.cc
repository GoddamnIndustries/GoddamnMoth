#include "GeomEdge.hh"
//#include "../../libLinearAlgebra/src/Test.h"

template <typename T>
int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/// @todo Move to header one day.
struct geom_e2d_list : public geom_p2d
{
    geom_e2d_list* next = this;
    geom_e2d_list* next_in = nullptr;
    geom_e2d_list* next_out = nullptr;
    geom_e2d edge() const { return { *this, *next }; }

    void insert(const geom_p2d& p)
    {
        auto* n = new geom_e2d_list{};
        (geom_p2d&)*n = p;
        n->next = next;
        next = n;
    }


    static int orientation(const geom_e2d_list* E);

    static bool contains(const geom_e2d_list* E, const geom_p2d& p);

    static void clip(geom_e2d_list* E1, geom_e2d_list* E2);
};  // struct geom_e2d_list

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

int
geom_e2d_list::orientation(const geom_e2d_list* E)
{
    geom_real_t E_square = 0.0;
    const geom_e2d_list* E_i = E;
    do
    {
        const geom_e2d_list* E_ip = E->next;
        E_square += geom_p2d::det(*E_i, *E_ip);
    } while ((E_i = E_i->next) != E);
    return sgn(E_square);
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

void
geom_collide(const geom_e2d& e1, const geom_e2d& e2, geom_c2d_list* collision)
{
    geom_real_t det;
    const geom_p2d v1 = e1.t - e1.s;
    const geom_p2d v2 = e2.t - e2.s;
    det = geom_p2d::det(v1, v2);
    if(det != 0.0)
    {
        // Edges are not collinear.
        const geom_p2d b = e2.s - e1.s;
        const geom_real_t alpha = geom_p2d::det(b, v2) / det;
        const geom_real_t gamma = geom_p2d::det(b, v1) / det;
        if (alpha >= 0 && alpha <= 1 && gamma >= 0 && gamma <= 1)
        {
            // Edges intersect on point.
            collision->e.s = e1.s + alpha * v1;
            collision->e.t = collision->e.s;
            collision->type = GEOM_C2D_INTERSECT;

            // Edges intersection direction.
            const geom_p2d n1 = collision->e.s - e1.s;
            const geom_p2d n2 = collision->e.s - e2.s;
            if (geom_p2d::dot(n1, geom_p2d::normal(n2)) > 0.0)
            {
                collision->direction = GEOM_C2D_DIRECTION_IN;
            }
            else
            {
                collision->direction = GEOM_C2D_DIRECTION_OUT;
            }
        }
    }
    else
    {
        // Edges are collinear.
        const geom_p2d q1 = e2.s - e1.t;
        const geom_p2d q2 = e2.t - e1.s;
        det = geom_p2d::det(q1, q2);
        if (det == 0.0)
        {
            // Edges are on one line.
            if (geom_p2d::len(v1) + geom_p2d::len(v2) != fabs(geom_p2d::len(q2) - geom_p2d::len(q1)))
            {
                collision->e.s = geom_p2d::max(geom_p2d::min(e2.s, e2.t), geom_p2d::min(e1.s, e1.t));
                collision->e.t = geom_p2d::min(geom_p2d::max(e2.s, e2.t), geom_p2d::max(e1.s, e1.t));
                collision->direction = GEOM_C2D_DIRECTION_NONE;
                if (geom_p2d::len(q1) * geom_p2d::len(q2) != 0.0)
                {
                    collision->type = GEOM_C2D_TOUCH;
                }
                else
                {
                    collision->type = GEOM_C2D_INTERSECT;
                }
            }
        }
    }
}

void
geom_collide(const geom_e2d_list* E1, const geom_e2d& e2, geom_c2d_list* collision)
{
    // Collide each edge of the E1 with e2.
    const geom_e2d_list* E1_i = E1;
    do
    {
        const geom_e2d e1 = E1_i->edge();
        geom_collide(e1, e2, collision);
        if (collision->type != GEOM_C2D_NONE)
        {
            if (collision->type == GEOM_C2D_INTERSECT)
            {
                const geom_p2d v1 = collision->e.s - e1.s;
                const geom_p2d v2 = collision->e.s - e1.s;
                if (geom_p2d::dot(v1, geom_p2d::normal(v2)) > 0.0)
                {
                    collision->direction = GEOM_C2D_DIRECTION_IN;
                }
                else
                {
                    collision->direction = GEOM_C2D_DIRECTION_OUT;
                }
            }

            // Create new node for the next collision.
            collision->next = new geom_c2d_list();
            collision = collision->next;
        }
    } while ((E1_i = E1_i->next) != E1);
}

void
geom_collide(const geom_e2d_list* E1, const geom_e2d_list* E2, geom_c2d_list* collision)
{
    // Collide each edge of the E2 with E1.
    const geom_e2d_list* E2_j = E2;
    do
    {
        const geom_e2d e2 = E2_j->edge();
        geom_collide(E1, e2, collision);

        // Skip to the last entry of the collision.
        while ((collision = collision->next)->type != GEOM_C2D_NONE);
    } while ((E2_j = E2_j->next) != E2);
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

void
geom_e2d_list::clip(geom_e2d_list* E1, geom_e2d_list* E2)
{
    geom_e2d_list* E1_i = E1;
    do
    {
        geom_e2d_list* E2_j = E2;
        do
        {
        } while ((E2_j = E2_j->next) != E2);
    } while ((E1_i = E1_i->next) != E1);
}


