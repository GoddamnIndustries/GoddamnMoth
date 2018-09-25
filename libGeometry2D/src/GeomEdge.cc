#include "GeomEdge.hh"
#include "../../libLinearAlgebra/src/Test.h"

#include <sstream>
#include <cstdarg>
#if (defined (__APPLE__) && defined (__MACH__))
#define __unix__ 1
#endif
#ifdef __unix__
#include <unistd.h>
#else
#include <windows.h>
#endif

#undef min
#undef max

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

std::ostream& operator<<(std::ostream& stream, const geom_p2d& p)
{
    return stream << "(" << p.x << ", " << p.y << ")";
}
std::istream& operator>>(std::istream& stream, geom_p2d& p)
{
    return stream;
}

std::string geom_p2d::str(const geom_p2d& p)
{
    std::stringstream stream;
    stream << p;
    return stream.str();
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
    abort();
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

std::string geom_e2d_list::plt(const geom_e2d_list* poly, ...)
{
    std::stringstream plt_stream;
    va_list poly_list{};
    va_start(poly_list, poly);
    const geom_e2d_list* head = poly;
    do {
        plt_stream << poly->point.x << " " << poly->point.y << " ";
        plt_stream << std::endl;
    } while (geom_e2d_list::move(poly, head));
    plt_stream << poly->point.x << " " << poly->point.y << " ";
    plt_stream << std::endl;
    plt_stream << "e";
    plt_stream << std::endl;
    va_end(poly_list);
    std::string plt_input = plt_stream.str();

#ifdef __unix__
    int plt_pipe[2] = {};
    pipe(plt_pipe);

    pid_t plt_child = fork();
    if (plt_child == 0) {
        close(plt_pipe[1]);
        dup2(plt_pipe[0], STDIN_FILENO);
        close(plt_pipe[0]);

        const char* plt_exe = "/usr/local/bin/gnuplot";
        const char* plt_args[] = {plt_exe, "-e", "set xrange[-3:3]; set yrange[-3:3]; plot '-' with lines; pause -1;", nullptr};
        execvp(plt_exe, const_cast<char* const*>(plt_args));
        exit(-1);
    } else {
        close(plt_pipe[0]);
        write(plt_pipe[1], plt_input.data(), plt_input.size() * sizeof(plt_input[0]));
        close(plt_pipe[1]);

        int plt_status = 0;
        waitpid(plt_child, &plt_status, 0);
    }
#else
    SECURITY_ATTRIBUTES plt_security_attrs{};
    plt_security_attrs.nLength = sizeof(plt_security_attrs);
    plt_security_attrs.bInheritHandle = TRUE;

    HANDLE plt_pipe_input_read{};
    HANDLE plt_pipe_input_write{};
    CreatePipe(&plt_pipe_input_read, &plt_pipe_input_write, &plt_security_attrs, 0);
    SetHandleInformation(plt_pipe_input_write, HANDLE_FLAG_INHERIT, 0);

    PROCESS_INFORMATION plt_process_info{};
    STARTUPINFO plt_startup_info{};
    plt_startup_info.cb = sizeof(plt_startup_info);
    plt_startup_info.hStdInput = plt_pipe_input_read;
    plt_startup_info.dwFlags |= STARTF_USESTDHANDLES;

    CHAR plt_cmd[] = "gnuplot.exe -e \"set xrange[-3:3]; set yrange[-3:3]; plot '-' with lp; pause -1;\"";
    if (CreateProcessA(nullptr,
                       plt_cmd,
                       nullptr, nullptr, TRUE, 0, nullptr, nullptr,
                       &plt_startup_info,
                       &plt_process_info)) {
        CloseHandle(plt_process_info.hProcess);
        CloseHandle(plt_process_info.hThread);
        WriteFile(plt_pipe_input_write, plt_input.data(), plt_input.size() * sizeof(plt_input[0]), nullptr, nullptr);
        DebugBreak();
        CloseHandle(plt_pipe_input_write);
    }
#endif
    return "";
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

void
geom_collide(geom_e2d e1, geom_e2d e2, geom_c2d_list* collision)
{
    geom_real_t det;
    geom_p2d v1 = e1.t - e1.s;
    geom_p2d v2 = e2.t - e2.s;
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
        if (geom_p2d::dot(v1, v2) < 0) {
            std::swap(e2.s, e2.t);
            v2 = -v2;
        }

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

#if 0
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
#endif

