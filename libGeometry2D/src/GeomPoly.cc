#include "libGeometry2D/src/GeomPoly.hh"
#include "libCommon/src/CommTest.hh"
#include "GeomPoly.hh"


#include <sstream>
#include <cstdarg>

#if defined (__APPLE__) && defined (__MACH__)
#undef __unix__
#define __unix__ 1
#endif  // defined (__APPLE__) && defined (__MACH__)
#ifdef __unix__
#include <unistd.h>
#else   // ifdef __unix__
#include <windows.h>
#endif  // ifdef __unix__

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

std::ostream& operator<<(std::ostream& stream, const moth_poly2d& poly)
{
    stream << "(";
    moth_poly2d_iter iter = poly.iter();
    do {
        stream << iter.point();
        if (iter.next() != poly.iter()) {
            stream << ", ";
        }
    } while ((++iter) != poly.iter());
    stream << ")";
    return stream;
}

std::istream& operator>>(std::istream& stream, moth_poly2d& poly)
{
    abort();
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#include <map>

struct moth_poly2d_event_t
{
    moth_poly2d_iter pE{};
    bool left{};
};  // struct moth_poly2d_event_t

struct moth_poly2d_event_queue : public std::vector<moth_poly2d_event_t>
{
public:
    explicit moth_poly2d_event_queue(const moth_poly2d& poly)
    {
        moth_poly2d_iter pE_cur = poly.iter();
        do {
            moth_poly2d_event_t eP_cur{pE_cur};
            moth_poly2d_event_t eP_nxt{pE_cur.next()};
            if (eP_cur.pE.point() < eP_nxt.pE.point()) {
                eP_cur.left = true;
            } else {
                eP_nxt.left = true;
            }
            push_back(eP_cur);
            push_back(eP_nxt);
        } while ((++pE_cur) != poly.iter());
        std::sort(begin(), end(), [](const moth_poly2d_event_t& eP1, const moth_poly2d_event_t& eP2) {
            return eP1.pE.point() < eP2.pE.point();
        });
    }
};  // struct moth_poly2d_event_queue

struct moth_poly2d_sweep_line_segment_t
{

};  // struct moth_poly2d_sweep_line_segment_t






bool moth_poly2d::simple(const moth_poly2d& poly)
{
    return false;
}







































// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

MOTH_HOST
std::string moth_poly2d::str(const moth_poly2d& poly)
{
    std::stringstream stream;
    stream << poly;
    std::string poly_str = stream.str();
    return poly_str;
}

MOTH_HOST
std::string moth_poly2d::plt(const moth_poly2d& poly)
{
    std::stringstream plt_stream;
    moth_poly2d_iter iter = poly.iter();
    do {
        plt_stream << iter.point().x << " " << iter.point().y << " ";
        //plt_stream << iter.point().u << " " << iter.point().v << " ";
        plt_stream << std::endl;
    } while ((++iter) != poly.iter());
    plt_stream << iter.point().x << " " << iter.point().y << " ";
    //plt_stream << iter.point().u << " " << iter.point().v << " ";
    plt_stream << std::endl;
    plt_stream << "e";
    plt_stream << std::endl;
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

    CHAR plt_cmd[] = "gnuplot.exe -p -e \"set xrange[0:10]; set yrange[0:2]; plot '-' with vectors\"";
    if (CreateProcessA(nullptr,
                       plt_cmd,
                       nullptr, nullptr, TRUE, 0, nullptr, nullptr,
                       &plt_startup_info,
                       &plt_process_info)) {

        WriteFile(plt_pipe_input_write, plt_input.data(), plt_input.size() * sizeof(plt_input[0]), nullptr, nullptr);
        CloseHandle(plt_pipe_input_write);

        WaitForSingleObject(plt_process_info.hProcess, INFINITE);
        CloseHandle(plt_process_info.hProcess);
        CloseHandle(plt_process_info.hThread);
    }
#endif
    return "";
}

// >>>>----------------------------------------------------------------------------<<<< //
// >>>>----------------------------------------------------------------------------<<<< //

COMM_UNIT_TEST()
{
    moth_poly2d square = geom_poly2d_primitives::rect({0.0,0.0}, {1.0,1.0});
    COMM_UNIT_VERIFY_T(moth_poly2d::area(square) == 1.0);
    COMM_UNIT_VERIFY_T(moth_poly2d::len(square) == 4.0);

    moth_poly2d circle = geom_poly2d_primitives::circle({0.0,0.0}, 1, 4);
    COMM_UNIT_VERIFY_T(moth_poly2d::area(circle) == 2.0);
};

COMM_UNIT_TEST()
{
    moth_poly2d square = geom_poly2d_primitives::rect({0.0,0.0}, {1.0,1.0});
    COMM_UNIT_VERIFY_T(moth_poly2d::str(square) ==
                       "((0, 0), (1, 0), (1, 1), (0, 1))");
};

#if 0
COMM_UNIT_TEST()
{
    moth_poly2d box1;
    moth_poly2d box2;

    /* Simplest possible intersection.
     * ( And with reversed second box. )*/
    box1 = geom_poly2d_primitives::rect({0.0,0.0}, {2.0,2.0});
    box2 = geom_poly2d_primitives::rect({1.0,1.0}, {3.0,3.0});
    COMM_UNIT_VERIFY_T(moth_poly2d::str(box1 + box2) ==
                       "((0, 0), (2, 0), (2, 1), (3, 1),"
                       " (3, 3), (1, 3), (1, 2), (0, 2))");
    box1 = geom_poly2d_primitives::rect({0.0,0.0}, {2.0,2.0});
    box2 = geom_poly2d_primitives::rect({3.0,3.0}, {1.0,1.0});
    COMM_UNIT_VERIFY_T(moth_poly2d::str(box1 + box2) ==
                       "((0, 0), (2, 0), (2, 1), (3, 1),"
                       " (3, 3), (1, 3), (1, 2), (0, 2))");

    box1 = geom_poly2d_primitives::rect({0.0,0.0}, {4.0,4.0});
    box2 = geom_poly2d_primitives::rect({2.0,1.0}, {6.0,3.0});
    COMM_UNIT_VERIFY_T(moth_poly2d::str(box1 + box2) ==
                       "((0, 0), (4, 0), (4, 1), (6, 1),"
                       " (6, 3), (4, 3), (4, 4), (0, 4))");
    box1 = geom_poly2d_primitives::rect({0.0,0.0}, {4.0,4.0});
    box2 = geom_poly2d_primitives::rect({6.0,3.0}, {2.0,1.0});
    COMM_UNIT_VERIFY_T(moth_poly2d::str(box1 + box2) ==
                       "((0, 0), (4, 0), (4, 1), (6, 1),"
                       " (6, 3), (4, 3), (4, 4), (0, 4))");
};
#endif
