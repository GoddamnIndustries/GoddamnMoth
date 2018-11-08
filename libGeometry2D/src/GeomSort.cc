#include "libGeometry2D/src/GeomSort.hh"

#include <fstream>
#include <omp.h>

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

#define MOTH_OMP 1

MOTH_HOST MOTH_CORE
template<typename moth_task1_t, typename moth_task2_t>
inline void moth_omp_tasks(bool cond, moth_task1_t&& task1,
                                      moth_task2_t&& task2)
{
    if (cond) {
        {
#pragma omp task
            {
                task1();
            }
#pragma omp task
            {
                task2();
            }
        }
    } else {
        task1();
        task2();
    }
}

MOTH_HOST MOTH_CORE
inline moth_size_t moth_omp_get_num_threads()
{
    return static_cast<moth_size_t>(omp_get_num_threads());
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * Recurring 2D Hilbert sort of the given point set
 * starting with the given quadrant -- O(n*log(n)).
 */
MOTH_HOST MOTH_CORE
static void moth_sort_recursive(moth_p2d* pP_beg,
                                moth_p2d* pP_end,
                                const moth_p2d& p_min,
                                const moth_p2d& p_max,
                                moth_size_t nOrient,
                                moth_size_t nThreads)
{
    assert(pP_beg != nullptr);
    assert(pP_end != nullptr);
    if (pP_beg + 1 < pP_end) {
        moth_p2d p_cnr{0.5 * (p_min + p_max)};

        /* Separate quoters based on the orientation. */
        moth_p2d* pP1{pP_beg};
        moth_p2d* pP3{};
        switch (nOrient) {
            /* 4-|3
             * 1-|2 Orientation. */
            case 0:
                /* Separate lower and upper halves. */
                pP3 = pP1;
                for (moth_p2d* pP3_end{pP_end}; pP3 != pP3_end;) {
                    if (pP3->y > p_cnr.y) {
                        std::swap(*pP3, *--pP3_end);
                    } else {
                        ++pP3;
                    }
                }
                moth_omp_tasks(nThreads >= 2, [=]() {
                    /* Separate left and right quadrants of the lower half
                     * and recursively process the quadrants. */
                    moth_p2d* pP2 = pP1;
                    for (moth_p2d* pP2_end{pP3}; pP2 != pP2_end;) {
                        if (pP2->x > p_cnr.x) {
                            std::swap(*pP2, *--pP2_end);
                        } else {
                            ++pP2;
                        }
                    }
                    moth_size_t nThreads2 = moth_fdiv2(nThreads);
                    moth_omp_tasks(nThreads2 >= 2, [=]() {
                        moth_size_t nThreads4 = moth_fdiv2(nThreads2);
                        moth_sort_recursive(pP1, pP2, {p_min.x, p_min.y}, {p_cnr.x, p_cnr.y}, 1, nThreads4);
                    }, [=]() {
                        moth_size_t nThreads4 = moth_cdiv2(nThreads2);
                        moth_sort_recursive(pP2, pP3, {p_cnr.x, p_min.y}, {p_max.x, p_cnr.y}, 0, nThreads4);
                    });
                }, [=]() {
                    /* Separate right and left quadrants of the upper half
                     * and recursively process the quadrants. */
                    moth_p2d* pP4 = pP3;
                    for (moth_p2d* pP4_end{pP_end}; pP4 != pP4_end;) {
                        if (pP4->x < p_cnr.x) {
                            std::swap(*pP4, *--pP4_end);
                        } else {
                            ++pP4;
                        }
                    }
                    moth_size_t nThreads2 = moth_cdiv2(nThreads);
                    moth_omp_tasks(nThreads2 >= 2, [=]() {
                        moth_size_t nThreads4 = moth_fdiv2(nThreads2);
                        moth_sort_recursive(pP3, pP4, {p_cnr.x, p_cnr.y}, {p_max.x, p_max.y}, 0, nThreads4);
                    }, [=]() {
                        moth_size_t nThreads4 = moth_cdiv2(nThreads2);
                        moth_sort_recursive(pP4, pP_end, {p_min.x, p_cnr.y}, {p_cnr.x, p_max.y}, 2, nThreads4);
                    });
                });
                break;

            /* 2--3
             * 1||4 Orientation. */
            case 1:
                /* Separate left and right halves. */
                pP3 = pP1;
                for (moth_p2d* pP3_end{pP_end}; pP3 != pP3_end;) {
                    if (pP3->x > p_cnr.x) {
                        std::swap(*pP3, *--pP3_end);
                    } else {
                        ++pP3;
                    }
                }
                moth_omp_tasks(nThreads >= 2, [=]() {
                    /* Separate lower and upper quadrants of the left half
                     * and recursively process the quadrants. */
                    moth_p2d* pP2 = pP1;
                    for (moth_p2d* pP2_end{pP3}; pP2 != pP2_end;) {
                        if (pP2->y > p_cnr.y) {
                            std::swap(*pP2, *--pP2_end);
                        } else {
                            ++pP2;
                        }
                    }
                    moth_size_t nThreads2 = moth_fdiv2(nThreads);
                    moth_omp_tasks(nThreads2 >= 2, [=]() {
                        moth_size_t nThreads4 = moth_fdiv2(nThreads2);
                        moth_sort_recursive(pP1, pP2, {p_min.x, p_min.y}, {p_cnr.x, p_cnr.y}, 0, nThreads4);
                    }, [=]() {
                        moth_size_t nThreads4 = moth_cdiv2(nThreads2);
                        moth_sort_recursive(pP2, pP3, {p_min.x, p_cnr.y}, {p_cnr.x, p_max.y}, 1, nThreads4);
                    });
                }, [=]() {
                    /* Separate upper and lower quadrants of the right half
                     * and recursively process the quadrants. */
                    moth_p2d* pP4 = pP3;
                    for (moth_p2d* pP4_end{pP_end}; pP4 != pP4_end;) {
                        if (pP4->y < p_cnr.y) {
                            std::swap(*pP4, *--pP4_end);
                        } else {
                            ++pP4;
                        }
                    }
                    moth_size_t nThreads2 = moth_cdiv2(nThreads);
                    moth_omp_tasks(nThreads2 >= 2, [=]() {
                        moth_size_t nThreads4 = moth_fdiv2(nThreads2);
                        moth_sort_recursive(pP3, pP4, {p_cnr.x, p_cnr.y}, {p_max.x, p_max.y}, 1, nThreads4);
                    }, [=]() {
                        moth_size_t nThreads4 = moth_cdiv2(nThreads2);
                        moth_sort_recursive(pP4, pP_end, {p_cnr.x, p_min.y}, {p_max.x, p_cnr.y}, 3, nThreads4);
                    });
                });
                break;

            /* 4||1
             * 3--2 Orientation. */
            case 2:
                /* Separate right and left halves. */
                pP3 = pP1;
                for (moth_p2d* pP3_end{pP_end}; pP3 != pP3_end;) {
                    if (pP3->x < p_cnr.x) {
                        std::swap(*pP3, *--pP3_end);
                    } else {
                        ++pP3;
                    }
                }
                moth_omp_tasks(nThreads >= 2, [=]() {
                    /* Separate upper and lower quadrants of the right half
                     * and recursively process the quadrants. */
                    moth_p2d* pP2 = pP1;
                    for (moth_p2d* pP2_end{pP3}; pP2 != pP2_end;) {
                        if (pP2->y < p_cnr.y) {
                            std::swap(*pP2, *--pP2_end);
                        } else {
                            ++pP2;
                        }
                    }
                    moth_size_t nThreads2 = moth_fdiv2(nThreads);
                    moth_omp_tasks(nThreads2 >= 2, [=]() {
                        moth_size_t nThreads4 = moth_fdiv2(nThreads2);
                        moth_sort_recursive(pP1, pP2, {p_cnr.x, p_cnr.y}, {p_max.x, p_max.y}, 3, nThreads4);
                    }, [=]() {
                        moth_size_t nThreads4 = moth_cdiv2(nThreads2);
                        moth_sort_recursive(pP2, pP3, {p_cnr.x, p_min.y}, {p_max.x, p_cnr.y}, 2, nThreads4);
                    });
                }, [=]() {
                    /* Separate lower and upper quadrants of the left half
                     * and recursively process the quadrants. */
                    moth_p2d* pP4 = pP3;
                    for (moth_p2d* pP4_end{pP_end}; pP4 != pP4_end;) {
                        if (pP4->y > p_cnr.y) {
                            std::swap(*pP4, *--pP4_end);
                        } else {
                            ++pP4;
                        }
                    }
                    moth_size_t nThreads2 = moth_cdiv2(nThreads);
                    moth_omp_tasks(nThreads2 >= 2, [=]() {
                        moth_size_t nThreads4 = moth_fdiv2(nThreads2);
                        moth_sort_recursive(pP3, pP4, {p_min.x, p_min.y}, {p_cnr.x, p_cnr.y}, 2, nThreads4);
                    }, [=]() {
                        moth_size_t nThreads4 = moth_cdiv2(nThreads2);
                        moth_sort_recursive(pP4, pP_end, {p_min.x, p_cnr.y}, {p_cnr.x, p_max.y}, 0, nThreads4);
                    });
                });
                break;

            /* 2|-1
             * 3|-4 Orientation. */
            case 3:
                /* Separate upper and lower halves. */
                pP3 = pP1;
                for (moth_p2d* pP3_end{pP_end}; pP3 != pP3_end;) {
                    if (pP3->y < p_cnr.y) {
                        std::swap(*pP3, *--pP3_end);
                    } else {
                        ++pP3;
                    }
                }
                moth_omp_tasks(nThreads >= 2, [=]() {
                    /* Separate left and right quadrants of the upper half
                     * and recursively process the quadrants. */
                    moth_p2d* pP2 = pP1;
                    for (moth_p2d* pP2_end{pP3}; pP2 != pP2_end;) {
                        if (pP2->x < p_cnr.x) {
                            std::swap(*pP2, *--pP2_end);
                        } else {
                            ++pP2;
                        }
                    }
                    moth_size_t nThreads2 = moth_fdiv2(nThreads);
                    moth_omp_tasks(nThreads2 >= 2, [=]() {
                        moth_size_t nThreads4 = moth_fdiv2(nThreads2);
                        moth_sort_recursive(pP1, pP2, {p_cnr.x, p_cnr.y}, {p_max.x, p_max.y}, 2, nThreads4);
                    }, [=]() {
                        moth_size_t nThreads4 = moth_cdiv2(nThreads2);
                        moth_sort_recursive(pP2, pP3, {p_min.x, p_cnr.y}, {p_cnr.x, p_max.y}, 3, nThreads4);
                    });
                }, [=]() {
                    /* Separate right and left quadrants of the lower half
                     * and recursively process the quadrants. */
                    moth_p2d* pP4 = pP3;
                    for (moth_p2d* pP4_end{pP_end}; pP4 != pP4_end;) {
                        if (pP4->x > p_cnr.x) {
                            std::swap(*pP4, *--pP4_end);
                        } else {
                            ++pP4;
                        }
                    }
                    moth_size_t nThreads2 = moth_cdiv2(nThreads);
                    moth_omp_tasks(nThreads2 >= 2, [=]() {
                        moth_size_t nThreads4 = moth_fdiv2(nThreads2);
                        moth_sort_recursive(pP3, pP4, {p_min.x, p_min.y}, {p_cnr.x, p_cnr.y}, 3, nThreads4);
                    }, [=]() {
                        moth_size_t nThreads4 = moth_cdiv2(nThreads2);
                        moth_sort_recursive(pP4, pP_end, {p_cnr.x, p_min.y}, {p_max.x, p_cnr.y}, 1, nThreads4);
                    });
                });
                break;

            default:
                std::cerr << "Invalid orientation of the Hilbert sort." << std::endl;
                std::abort();
        }
    }
}

/**
 * 2D Hilbert sort of the given point set -- O(n*log(n)).
 */
MOTH_HOST MOTH_CORE
void moth_sort(moth_p2d* pP_beg, moth_p2d* pP_end)
{
    assert(pP_beg != nullptr);
    assert(pP_end != nullptr);
    if (pP_beg + 1 < pP_end) {
        /* Find the bounding box. */
        moth_p2d p_min{*pP_beg};
        moth_p2d p_max{*pP_beg};
        for (moth_p2d* pP_cur = pP_beg + 1;
                       pP_cur != pP_end; ++pP_cur) {
            p_min = moth_p2d::min(p_min, *pP_cur);
            p_max = moth_p2d::max(p_max, *pP_cur);
        }

        /* Perform the sort. */
#if MOTH_OMP
#pragma omp parallel
        {
#pragma omp single nowait
            {
                moth_sort_recursive(pP_beg, pP_end, p_min, p_max, 1,
                                    moth_omp_get_num_threads());
            }
        }
#else
        moth_sort_recursive(pP_beg, pP_end, p_min, p_max, 1, 1);
#endif
    }
#if 0
    std::string file_path("res/st-" + std::to_string(99999) + ".txt");
    std::ofstream file(file_path);
    for (moth_p2d* pP_cur = pP_beg; pP_cur != pP_end; ++pP_cur) {
        file << *pP_cur << std::endl;
    }
#endif
}

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * 3D Hilbert sort of the given point set -- O(n*log(n)).
 */
MOTH_HOST MOTH_CORE
void moth_sort(moth_p3d* pP_beg, moth_p3d* pP_end)
{
    std::cerr << "Nor implemented." << std::endl;
    std::abort();
}
