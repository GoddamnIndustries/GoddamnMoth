#include "libGeometry2D/src/GeomPoint.hh"

// ------------------------------------------------------------------------------------ //
// ------------------------------------------------------------------------------------ //

/**
 * 2D Hilbert sort of the given point set -- O(n*log(n)).
 */
MOTH_HOST MOTH_CORE
extern void moth_sort(moth_p2d* pP_beg, moth_p2d* pP_end);

/**
 * 3D Hilbert sort of the given point set -- O(n*log(n)).
 */
MOTH_HOST MOTH_CORE
extern void moth_sort(moth_p3d* pP_beg, moth_p3d* pP_end);
