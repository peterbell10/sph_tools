#ifndef _H_PROJ_H_
#define _H_PROJ_H_

#include <stddef.h>

#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))

#define H_THRESHOLD_FOR_DIRECT_BINNING 0.1
#define H_THRESHOLD_FOR_REFINEMENT     1.0


#ifdef __cplusplus
extern "C" {
#endif

void h_proj_3d_core(const double * const r,
                    const double * const phi,
		    const double * const h,
		    const size_t L,
		    const size_t N,
		    double * const A);

double kernel(const double x);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif
