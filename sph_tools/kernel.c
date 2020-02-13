/* SPH Smoothing Kernel (whose?) */

#include <math.h>

#include "h_proj.h"

double kernel(const double x) {

/* x == r/h
 * Compact over 0 to 2*h. This is different from Gadget's, which is
 * compact over 0 to h. */

 double W=0.0;
 const double wnorm = 0.25/M_PI;

 if ( (x > 0.0) && (x < 2.0) ) {
  if (x <= 1.0) {
  /* W=wnorm*(4.0-6.0*x*x+3.0*x*x*x); */
  W = wnorm*(4.0-x*x*(6.0-3.0*x));
  } else {
   W = wnorm*(2.0-x)*(2.0-x)*(2.0-x);
  }
 }
 return W;
}
