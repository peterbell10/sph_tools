/* h_proj_3d_core: Distributes the values of phi at the locations r onto the LxLxL grid A.
 * sum(sum(sum(A))) == sum(phi)
 *
 * ARGUMENTS
 *  Input, not modified
 *   r   [3,N] array of positions, on the span [0,L]
 *   phi N-vector of values to distribute
 *   h   N-vector of smoothing radii, with the same length unit as r
 *   L   The number of nodes per side of the output array
 *   N   The number of particles (see r, phi, h)
 *  Output, initialised here but allocated elsewhere
 *   A           [LxLxL] array.  The mesh containing the distributed values.
 *
 * NOTES
 * - r is on the span 0:L.  Whatever r is originally, it is put on the
 *   span 0:1 in h_proj_3d, and { put on the span 0:L
 *
 * - This method uses the scatter interpretation of SPH.  Each particle's z
 *   component is distributed to the cells within 2h of it's location. 
 *
 * AUTHOR: Eric Tittley
 *
 * TODO
 *  Doesn't deal with repeating boundary conditions.  Assumes they are not.
 */

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

#include "h_proj.h"

static int ceilDiv(int q, int d) {
 return (q + d - 1)  / d;
}

void processParticlesInSlice(const double * const r,
                             const double * const phi,
		             const double * const h,
		             const size_t L,
		             const size_t N,
		             double * const A,
			     const int iz_lo,
			     const int iz_hi);

void h_proj_3d_core(const double * const r,
                    const double * const phi,
		    const double * const h,
		    const size_t L,
		    const size_t N,
		    double * const A) {
		    
 /* Initialise the output array */
 size_t i; /* Particle indexing */
#pragma omp parallel for
 for(i=0;i<(L*L*L);++i) {
  A[i]=0.0;
 }

 #pragma omp parallel
 {
  const int nThreads = omp_get_num_threads();
  const int myId = omp_get_thread_num();
  const int cellsPerChunk = ceilDiv((int)L, nThreads);

  const int iz_lo = myId * cellsPerChunk;
  const int iz_hi = MIN((int)L, (myId + 1) * cellsPerChunk);
  processParticlesInSlice(r,phi,h,L,N,A,iz_lo,iz_hi);
 }
}


void processParticlesInSlice(const double * const r,
                             const double * const phi,
		             const double * const h,
		             const size_t L,
		             const size_t N,
		             double * const A,
			     const int iz_lo,
			     const int iz_hi)
{
 /* Loop over all particles */
 size_t i;
 double z_lo = (double)iz_lo;
 double z_hi = (double)iz_hi;
 for(i=0;i<N;++i) {
  size_t j = 3*i;
  double z=r[j+2];
  double hi=h[i];
  double twohi=2.0*hi;
  if( (z + twohi >= z_lo) && (z - twohi < z_hi) ) {
   double x=r[j+0];
   double y=r[j+1];

   /* See the Appendix for an explanation of the following line */	     
   double phih3=phi[i]/(hi*hi*hi);
   if(hi<=H_THRESHOLD_FOR_DIRECT_BINNING) {
    size_t ix=(size_t)floor(x);
    size_t iy=(size_t)floor(y);
    size_t iz=(size_t)floor(z);
    size_t index=ix+iy*L+iz*L*L;
    if (iz >= iz_lo && iz < iz_hi) {
     A[index]=A[index] + phi[i];
    }
   } else if (hi<=H_THRESHOLD_FOR_REFINEMENT) {
    /* ****** Refine cells ****** */

    /* The limits of cells the particle contributes to.
     * These can be outside the range [0,L-1]. We'll deal with that case when we assign to cells.*/
    int ix_min=(int)floor(x-twohi);
    int iy_min=(int)floor(y-twohi);
    int iz_min=(int)floor(z-twohi);
    if(iz_min<iz_lo) iz_min=iz_lo;

    int ix_max=(int)floor(x+twohi);
    int iy_max=(int)floor(y+twohi);
    int iz_max=(int)floor(z+twohi);
    if(iz_max>iz_hi-1) iz_max=iz_hi-1;

    /* Size of the region the particle spreads itself over */
    size_t x_cell_range=(size_t)(ix_max-ix_min+1);
    size_t y_cell_range=(size_t)(iy_max-iy_min+1);
    size_t z_cell_range=(size_t)(iz_max-iz_min+1);

    /* Needs to be cubic */
    size_t cell_range=MAX(x_cell_range,y_cell_range);
    cell_range=MAX(cell_range,z_cell_range);
    
    /* Deterimine the size of the refinement */
    double h_ref=hi;
    size_t L_ref=cell_range;
    while(h_ref<2.0) {
     h_ref=h_ref*2.0;
     L_ref=L_ref*2;
    }
    size_t refinementFactor=L_ref/cell_range;
    double * A_ref = (double *)calloc(L_ref*L_ref*L_ref,sizeof(double));
    if(A_ref==NULL) {
     printf("ERROR: %s: %i: Unable to allocate memory\n",__FILE__,__LINE__);
     return;
    }
    double r_ref[3];
    /* r is in the range [0 L)
     * Need to put it in [0 L_ref) */
    r_ref[0]=(r[j+0]-(double)ix_min)*((double)refinementFactor);
    r_ref[1]=(r[j+1]-(double)iy_min)*((double)refinementFactor);
    r_ref[2]=(r[j+2]-(double)iz_min)*((double)refinementFactor);
    /* Iteratively call the routine that does the work */
    processParticlesInSlice(r_ref,&(phi[i]),&h_ref,L_ref,1,A_ref,0,L_ref);
    
    /* Assign the values in the refined grid to the appropriate base grid cells */
    size_t ix_ref,iy_ref,iz_ref;
    for (iz_ref=0;iz_ref<L_ref;++iz_ref) {
     int iz=iz_min+(int)(iz_ref/refinementFactor); /* Integer division */
     if(iz<0) iz=iz+(int)L;       /* Assume repeating boundary conditions */
     if(iz>=(int)L) iz=iz-(int)L; /* Assume repeating boundary conditions */
     size_t iz_L    = (size_t)iz*L*L;
     size_t iz_ref_L=     iz_ref*L_ref*L_ref;
     for (iy_ref=0;iy_ref<L_ref;++iy_ref) {
      int iy=iy_min+(int)(iy_ref/refinementFactor); /* Integer division */
      if(iy<0) iy=iy+(int)L;       /* Assume repeating boundary conditions */
      if(iy>=(int)L) iy=iy-(int)L; /* Assume repeating boundary conditions */
      size_t iy_L    = iz_L     + (size_t)iy*L;
      size_t iy_ref_L= iz_ref_L +   iy_ref  *L_ref;
      for (ix_ref=0;ix_ref<L_ref;++ix_ref) {
       int ix=ix_min+(int)(ix_ref/refinementFactor); /* Integer division */
       if(ix<0) ix=ix+(int)L;       /* Assume repeating boundary conditions */
       if(ix>=(int)L) ix=ix-(int)L; /* Assume repeating boundary conditions */
       size_t index     =(size_t)ix + iy_L;
       size_t index_ref =    ix_ref + iy_ref_L;
       A[index]=A[index]+A_ref[index_ref];
      }
     }
    }
    /* Clean up */
    free(A_ref);
   } else {
    int ix_min=(int)floor(x-twohi);
    if(ix_min<0) ix_min=0;
    int iy_min=(int)floor(y-twohi);
    if(iy_min<0) iy_min=0;
    int iz_min=(int)floor(z-twohi);
    if(iz_min<iz_lo) iz_min=iz_lo;

    int ix_max=(int)floor(x+twohi);
    if(ix_max>(int)(L-1)) ix_max=(int)(L-1);
    int iy_max=(int)floor(y+twohi);
    if(iy_max>(int)(L-1)) iy_max=(int)(L-1);
    int iz_max=(int)floor(z+twohi);
    if(iz_max>iz_hi-1) iz_max=iz_hi-1;

    double xPhalf=x-0.5;
    double yPhalf=y-0.5;
    double zPhalf=z-0.5;

    size_t iz,iy,ix;
    for(iz=(size_t)iz_min;iz<=(size_t)iz_max;++iz) {
     double dz2 = ((double)iz-zPhalf)*((double)iz-zPhalf);
     size_t iz_L = iz*L*L;
     for(iy=(size_t)iy_min;iy<=(size_t)iy_max;++iy) {
      double dy2 = dz2 + ((double)iy-yPhalf)*((double)iy-yPhalf);
      size_t iy_L = iy*L + iz_L; 
      for(ix=(size_t)ix_min;ix<=(size_t)ix_max;++ix) {
       double dist=sqrt( dy2 + ((double)ix-xPhalf)*((double)ix-xPhalf) )/hi;
       size_t index=ix+iy_L;
       A[index]=A[index]+kernel(dist)*phih3;
      }
     }
    }
   }
  } /* end if particle in slice */
 }
}

/*
 * ----------------------------------------
 * ------- APPENDIX -----------------------
 * ----------------------------------------
 *  A note about the normalisation of the weighted summation
 *  SPH assumes VolumeIntegrate( dr' W(r-r',h) ) == 1
 *  Let's take the kernel() to be w().
 *  LinearIntegrate( 4*pi*x^2*w(x) ) = 1 for x = 0 to 2, x == r/h
 *  But, LinearIntegrate( 4*pi*r^2*w(rh/h)) = h^3 for r = 0 to 2h
 *  So, we must normalise the weighting factor by h^3.
 *  Hence, W(x)==w(x)/h^3
 *
 *  In SPH:
 *   <A(r)> = sum A(r_i) * W(dist/h) / n_i
 *  Translating to this routine,
 *   <A(r)> = sum A(r_i) * w(dist/h) / (h_i^3 * n_i)
 *  In the limit of Nsph particles within 2*h_i,
 *   n_i = Nsph / ( 4/3 pi (2*h_i)^3 )
 *  So
 *   h_i^3 * n_i = 3 Nsph / (32 pi)
 *  Which implies the normalisation factor should be
 *   32 pi / ( 3 Nsph )
 *  But that does not work at all
 *  WHAT'S WRONG?
 *  Answer: The SPH formula:
 *   <A(r)> = sum A(r_i) * W(dist/h) / n_i
 *  is for finding the _average_ value of A at r.
 *  THAT IS NOT WHAT THIS ROUTINE IS FOR.
 *  This routine distrubutes the values of A (phi) onto a grid, such that 
 *   sum(sum(sum(A))) == sum(phi)
 *  Hence, no averaging and no /n(i).
 */
