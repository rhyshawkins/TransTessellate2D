
#include <math.h>

#include "solver2x2.h"


int
Solve2x2(double a11,
         double a12,
         double a21,
         double a22,
         double b1,
         double b2,
         double *x1,
         double *x2)
{
  double ur11;
  double cr21;
  double ur12;
  double cr22;

  int rswap;
  int xswap;

  double ur11r;
  double lr21;
  double ur22;

  double br1, br2;
  
  double xr2;
  double xr1;
  
  //
  // Pivot to place largest value at 11
  //
  if (fabs(a11) > fabs(a12)) {
    if (fabs(a11) > fabs(a21)) {
      if (fabs(a11) > fabs(a22)) {
        // 0
        ur11 = a11;
        cr21 = a21;
        ur12 = a12;
        cr22 = a22;

        rswap = 0;
        xswap = 0;
        
      } else {
        // 3
        ur11 = a22;
        cr21 = a12;
        ur12 = a21;
        cr22 = a11;

        rswap = 1;
        xswap = 1;
        
      }
    } else {
      if (fabs(a21) > fabs(a22)) {
        // 2
        ur11 = a21;
        cr21 = a11;
        ur12 = a22;
        cr22 = a12;

        rswap = 1;
        xswap = 0;
        
      } else {
        // 3
        ur11 = a22;
        cr21 = a12;
        ur12 = a21;
        cr22 = a11;

        rswap = 1;
        xswap = 1;
      }
    }
  } else {
    if (fabs(a12) > fabs(a21)) {
      if (fabs(a12) > fabs(a22)) {
        // 1
        ur11 = a12;
        cr21 = a22;
        ur12 = a11;
        cr22 = a21;

        rswap = 0;
        xswap = 1;
        
      } else {
        // 3
        ur11 = a22;
        cr21 = a12;
        ur12 = a21;
        cr22 = a11;

        rswap = 1;
        xswap = 1;
      }
    } else {
      if (fabs(a21) > fabs(a22)) {
        // 2 
        ur11 = a21;
        cr21 = a11;
        ur12 = a22;
        cr22 = a12;

        rswap = 1;
        xswap = 0;
      } else {
        // 3
        ur11 = a22;
        cr21 = a12;
        ur12 = a21;
        cr22 = a11;

        rswap = 1;
        xswap = 1;
      }
    }
  }

  ur11r = 1.0/ur11;
  lr21 = ur11r*cr21;
  ur22 = cr22 - ur12*lr21;

  if (rswap) {
    br1 = b2;
    br2 = b1;
  } else {
    br1 = b1;
    br2 = b2;
  }

  br2 -= lr21*br1;

  xr2 = br2/ur22;
  xr1 = br1*ur11r - xr2*(ur11r*ur12);
  
  if (xswap) {
    *x1 = xr2;
    *x2 = xr1;
  } else {
    *x1 = xr1;
    *x2 = xr2;
  }
    
  return 1;
}
