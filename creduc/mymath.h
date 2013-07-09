#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/types.h>
//
// Overview
//
template <class T> T *intep(int32_t &np, float32_t *x, float32_t *y, int32_t &nxp, T *xp);
template <class T> T *invert4x4(T *a);
template <class T> void Shipley_Inversion(T **M,int m);

//
// Definitions
//
template <class T> T *intep(int32_t &np, float32_t *x, float32_t *y, int32_t &nxp, T *xp){
  //
  int32_t np1 = np - 1;
  T *yp  = new T [nxp];
  T *lp1 = new T [np1];
  T *lp2 = new T [np1];
  T *fp1 = new T [np1];
  T *fp2 = new T [np1];
  //
  // Compute interpolation coefficients
  //
  for(int32_t ii=0;ii<=np-2;++ii) lp1[ii] = 1.0 / (x[ii] - x[ii+1]);
  for(int32_t ii=0;ii<=np-2;++ii) lp2[ii] = 1.0 / (x[ii+1] - x[ii]);
  //
  for(int32_t ii=1;ii<=np-2;++ii) fp1[ii] = (y[ii+1] - y[ii-1]) / (x[ii+1] - x[ii-1]);
  fp1[0] = (y[1] - y[0]) / (x[1] - x[0]);
  //
  fp2[np-2] = (y[np-1] - y[np-2]) / (x[np-1] - x[np-2]);
  for(int32_t ii=0;ii<=np-3;++ii) fp2[ii] = fp1[ii+1];
  //
  // output dependent part
  //
  int32_t i, i1;
  T xpi, xpi1, l1, l2;
  for(int32_t ii=1;ii<=np;++ii){
    i = ii - 1;
    i1 = i - 1;
    if(i < 0) i = 0;
    if(i1 < 0) i1 = 0;
    //
    for(int32_t jj = 0;jj<=nxp-1;++jj){
      if((xp[jj] <= x[i]) && (xp[jj] >= x[i1])){
	//
	xpi = xp[jj] - x[i1];
	xpi1 = xp[jj] - x[i];
	//
	l1 = xpi1 * lp1[i1];
	l1*= l1;
	//
	l2 = xpi * lp2[i1];
	l2*= l2;
	//
	yp[jj] = y[i1] * (1.0 - (2.0 * lp1[i1] * xpi)) * l1 + 
	  y[i] * (1.0 - (2.0 * lp2[i1] * xpi1)) * l2 + 
	  fp2[i1] * xpi1 * l2 + fp1[i1] * xpi * l1;
      }
    }
  }
  //
  // any point outside bounds? -> linear extrapolation
  //
  T a0 = fp1[0]; //(y[0] - y[1]) / (x[0] - x[1]);
  T b0 = y[0] - a0 * x[0];
  //
  T a1 = fp2[np-2]; // (y[np-2] - y[np-1]) / (x[np-2] - x[np-1]);
  T b1 = y[np-1] - a1 * x[np-1];
  //
  for(int32_t ii = 0;ii<=nxp-1;++ii){
    if(xp[ii] < x[0]) yp[ii] = a0 * xp[ii] + b0;
    if(xp[ii] > x[np-1]) yp[ii] = a1 * xp[ii] + b1;
  }
  //
  // Clean-up
  //
  delete [] lp1;
  delete [] lp2;
  delete [] fp1;
  delete [] fp2;
  //
  return yp;
}
/*
  Invert 4x4 matrix. Adapted from the fortran code NICOLE.
 */
template <class T> T *invert4x4(T *a){
  // 
  // Init result
  //
  T *b = new T [16];
  //
  T fabsmx = fabs(a[0]);
  memset(&b[0], 0, 16 * sizeof(T));
  //
  //for(uint8_t ii=1;ii<16;ii++) fabsmx = max(fabs(a[ii]),fabsmx);
  //fabsmx = 1.0 / fabsmx;
  //
  // operations
  //
  b[0] = a[5] * a[10] * a[15] + a[9] * a[14] * a[7] 
    + a[13] * a[6] * a[11] - a[5] * a[14] * a[11] 
    - a[9] * a[6] * a[15] - a[13] * a[10] * a[7];
  b[1] = a[9] * a[2] * a[15] + a[13] * a[10] * a[3] 
    + a[1] * a[14] * a[11] - a[9] * a[14] * a[3] 
    - a[13] * a[2] * a[11] - a[1] * a[10] * a[15];
  b[2] = a[13] * a[2] * a[7] + a[1] * a[6] * a[15]
    + a[5] * a[14] * a[3] - a[13] * a[6] * a[3]
    - a[1] * a[14] * a[7] - a[5] * a[2] * a[15];
  b[3] = a[1] * a[10] * a[7] + a[5] * a[2] * a[11] 
       + a[9] * a[6] * a[3] - a[1] * a[6] * a[11] 
    - a[5] * a[10] * a[3] - a[9] * a[2] * a[7];
  b[4] = a[6] * a[15] * a[8] + a[10] * a[7] * a[12] 
       + a[14] * a[11] * a[4] - a[6] * a[11] * a[12] 
    - a[10] * a[15] * a[4] - a[14] * a[7] * a[8];
  b[5] = a[10] * a[15] * a[0] + a[14] * a[3] * a[8] 
       + a[2] * a[11] * a[12] - a[10] * a[3] * a[12] 
    - a[14] * a[11] * a[0] - a[2] * a[15] * a[8];
  b[6] = a[14] * a[7] * a[0] + a[2] * a[15] * a[4] 
       + a[6] * a[3] * a[12] - a[14] * a[3] * a[4] 
    - a[2] * a[7] * a[12] - a[6] * a[15] * a[0];
  b[7] = a[2] * a[7] * a[8] + a[6] * a[11] * a[0] 
       + a[10] * a[3] * a[4] - a[2] * a[11] * a[4] 
    - a[6] * a[3] * a[8] - a[10] * a[7] * a[0];
  b[8] = a[7] * a[8] * a[13] + a[11] * a[12] * a[5] 
       + a[15] * a[4] * a[9] - a[7] * a[12] * a[9] 
    - a[11] * a[4] * a[13] - a[15] * a[8] * a[5];
  b[9] = a[11] * a[0] * a[13] + a[15] * a[8] * a[1] 
       + a[3] * a[12] * a[9] - a[11] * a[12] * a[1] 
    - a[15] * a[0] * a[9] - a[3] * a[8] * a[13];
  b[10] = a[15] * a[0] * a[5] + a[3] * a[4] * a[13] 
       + a[7] * a[12] * a[1] - a[15] * a[4] * a[1] 
    - a[3] * a[12] * a[5] - a[7] * a[0] * a[13];
  b[11] = a[3] * a[8] * a[5] + a[7] * a[0] * a[9] 
       + a[11] * a[4] * a[1] - a[3] * a[4] * a[9] 
    - a[7] * a[8] * a[1] - a[11] * a[0] * a[5];
  b[12] = a[4] * a[13] * a[10] + a[8] * a[5] * a[14] 
       + a[12] * a[9] * a[6] - a[4] * a[9] * a[14] 
    - a[8] * a[13] * a[6] - a[12] * a[5] * a[10];
  b[13] = a[8] * a[13] * a[2] + a[12] * a[1] * a[10] 
       + a[0] * a[9] * a[14] - a[8] * a[1] * a[14] 
    - a[12] * a[9] * a[2] - a[0] * a[13] * a[10];
  b[14] = a[12] * a[5] * a[2] + a[0] * a[13] * a[6] 
       + a[4] * a[1] * a[14] - a[12] * a[1] * a[6] 
    - a[0] * a[5] * a[14] - a[4] * a[13] * a[2];
  b[15] = a[0] * a[5] * a[10] + a[4] * a[9] * a[2] 
       + a[8] * a[1] * a[6] - a[0] * a[9] * a[6] 
    - a[4] * a[1] * a[10] - a[8] * a[5] * a[2];
  
  //
  // Determinant
  //
  T idet = 1.0 / (a[0] * b[0] + a[4] * b[1] + a[8] * b[2] + a[12] * b[3]);
  //
  for(short ii=0;ii<16;++ii) b[ii] *= idet;
  //
  return b;
}

template <class T> void Shipley_Inversion(T **MI,int m){
  //
  for(int k=0;k<m;k++){
    MI[k][k]=-1.0/MI[k][k];         // the pivot element 
    for(int i=0;i<m;++i) if(i!=k) MI[i][k]*=MI[k][k];//the pivot column 
    for(int i=0;i<m;++i)           //elements not in a pivot row or column
      if(i!=k)
        for(int j=0;j<m;++j)
          if(j!=k)
            MI[i][j]+=MI[i][k]*MI[k][j];
    for(int i=0;i<m;++i)           //elements in a pivot row
      if(i!=k)
        MI[k][i]*=MI[k][k];
  }
  for(int i=0;i<m;++i)
    for(int j=0;j<m;++j) MI[i][j]=-MI[i][j];
  return;
}
