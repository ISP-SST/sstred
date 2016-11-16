#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <sys/types.h>
#include <string.h>

template <class T> double *cent_der(int n, T *x, T *y);
template <class T> double *ccbezier2(int n, T *x, T *y, int np, double *xp);
template <class T> double *ccbezier3(int n, T *x, T *y, int np, double *xp);


template <class T> double *cent_der(int n, T *x, T *y){
  double *yprime = new double [n];
  memset(yprime, 0, n*sizeof(double));
  
  double dx = (x[1] - x[0]);
  double der = (y[1] - y[0]) / dx;
  yprime[0] = der;

  for (int k = 1;k<n-1;k++){
    double dx1 = x[k+1] - x[k];
    double der1 = (y[k+1] - y[k]) / dx1;


    if((der * der1) > 0.0){
      double lambda = (1.0 + dx1 / (dx1 + dx)) / 3.0;
      yprime[k] = (der / (lambda * der1 + (1.0 - lambda) * der)) * der1;
    }
    der = der1;
    dx = dx1;
    
  }
  
  yprime[n-1] = der;



  return yprime;
}


template <class T> double *ccbezier2(int n, T *xin, T *yin, int np, double *xp){
  //
  // Reverse arrays?
  //
  T *x, *y;
  short allocated = 0;
  if((xin[1]-xin[0]) > 0.0){
    x = (T*)xin;
    y = (T*)yin;
  }else{
    x = new T [n];
    y = new T [n];
    for(int k=0;k<n;k++){
      x[k] = xin[n-k-1];
      y[k] = yin[n-k-1];
     }
    allocated = 1;
  }


  //
  // Centered derivatives
  //
  double *yp = cent_der<T>(n, x, y);
  double *res = new double [np];


  //
  // Loop intervals
  //
  for(int k=0;k<n-1;k++){
    double dx = x[k+1] - x[k];
    double cntrl = 0.5 * (y[k] + 0.5 * dx * yp[k] + y[k+1] - dx * 0.5 * yp[k+1]);
    
    for(int j=0;j<np;j++){
      if((xp[j] >= x[k]) && (xp[j] < x[k+1])){
	double u = (xp[j] - x[k]) / dx;
	double u1 = 1.0 - u;
	res[j] =  y[k] * u1 * u1 + y[k+1] * u*u + 2.0 * cntrl * u * u1;
      }
    }
  }
  
 //
  // Out-of bound values, linear extrapolation
  //
  double aa0 = (y[1]-y[0])/(x[1]-x[0]);
  double bb0 = y[0]-aa0*x[0];
  double aa1 = (y[n-1]-y[n-2])/(x[n-1]-x[n-2]);
  double bb1 = y[n-1]-aa1*x[n-1];
  
  for(int k=0;k<np;k++){
    //    if(xp[k] < x[0]) res[k] = y[0];
    //    if(xp[k] >= x[n-1]) res[k] = y[n-1];
    if(xp[k] < x[0]) res[k] = xp[k]*aa0+bb0;
    if(xp[k] >= x[n-1]) res[k] = xp[k]*aa1+bb1;
  }

  if(allocated){
    delete [] x;
    delete [] y;
  }

  delete [] yp;
  return res;
}
template <class T> double *ccbezier3(int n, T *xin, T *yin, int np, double *xp){
  //
  // Reverse arrays?
  //
  T *x, *y;
  short allocated = 0;
  if((xin[1]-xin[0]) > 0.0){
    x = (T*)xin;
    y = (T*)yin;
  }else{
    x = new T [n];
    y = new T [n];
    for(int k=0;k<n;k++){
      x[k] = xin[n-k-1];
      y[k] = yin[n-k-1];
     }
    allocated = 1;
  }


  //
  // Centered derivatives
  //
  double *yp = cent_der<T>(n, x, y);
  double *res = new double [np];


  //
  // Loop intervals
  //
  for(int k=0;k<n-1;k++){
    double dx = x[k+1] - x[k];
    double cntrl1 = y[k] +  dx * yp[k] /3.0;
    double cntrl2 = y[k+1] - dx * yp[k+1]/3.0;
    for(int j=0;j<np;j++){
      if((xp[j] >= x[k]) && (xp[j] < x[k+1])){
	double u = (xp[j] - x[k]) / dx;
	double u1 = 1.0 - u;
	res[j] =  y[k] * u1 * u1 * u1 + y[k+1] * u*u*u + 3.0 * cntrl1 * u * u1*u1 + 3.0 * cntrl2 * u*u*u1;
      }
    }
  }
  
  //
  // Out-of bound values, linear extrapolation
  //
  double aa0 = (y[1]-y[0])/(x[1]-x[0]);
  double bb0 = y[0]-aa0*x[0];
  double aa1 = (y[n-1]-y[n-2])/(x[n-1]-x[n-2]);
  double bb1 = y[n-1]-aa1*x[n-1];
  
  for(int k=0;k<np;k++){
    //    if(xp[k] < x[0]) res[k] = y[0];
    //    if(xp[k] >= x[n-1]) res[k] = y[n-1];
    if(xp[k] < x[0]) res[k] = xp[k]*aa0+bb0;
    if(xp[k] >= x[n-1]) res[k] = xp[k]*aa1+bb1;
  }

  if(allocated){
    delete [] x;
    delete [] y;
  }

  delete [] yp;
  return res;
}
