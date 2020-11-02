#include <cmath>
#include <cstring>
#include <cstdio>
#include <omp.h>
#include <iostream>

//#include "Arrays.hpp"

using namespace std;

// ********************************************************************************************************** //

// ---------------------------------------------------------------------- //
// Nearest neighbor 2D interpolation. It assumes that the input image is
// given in a regular grid, but the output does not need to be given in a
// regular grid. That way we can use rotation.
//
// Coded by J. de la Cruz Rodriguez (ISP-SU 2019)
// ---------------------------------------------------------------------- //

void nearest2D(int const ny, int const nx, double* __restrict__ y, double* __restrict__ x,
	       double*  d_in, int const ny1, int const nx1, double*  yy_in,
	       double*  xx_in,  double*  res_in, int const nthreads, double const missing)
{

  // --- map data into 2D array ---/

  //mem::Array<double,2>  d(d_in,    ny, nx);
  //mem::Array<double,2> xx(xx_in,  ny1, nx1);
  //mem::Array<double,2> yy(yy_in,  ny1, nx1);
  //mem::Array<double,2> res(res_in,ny1, nx1);

  // --- Need to pre-define these types so Clang swallows it --- //
  using ft = double[ny][nx];
  using ft1= double[ny1][nx1];
  
  ft  &d   = *reinterpret_cast<ft*>(d_in);
  ft1 &xx  = *reinterpret_cast<ft1*>(xx_in);
  ft1 &yy  = *reinterpret_cast<ft1*>(yy_in);
  ft1 &res = *reinterpret_cast<ft1*>(res_in);

  
  // --- Now find the closest y and x values for each interpolated value ---//

  int jj=0, ii=0, iy=0, ix=0, tt=0;
  double dy0=0.0, dy1=0.0, dx1=0.0, dx0=0.0, inty=0.0, intx=0.0;

#pragma omp parallel default(shared) firstprivate(jj,ii,ix,iy,tt,dy0,dy1,dx1,dx0) num_threads(nthreads)  
  {
#pragma omp for schedule(dynamic,10) collapse(2)
    for( jj=0; jj<ny1; ++jj){
      for( ii=0; ii<nx1; ++ii){
	//inty = std::min<double>(std::max<double>(yy[jj][ii],y[0]*1.0000000001),y[ny-1]*0.999999999999);
	//intx = std::min<double>(std::max<double>(xx[jj][ii],x[0]*1.0000000001),x[nx-1]*0.999999999999);
	
	intx = xx[jj][ii], inty = yy[jj][ii];
	if((intx < x[0]) || (intx > x[nx-1]) || (inty < y[0]) || (inty > y[ny-1])) res[jj][ii] = missing;
	else{
	  
	  // --- Search y --- //
	  iy = 0;
	  
	  
	  for( tt=1; tt<ny; ++tt){
	    if((inty > y[tt-1]) && (inty <= y[tt])){
	      dy0 = std::abs(inty-y[tt-1]);
	      dy1 = std::abs(inty-y[tt]);
	      
	      if(dy0 < dy1) iy = tt-1;
	      else iy = tt;
	      
	      break;
	    } // if
	  } //tt
	  
	  // --- now search xx --- //
	  
	  ix = 0;
	  for( tt=1; tt<nx; ++tt){
	    if((intx > x[tt-1]) && (intx <= x[tt])){
	      dx0 = std::abs(intx-x[tt-1]);
	      dx1 = std::abs(intx-x[tt]);
	      
	      if(dx0 < dx1) ix = tt-1;
	      else ix = tt;
	      
	      break;
	    } // if
	  } // tt
	  
	  
	  // --- Assign result --- //
	  
	  res[jj][ii] = d[iy][ix];
	  
	}
      } // ii 
    } // jj
  } // pragma
}

// ********************************************************************************************************** //

// ---------------------------------------------------------------------- //
// 2D bilinear interpolation. It assumes that the input image is
// given in a regular grid, but the output does not need to be given in a
// regular grid. That way we can use rotation.
//
// Coded by J. de la Cruz Rodriguez (ISP-SU 2019)
// ---------------------------------------------------------------------- //

void bilint2D(int const ny, int const nx, double* __restrict__ y, double* __restrict__ x,
	      double*  d_in, int const ny1, int const nx1, double*  yy_in,
	      double*  xx_in, double* res_in, int const nthreads, double const missing)
{
  
  double* __restrict__ c_in = new double [4*(ny-1)*(nx-1)]();

  // --- Need to pre-define these types so Clang++ swallows the "reinterpret_cast" --- //
  
  using ft = double[ny][nx];
  using ft1= double[ny1][nx1];
  using ft2= double[ny-1][nx-1][4];
  
  ft  &d   = *reinterpret_cast<ft*>(d_in);

  ft1 &xx  = *reinterpret_cast<ft1*>(xx_in);
  ft1 &yy  = *reinterpret_cast<ft1*>(yy_in);
  ft1 &res = *reinterpret_cast<ft1*>(res_in);
  
  ft2 &c   = *reinterpret_cast<ft2*>(c_in);

  
  
  // --- compute interpolation coefficients --- //

  int const nxx = nx-1;
  int const nyy = ny-1;

  const double* __restrict__ cc = NULL;
  double dydx=0, x1=0,x2=0,dy=0,y1=0,y2=0,dx=0,ixx=0,iyy=0;
  int ii=0, jj=0,idx=0,idy=0,tt=0;
  
#pragma omp parallel default(shared) firstprivate(jj,ii,dy,y1,y2,dx,x1,x2,dydx) num_threads(nthreads)  
  {
#pragma omp for schedule(static) collapse(2)
    for( jj=0; jj<nyy; ++jj){
      for( ii=0; ii<nxx; ++ii){
	
	dy = y[jj] - y[jj+1];
	dx = x[ii] - x[ii+1];
	
	dydx = dx*dy;

	y1 = y[jj]/dydx;
	y2 = y[jj+1]/dydx;
	x1 = x[ii]/dydx;
	x2 = x[ii+1]/dydx;
	
	c[jj][ii][0] =( d[jj][ii]*x2*y2 - d[jj+1][ii]*x2*y1 - d[jj][ii+1]*x1*y2 + d[jj+1][ii+1]*x1*y1);
	c[jj][ii][1] =(-d[jj][ii]*y2    + d[jj+1][ii]*y1    + d[jj][ii+1]*y2    - d[jj+1][ii+1]*y1   );
	c[jj][ii][2] =(-d[jj][ii]*x2    + d[jj+1][ii]*x2    + d[jj][ii+1]*x1    - d[jj+1][ii+1]*x1   );
	c[jj][ii][3] =( d[jj][ii]       - d[jj+1][ii]       - d[jj][ii+1]       + d[jj+1][ii+1]      );
      }//ii
    }//jj
  }// pragma

  // --- Do interpolation --- //
  
#pragma omp parallel default(shared) firstprivate(jj,ii,ixx,iyy,idx,idy,tt,cc) num_threads(nthreads)  
  {
#pragma omp for schedule(static) collapse(2)
    for( jj=0; jj<ny1; ++jj){
      for( ii=0; ii<nx1; ++ii){

	// --- If the pixel is out of bounds, force coefficients at the edge --- //
	ixx = xx[jj][ii], iyy = yy[jj][ii];
	
	if((ixx < x[0]) || (ixx > x[nx-1]) || (iyy < y[0]) || (iyy > y[ny-1])) res[jj][ii] = missing;
	else{ // pixel inside data bounds
	  
	  //ixx = std::min<double>(std::max<double>(xx[jj][ii],x[0]*1.0000000001),x[nx-1]*0.999999999999);
	  //iyy = std::min<double>(std::max<double>(yy[jj][ii],y[0]*1.0000000001),y[ny-1]*0.999999999999);
	  idx = 0, idy = 0;
	  
	  
	  // --- bracket coordinates --- //
	  
	  for( tt=0; tt<nyy; ++tt){
	    if((iyy >= y[tt]) && (iyy < y[tt+1])){
	      idy = tt;
	      break;
	    } // if
	  }
	  
	  for( tt=0; tt<nxx; ++tt){
	    if((ixx >= x[tt]) && (ixx < x[tt+1])){
	      idx = tt;
	      break;
	    } // if
	  }
	  
	  
	  // --- Apply interpolation coefficients --- //
	  
	  cc = static_cast<const double* __restrict__>(&c[idy][idx][0]);
	  res[jj][ii] = cc[0] + cc[1]*ixx + cc[2]*iyy + cc[3]*iyy*ixx;

	} // else
      } // ii
    } // jj
  }// pragma

  delete [] c_in;
  
}

// ********************************************************************************************************** //
