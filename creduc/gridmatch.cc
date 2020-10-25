/*
  Function ana_stretch taken from LMSAL's package and conveted to operate only with DOUBLE variables.
  All IDL DLM module variables have been removed

  Some bugfixes:
  1) in ana_gridmatch, gx and gy are expected as int32 but dsgridnest.pro pases floats...
  2) all operations should be done as float64 otherwise we have found situations where NaNs appear.

  Ported by J. de la Cruz Rodriguez (ISP-SU, 2019)

 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#undef ABS
#define ABS(x) ((x) >= 0? (x): -(x))
/* lesser of two arguments */
#undef MIN
#define MIN(x,y) ((x) < (y)? (x): (y))
/* greater of two arguments */
#undef MAX
#define MAX(x,y) ((x) > (y)? (x): (y))





void ana_stretch_full_matrix(const int ny, const int nx, int npy, int npx, const double *gr,  double *out_x, double *out_y)
{
   double *bb;

  int	iq, type, n, m, nxg, nyg, jy, j1, j2, nm1, nm2, mm1, mm2;
  int	nxgm, nygm, ix, iy, i1, i2, i3, i4, j3, j4, jx;
  double	xd, yd, xinc, yinc, y, x, xs, dy, dy1;
  double	dx, dx1, dx0, dx2, dx3, dx4, fn, fm, xq;
  double	w1, w2, w3, w4, xl, yl, c1, c2, c3, c4, b1, b2, b3, b4;
  const double 	*jpbase, *jbase, *xgbase;


  n = nx;//argv[0]->value.arr->dim[0];
  m = ny;//argv[0]->value.arr->dim[1];
  nm1 = n - 1;
  mm1 = m - 1;
  nm2 = n -2;
  mm2 = m -2;
  fn = n;
  fm = m;

  xgbase = gr;//(double *) arg_struct[1].value.arr->data;
  nxg = npx;//arg_struct[1].value.arr->dim[1];
  nyg = npy;//arg_struct[1].value.arr->dim[2];
  nxgm = nxg - 1;
  nygm = nyg - 1;
  /* linearly interpolate the displacement grid values over array */
  /* similar to regrid3 in inner part */
  xd = (double) n/nxg; // number of pixels in the image per grid cell
  xinc = 1.0/xd; // maps the size of a pixel of the image in the matrix
  xs = xinc + (xd - 1.0)/(2.0*xd);
  yd = (double) m/nyg;
  yinc = 1.0/yd;
  y = yinc + (yd - 1.0)/(2.0*yd);
  for (iy = 0; iy < m; iy++) {
    x = xs;
    jy = y;
    dy = y - jy;
    dy1 = 1.0 - dy;
    if (jy < 1)
      j1 = j2 = 0;
    else if (jy >= nyg)
      j1 = j2 = nygm;
    else {
      j1 = jy - 1;
      j2 = j1 + 1;
    }
    jbase  = xgbase + j1 * 2 * nxg;
    jpbase = xgbase + j2 * 2 * nxg;
    for (ix = 0; ix < n; ix++) {
      jx = x;
      dx = x - jx;
      dx1 = 1.0 - dx;
      if (jx < 1)
	i1 = i2 = 0;
      else if (jx >= nxg)
	i1 = i2 = nxgm;
      else {
	i1 = jx - 1;
	i2 = i1 + 1;
      }
      w1 = dy1*dx1;
      w2 = dy1*dx;
      w3 = dy*dx1;
      w4 = dy*dx;
      i1 = 2*i1;
      i2 = 2*i2;
      xl = w1 * *(jbase+i1) + w2 * *(jbase+i2) + w3 * *(jpbase+i1)
	+ w4 * *(jpbase+i2) + (double) ix;
      i1 += 1;  i2 += 1;
      yl = w1 * *(jbase+i1) + w2 * *(jbase+i2) + w3 * *(jpbase+i1)
	+ w4 * *(jpbase+i2) + (double) iy;
      
      /* xl, yl is the place, now do a cubic interpolation for value */
      *out_x++ = xl;
      *out_y++ = yl;
      
      x += xinc;
    }
    y += yinc;
  }
}


void ana_stretch(const int ny, const int nx,  double *im, int npy, int npx, const double *gr,  double *out)
/* the call is MS = STRETCH( M2, DELTA)
   where M2 is the original array to be destretched, MS is the result, and
   DELTA is a displacement grid as generated by GRIDMATCH */
{
   double   *base; double *bb;

  int	iq, type, n, m, nxg, nyg, jy, j1, j2, nm1, nm2, mm1, mm2;
  int	nxgm, nygm, ix, iy, i1, i2, i3, i4, j3, j4, jx;
  double	xd, yd, xinc, yinc, y, x, xs, dy, dy1;
  double	dx, dx1, dx0, dx2, dx3, dx4, fn, fm, xq;
  double	w1, w2, w3, w4, xl, yl, c1, c2, c3, c4, b1, b2, b3, b4;
  const double 	*jpbase, *jbase, *xgbase;

  base = im; 

  n = nx;//argv[0]->value.arr->dim[0];
  m = ny;//argv[0]->value.arr->dim[1];
  nm1 = n - 1;
  mm1 = m - 1;
  nm2 = n -2;
  mm2 = m -2;
  fn = n;
  fm = m;

  xgbase = gr;//(double *) arg_struct[1].value.arr->data;
  nxg = npx;//arg_struct[1].value.arr->dim[1];
  nyg = npy;//arg_struct[1].value.arr->dim[2];
  nxgm = nxg - 1;
  nygm = nyg - 1;
  /* linearly interpolate the displacement grid values over array */
  /* similar to regrid3 in inner part */
  xd = (double) n/nxg;
  xinc = 1.0/xd;
  xs = xinc + (xd - 1.0)/(2.0*xd);
  yd = (double) m/nyg;
  yinc = 1.0/yd;
  y = yinc + (yd - 1.0)/(2.0*yd);
  for (iy = 0; iy < m; iy++) {
    x = xs;
    jy = y;
    dy = y - jy;
    dy1 = 1.0 - dy;
    if (jy < 1)
      j1 = j2 = 0;
    else if (jy >= nyg)
      j1 = j2 = nygm;
    else {
      j1 = jy - 1;
      j2 = j1 + 1;
    }
    jbase  = xgbase + j1 * 2 * nxg;
    jpbase = xgbase + j2 * 2 * nxg;
    for (ix = 0; ix < n; ix++) {
      jx = x;
      dx = x - jx;
      dx1 = 1.0 - dx;
      if (jx < 1)
	i1 = i2 = 0;
      else if (jx >= nxg)
	i1 = i2 = nxgm;
      else {
	i1 = jx - 1;
	i2 = i1 + 1;
      }
      w1 = dy1*dx1;
      w2 = dy1*dx;
      w3 = dy*dx1;
      w4 = dy*dx;
      i1 = 2*i1;
      i2 = 2*i2;
      xl = w1 * *(jbase+i1) + w2 * *(jbase+i2) + w3 * *(jpbase+i1)
	+ w4 * *(jpbase+i2) + (double) ix;
      i1 += 1;  i2 += 1;
      yl = w1 * *(jbase+i1) + w2 * *(jbase+i2) + w3 * *(jpbase+i1)
	+ w4 * *(jpbase+i2) + (double) iy;
      
      /* xl, yl is the place, now do a cubic interpolation for value */
      
      i2 = xl;
      j2 = yl;
      if (i2 >= 1 && i2 < nm2) {/* normal interior */
	i1 = i2 - 1;
	i3 = i2 + 1;
	i4 = i2 + 2;
	dx0 = xl - i2;
      } else {			/* edge cases */
	i2 = MIN(i2, nm1);
	i2 = MAX(i2, 0);
	xq = MIN(xl, fn-1.0);
	xq = MAX(xq, 0);
	dx0 = xq - i2;
	i1 = MIN(i2-1, nm1);
	i1 = MAX(i1, 0);
	i3 = MIN(i2+1, nm1);
	i3 = MAX(i3, 0);
	i4 = MIN(i2+2, nm1);
	i1 = MAX(i4, 0);
      }
      dx1 = 1.0 - dx0;
      dx2 = -dx0*0.5;
      dx3 = dx0*dx2;
      dx4 = 3.0*dx0* dx3;
      c1 = dx2*dx1*dx1;
      c2 = 1.0 - dx4 + 5.0*dx3;
      c3 = dx4 - (dx2 + 4.0*dx3);
      c4 = dx3*dx1;
      if (j2 >= 1 && j2 < mm2) { /* normal interior */
	j1 = j2 - 1;
	j3 = j2 + 1;
	j4 = j2 + 2;
	dx0 = yl - j2;
      } else {			/* edge cases */
	j2 = MIN(j2, mm1);
	j2 = MAX( j2, 0);
	xq = MIN(yl, fm-1.0);
	xq = MAX(xq, 0);
	dx0 = xq - j2;
	j1 = MIN(j2 - 1, mm1);
	j1 = MAX(j1, 0);
	j3 = MIN(j2 + 1, mm1);
	j3 = MAX(j3, 0);
	j4 = MIN(j2 + 2, mm1);
	j1 = MAX(j4, 0);
      }
      dx1 = 1.0 - dx0;
      dx2 = -dx0*0.5;
      dx3 = dx0*dx2;
      dx4 = 3.0*dx0*dx3;
      b1 = dx2*dx1*dx1;
      b2 = 1.0 - dx4 + 5.0*dx3;
      b3 = dx4 - (dx2 + 4.0*dx3);
      b4 = dx3*dx1;
      /* low corner offset */
      iq = i1 + j1*n;
      i2 = i2 - i1;
      i3 = i3 - i1;
      i4 = i4 - i1;
      j4 = (j4 - j3)*n;
      j3 = (j3 - j2)*n;
      j2 = (j2 - j1)*n;

      bb = base+iq;
      xq = b1*(c1 * *(bb) + c2 * *(bb+i2)
	       + c3 * *(bb+i3) + c4 * *(bb+i4));
      bb += j2;
      xq += b2*(c1 * *(bb) + c2 * *(bb+i2)
		+ c3 * *(bb+i3) + c4 * *(bb+i4));
      bb += j3;
      xq += b3*(c1 * *(bb) + c2 * *(bb+i2)
		+ c3 * *(bb+i3) + c4 * *(bb+i4));
      bb += j4;
      xq += b4*(c1 * *(bb) + c2 * *(bb+i2)
		+ c3 * *(bb+i3) + c4 * *(bb+i4));
      *out++ = xq;
      
      x += xinc;
    }
    y += yinc;
  }
}

/*-----------------------------------------------------------------*/
inline void gwind0(double* __restrict__ gwx, double* __restrict__ gwy, double gwid, int nxa, int nxb, int nya, int nyb)
{
  double	wid, xcen, ycen, xq;
  int	i;
  
  wid = gwid*0.6005612;
  if (wid > 0) {
    xcen = (nxa + nxb)/2;
    ycen = (nya + nyb)/2; 
    for (i = nxa; i <= nxb; i++) {
      xq = (i - xcen)/wid;
      gwx[i] = exp(-(xq*xq));
    } 
    for (i = nya; i <= nyb; i++) {
      xq = (i - ycen)/wid;
      gwy[i] = exp(-(xq * xq));
    } 
  } else {
    for (i = nxa; i <= nxb; i++) 
      gwx[i] = 1.0;
    for (i = nya; i <= nyb; i++)
      gwy[i] = 1.0; 
  }
}

/*------------------------------------------------------------------------- */

inline double averag(double *p, int nxa, int nxb, int nya,
	     int nyb, int nxs, int nys, int idx,
	     int idy, double *gx, double *gy)
/* finds weighted average of array m over the block defined */
{
  int	nxc, nxd, nyc, nyd, jj;
  double	sum, sumg, sumx, sumgx;

  /* fix limits so sum doesn't run off edge of image */
  nxc = (nxa + idx < 0)? -idx: nxa;
  nyc = (nya + idy < 0)? -idy: nya;
  nxd = (nxb + idx > nxs)? nxs - idx: nxb;
  nyd = (nyb + idy > nys)? nys - idy: nyb;
  sum = sumg = sumgx = 0.0; 
  for (int i = nxc; i < nxd; i++)
    sumgx += gx[i];
  /* weighted sum in window */
  for (int j = nyc; j < nyd; j++) {
    sumx = 0.0;
    jj = idx + nxs*(j + idy);
    for (int i = nxc; i < nxd; i++)
      sumx += gx[i]*p[i + jj];
    sum += gy[j]*sumx;
    sumg += gy[j]*sumgx;
  } /* end of j loop */
  sum /= sumg;
  return sum;
}

/*------------------------------------------------------------------------- */

inline void unbias(double *m1, double *m2, int nxa, int nxb,
		   int nya, int nyb, int nxs, int nys,
		   double *gx, double *gy, double *av1, double *av2, double *cx,
		   double *cy, double *cxx, double *cxy, double *cyy,
		   int idelx, int idely)
{
  double	t0, t1, t2, t3, t4, t5;
  double averag(double *p, int nxa, int nxb, int nya,
		int nyb, int nxs, int nys, int idx,
		int idy, double *gx, double *gy);
  
  /*  find weighted means of m1 & m2 over the window 
      sets up quadratic fit to average of m2 as a fcn. of offsets */
  *av1 = averag(m1, nxa, nxb, nya, nyb, nxs, nys, 0, 0, gx, gy); 
  t0 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx, idely, gx, gy); 
  t1 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx + 1, idely, gx, gy); 
  t2 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx - 1, idely, gx, gy); 
  t3 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx, idely + 1, gx, gy); 
  t4 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx, idely - 1, gx, gy); 
  t5 = averag(m2, nxa, nxb, nya, nyb, nxs, nys, idelx + 1, idely + 1, gx, gy); 
  *av2 = t0; 
  *cx = 0.5*(t1 - t2); 
  *cy = 0.5*(t3 - t4); 
  *cxx = 0.5*(t1 - 2*t0 + t2); 
  *cyy = 0.5*(t3 - 2*t0 + t4); 
  *cxy = t5 + t0 - t1 - t3;
}

/*------------------------------------------------------------------------- */
inline int getmin(double *p, double *x0, double *y0)
{
  double	f11, f12, f13, f21, f22, f23, f31, f32, f33;
  double	fx, fy, t, fxx, fyy, fxy;

  /* find the min, p points to a 3x3 array */
  f11 = *p++;	f21 = *p++;	f31 = *p++;
  f12 = *p++;	f22 = *p++;	f32 = *p++;
  f13 = *p++;	f23 = *p++;	f33 = *p++;

  fx = 0.5 * ( f32 - f12 );	fy = 0.5 * ( f23 - f21 );
  t = 2.* ( f22 );
  fxx =  f32 + f12 - t;	fyy = f23 + f21 - t;
  /* find in which quadrant the minimum lies */
  if (f33 < f11) {
    if (f33 < f31) {
      if (f33 < f13)
	fxy = f33+f22-f32-f23;
      else
	fxy = f23+f12-f22-f13;
    } else {
      if (f31 < f13)
	fxy = f32+f21-f31-f22;
      else
	fxy = f23+f12-f22-f13;
    }
  } else { 
    if (f11 < f31) {
      if (f11 < f13) 
	fxy = f22+f11-f21-f12;
      else
	fxy = f23+f12-f22-f13;
    } else {
      if (f31 < f13)
	fxy = f32+f21-f31-f22;
      else
	fxy = f23+f12-f22-f13;
    }
  }
  t = -1./(fxx *fyy - fxy *fxy);
  *x0 = t * (fx * fyy - fy * fxy);
  *y0 = t * (fy * fxx - fx * fxy);
  if (ABS(*x0) >= 0.75 || ABS(*y0) >= 0.75) {
    *x0 = -fx/fxx;
    *y0 = -fy/fyy;
  }
  return 1;
}

/*------------------------------------------------------------------------- */

inline double resid(double *m1, double *m2, int idx, int idy, int nxa,
	    int nxb, int nya, int nyb, int nxs,
	    int nys, int ndmx, double *gx, double *gy, double bs)
{
  int	nxc, nxd, nyc, nyd, nx, ny;
  double 	*p1, *p2, *ps;
  double	sum, sumx, t, ndmx2;
  int	i, j;
  double   sumg;
  static int	mxc, mxd, myc, myd;
  static double	gsum;

  /*set up limits */
  nxc = nxa;
  if (nxc + idx < 0)
    nxc = -idx;
  nyc = nya;
  if (nyc + idy < 0)
    nyc = - idy;
  nxd = nxb;
  if (nxd + idx >= nxs)
    nxd = nxs - idx - 1;
  nyd = nyb;
  if (nyd + idy >= nys)
    nyd = nys - idy - 1;
  sum = sumg = 0.0;

  nx = nxd - nxc +1;
  p2 = gy + nyc;
  ps = gx + nxc;

  if (nxc != mxc || nxd != mxd || nyc != myc || nyd != myd) {
    /* sum gaussians over rectangle to get normalization */
    /* (only if limits change)*/
    j = nyd -nyc + 1;
    if (j <= 0 || nxd - nxc + 1 <= 0)
      return -1;		/* added 19feb95 LS */
    while (j) {
      i = nx;
      p1 = ps;
      while (i) {
	sumg += (*p1++) * (*p2);
	i--;
      }
      p2++;
      j--;
    }
    gsum = sumg;
    mxc = nxc;
    mxd = nxd;
    myc = nyc;
    myd = nyd;
  } else
    sumg = gsum;

  m1 += nyc*nxs + nxc;
  m2 += (nyc + idy)*nxs + nxc + idx;
  ny = nxs - nx;		/*residual increment after inner loop */
  p2 = gy + nyc;
  
  /* now the loop to compute the residual */
  j = nyd - nyc +1;
  ndmx2 = ndmx*ndmx;
  while (j) {
    i = nx;
    p1 = ps;
    sumx = 0.0;
    while (i) {
      t = *m1++ - *m2++;
      t = t + bs;
      t = t*t;
      t = MIN(t, ndmx2);
      sumx += (*p1++) * t;
      i--;
    }
    sum += (*p2++) * sumx;
    m1 += ny;
    m2 += ny;
    j--;
  }
  /*return normalized residual */
  sum /= sumg;
  return sum;
}

/*-----------------------------------------------------------------*/

inline void match_1(double *p1, double *p2, int nxa, int nxb,
		    int nya, int nyb, int nx, int ny,
		    double *gwx, double *gwy, double *xoffset, double *yoffset, int &stretch_clip, int &badmatch)
{
  int idelx, idely, i, j, k, ndmx = 1000, done[9]={};
  int	di, dj, in, jn, iter, dd, badflag = 0;
  double	av1, av2, cx, cy, cxx, cxy, cyy, avdif, t, res[9], buf[9], t1, t2;
  static int	itmax = 40;
  
  idelx = rint(*xoffset);
  idely = rint(*yoffset); 
  unbias(p1, p2, nxa, nxb, nya, nyb, nx, ny, gwx, gwy, &av1, &av2, &cx, &cy,
	 &cxx, &cxy, &cyy, idelx, idely);
  /* look at a 3x3 matrix of residuals centered at 0 offset, find the location
     of the minimum, if not at center, then look at new 3x3 centered
     on the edge minimum; repeat until min at center */
  iter = itmax;
  badflag = 0;
  while (iter--) {
    for (k = 0; k < 9; k++) {
      if (done[k] == 0) {
	i = idelx + (k % 3) - 1;
	j = idely + (k / 3) - 1;
	avdif = av2 +  i*cx + j*cy + i*i*cxx + i*j*cxy + j*j*cyy - av1;
	res[k] = resid(p1, p2, i, j, nxa, nxb, nya, nyb, nx, ny, ndmx,
		       gwx, gwy, avdif);
      }
    }
    t = res[0];
    i = 0;
    for (k = 1; k < 9; k++) 
      if (res[k] < t) {
	t = res[k];
	i = k;
      }
    if (t < 0) {		/* added LS 19feb95 */
      printf("match - ran out of data at edge\n");
      badflag = 1;
      break;
    }
    idelx += (i % 3) - 1;
    idely += (i / 3) - 1;
    /* check if we have gone too far */
    if (ABS(idelx) > stretch_clip || ABS(idely) > stretch_clip) {
      badflag++;
      break;
    }
    if (i == 4)
      break;			/* done if it is the center one */
    /* not in center, shuffle what we have to put the edge min in center */
    di = (i % 3) - 1;
    dj = (i / 3) - 1;
    dd = dj * 3 + di;
    for (k = 0; k < 9; k++) {
      in = k%3 + di;
      jn = k/3 + dj;
      if (in >= 0 && jn >= 0 && in < 3 && jn < 3) { /* in range */
	done[k] = 1;
	buf[k] = res[k + dd];
      } else 
	done[k] = 0;		/* not in range, mark undone */
    }
    for (k = 0; k < 9; k++)
      res[k] = buf[k];		/* put back in res array */
  } /* end of iter while */
  /* done or reached itmax, which ? */
  if (iter <= 0) {
    badflag++;
  }
  if (badflag) {
    badmatch++;
    *xoffset = *yoffset = 0;
    return;
  }
				/* must have been OK so far */
  getmin(res, &t1, &t2);
  *xoffset = idelx + t1;
  *yoffset = idely + t2;
}

/*-----------------------------------------------------------------*/

void ana_gridmatch(int const ny, int const nx, double* __restrict__ p1, double* __restrict__ p2,
		   int const nyg, int const nxg, int* __restrict__ gy, int* __restrict__ gx, int const dy, int const dx,
		   double const gwid, int stretch_clip, double* __restrict__ out)
/* the call is offsets = gridmatch(m1,m2,gx,gy,dx,dy,gwid,stretch_clip)
	where	m1 = reference input image
	m2 = image to compare with m1, m1 and m2 must be same size
	gx = array of x gridpoints
	gy = array of y gridpoints, gx and gy must have same size
	dx and dy are the window size, and gwid is the gaussian mask width
        stretch_clip = maximum allowed displacement before a bad match
          is declared

  Authors:  Richard A. Shine (original)
            Louis H. Strous (port from ANA to IDL)
            Lockheed-Martin Advanced Technology Center, Palo Alto, CA, USA
*/
{
  double	        xoffset, yoffset,*gwx, *gwy;
  int     dx2, dy2, nc,  i1, i2, j1, j2;
  int      dims[3] = {2,nxg,nxg};


  if (stretch_clip < 2)
    stretch_clip = 2;
  stretch_clip--;

  /* prepare the gaussian kernel */
  gwx = (double *) malloc((nx + ny)*sizeof(double));
  gwy = gwx + nx;
  nc = nxg*nyg;			/* the number of subimages */
  dx2 = dx/2;
  dy2 = dy/2;
  int badmatch = 0;
  
  while (nc--) {		/* loop over all subimages */
    i1 = *gx - dx2;		/* lower x coordinate */
    if (i1 < 0)			/* outside array? */
      i1 = 0;
    i1++;
    i2 = *gx++ + dx2;		/* upper x coordinate */
    if (i2 > nx)
      i2 = nx;
    j1 = *gy - dy2;		/* lower y coordinate */
    if (j1 < 0)
      j1 = 0;
    j1++;
    j2 = *gy++ + dy2;		/* upper y coordinate */
    if (j2 > ny)
      j2 = ny;
    xoffset = yoffset = 0;
    i1--; i2--; j1--; j2--;
    gwind0(gwx, gwy, gwid, i1, i2, j1, j2); /* get gaussian kernels */
    match_1(p1, p2, i1, i2, j1, j2, nx, ny, gwx, gwy,
	    &xoffset, &yoffset, stretch_clip, badmatch); /* get offsets */
    *out++ = xoffset;
    *out++ = yoffset;
  }

  free(gwx);
  //  return result;
}
