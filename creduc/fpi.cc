#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fftw3.h"
#include "fpi.h"
#include "mymath.h"


const double wng[5] = {0.11846344252809,0.23931433524968,0.28444444444444,0.23931433524968,0.11846344252809};
const double calp[5] = {0.04691007703067,0.23076534494716,0.50,0.76923465505284,0.95308992296933};
const int nang = 5;
const double pi2 = 2.0 * 3.14159265358979323846;

void init_fpi(fpi_t &fpi, int32_t line){
  fpi.calp = NULL;
  fpi.c1 = NULL;
  fpi.c2 = NULL;
  fpi.line = line;

  switch (fpi.line)
    {
    case 6302:
      fpi.w0   = 6302.4940;
      fpi.rhr  = 0.9350;
      fpi.rlr  = 0.8380;
      fpi.shr  = 787.0E4;
      fpi.slr  = 295.5E4;
      fpi.dhr  = 0.267;
      fpi.dlr  = 0.731;
      fpi.dw   = 0.005;
      break;
    case 8542:
      fpi.w0   = 8542.091;
      fpi.rhr  = 0.930E0;
      fpi.rlr  = 0.852E0;
      fpi.shr  = 787.0E4;
      fpi.slr  = 295.5E4;
      fpi.dhr  = 0.363;
      fpi.dlr  = 0.991;
      fpi.dw   = 0.010;
      break;
    case 6173:
      fpi.w0   = 6173.3356;
      fpi.rhr  = 0.9380;
      fpi.rlr  = 0.8469;
      fpi.shr  = 787.0E4;
      fpi.slr  = 295.5E4;
      fpi.dhr  = 0.267;
      fpi.dlr  = 0.731;
      fpi.dw   = 0.005;
      break;
    default:
      fprintf(stderr,"init_fpi : ERROR -> data not found for spectral line %d\n",fpi.line);
      exit(0);
      break;
    }
  fprintf(stderr,"init_fpi : Initializing FPI for spectral line %d\n",line);
  fpi.scl = true;

  //
  // Center HRE & LRE at w0
  //
  float nr = 0.5 + fpi.shr / (fpi.w0 * 0.5); //number of wavs HRE
  fpi.shr = (int)(nr) * fpi.w0 * 0.5;
  
  nr = 0.5 + fpi.slr / (fpi.w0 * 0.5); //number of wavs HRE
  fpi.slr = (int)(nr) * fpi.w0 * 0.5;
  
  
  //
  // compute quadrature for integral
  //
  fpi.calp = new double [nang];
  for(int jj=0;jj<nang;jj++) fpi.calp[jj] = cos(sqrt(calp[jj] * sqr<double>((double)0.5 / 165.0))) * pi2;

}

void init_fftw(fpi_t &fpi, int nmean, float *xl, float *yl){
  
  //
  // Use fine grid for the convolutions (dw defines the new wavelength step)
  //
  int np = ((xl[nmean-1] - xl[0]) / fpi.dw);
  np += 1;
  if((np/2)*2 == np) np += 1;
  fpi.npad = np;//*2-1;



  //
  // Allocate arrays
  //
  fpi.xlp = new double [fpi.npad];
  fpi.ylp = new double [fpi.npad];

  fpi.tr = new double [fpi.npad];
  fpi.tw = new double [fpi.npad];
  //
  fpi.ft = new double complex [fpi.npad/2+2];
  fpi.otf = new double complex [fpi.npad/2+2];

  //
  // Init plan for the spectrum and PSF
  //
  fprintf(stderr,"fft_init : initializing FFTW plans with npad=%d ...", fpi.npad);
  //
  fpi.fplan = fftw_plan_dft_r2c_1d(fpi.npad, (double*)fpi.ylp,(fftw_complex*)(fpi.ft), FFTW_MEASURE);
  fpi.otfplan=fftw_plan_dft_r2c_1d(fpi.npad, (double*)fpi.tr, (fftw_complex*)(fpi.otf),FFTW_MEASURE);
  fpi.bplan = fftw_plan_dft_c2r_1d(fpi.npad, (fftw_complex*)(fpi.otf),(double*)fpi.ylp,FFTW_MEASURE);
  fprintf(stderr,"done\n");

  //
  // Interpolate to finer grid
  //
  for(int ii=0;ii<fpi.npad;ii++) {
    fpi.xlp[ii] = ii*fpi.dw + xl[0];
    fpi.tw[ii] = (ii-fpi.npad/2)*fpi.dw + fpi.w0;
  }
  intep<double>(nmean,xl,yl,fpi.npad,fpi.xlp,fpi.ylp);

  //
  // Pad array with the largest value (continuum?)
  //
  float tmp = yl[0];
  for(int ii=1;ii<nmean;ii++) if(yl[ii] > tmp) tmp = yl[ii];
  for(int ii=np;ii<fpi.npad;ii++) fpi.ylp[ii] = tmp;

  //for(int ii=np;ii<np+np/2;ii++) fpi.ylp[ii] = yl[nmean-1]; //padding
  //for(int ii=np+np/2;ii<fpi.npad;ii++) fpi.ylp[ii] = yl[0]; //padding
  //for(int ii=0;ii<fpi.npad;ii++) fprintf(stderr,"%f %f\n",fpi.xlp[ii],fpi.ylp[ii]);
  

  //
  // execute FFT for the spectrum and store it
  // 
  fftw_execute(fpi.fplan);
}

void dual_fpi(fpi_t *fpi, double erh){
  //
  // reflectivities + error
  //
  fpi->orh = erh;
  double mrhr = erh;
  double mrlr = fpi->rlr; //+ fpi->erl;

  //
  // Finesse
  //
  double fhr = 4.0 * mrhr / sqr<double>(1.0 - mrhr);
  double flr = 4.0 * mrlr / sqr<double>(1.0 - mrlr);
  
  
  //
  // Phase difference for each ethalon
  //
  double phr = pi2 * fpi->shr; //* cos(angle) = 1.0 in this case
  double plr = pi2 * fpi->slr; //* cos(angle) = 1.0 in this case 
  

  //
  // Transmission peaks
  // 
  int npad2 = fpi->npad/2;
  for(int jj = 0;jj<=npad2;jj++){
    //
    fpi->tr[jj] =  1.0 / (1.0 + fhr * sqr<double>(sin(phr / fpi->tw[jj+npad2])));
    fpi->tr[jj] *= 1.0 / (1.0 + flr * sqr<double>(sin(plr / fpi->tw[jj+npad2])));

    //
  }
  for(int jj = 1;jj<=npad2;jj++){
    //
    //fpi->tr[jj] =  1.0 / (1.0 + fhr * sqr<double>(sin(phr / fpi->tw[jj-npad2])));
    //fpi->tr[jj] *= 1.0 / (1.0 + flr * sqr<double>(sin(plr / fpi->tw[jj-npad2])));
    fpi->tr[jj+npad2] = fpi->tr[npad2-jj+1]; //the profile is symmetric here.
    //
    
  }
  //
  // Normalize to the total number of elements (the area of the 
  // CRISP profile is accounted for when computing the ratio).
  //
  if(fpi->scl){
    for(int ww=0;ww<fpi->npad;ww++) fpi->tr[ww] *= (double)fpi->npad;
    fpi->scl=false;
  }

  //
  // FFT
  // 
  fftw_execute(fpi->otfplan);  
}

void bezier3_control(int n, double *x, double *y, double *&c1, double *&c2){
  
  // Allocate arrays?
  if(c1 == NULL) c1 = new double [n-1];
  if(c2 == NULL) c2 = new double [n-1];

  // Init der in the first interval
  c1[0] = y[0] +  (y[1] - y[0]) / 3.0; 

  double dx1, der1;
  for(int k=1;k<n-1;k++){
    int k1 = k-1;

    // Derivatives
    double dx = x[k] - x[k1];
    double der = (y[k] - y[k1]) / dx;
    dx1 = x[k+1] - x[k];
    der1 = (y[k+1] - y[k]) / dx1;

    // If Max/Min then set derivative to zero, else compute it.
    double yp = 0.0;
    if(der*der1 > 0){
      double lambda = (1.0 + dx1 / (dx1 + dx)) / 3.0;
      yp = (der * der1) / (lambda * der1 + (1.0 - lambda) * der);      
    }

    // Control points
    c1[k] = y[k] + dx1 * yp / 3.0;
    c2[k1] = y[k] - dx * yp / 3.0;
  }

  // Fill missing element in c2
  c2[n-2] = y[n-1] - dx1 * der1 / 3.0;


}

void bezier3(int n, double *x, double *y, int np, double *xp, double *yp,double *c1, double *c2){

  int pos = 0;
  for(int k=0;k<n-1;k++){
    double dx = x[k+1] - x[k];

    for(int jj=pos;jj<np;jj++){
      if(xp[jj] >= x[k] && xp[jj] < x[k+1]){
	pos += 1;
	double u = (xp[jj] - x[k]) / dx;
	double u1 = 1.0 - u;
	
	yp[jj] = y[k] * u1*u1*u1 + y[k+1]*u*u*u + 
	  3.0 * c1[k] * u * u1*u1 + 3.0 *c2[k]* u*u * u1;
      }
    }
  }

  //
  // out-of-bounds
  //
  double a0 = (y[1]-y[0])/(x[1]-x[0]);
  double b0 = y[0] - a0 * x[0];
  double a1 =  (y[n-1]-y[n-2])/(x[n-1]-x[n-2]);
  double b1 = y[n-1] - a1 * x[n-1];
  for(int k=0;k<np;k++){
    if(xp[k] >= x[n-1]) yp[k] = a1*xp[k] + b1;
    if(xp[k] < x[0]) yp[k] = a0 * xp[k] + b0;
  }

}

void clean_fpi(fpi_t &fpi){
  fftw_destroy_plan(fpi.otfplan);
  fftw_destroy_plan(fpi.fplan);
  fftw_destroy_plan(fpi.bplan);
  delete [] fpi.calp;
  delete [] fpi.tw;
  delete [] fpi.tr;
  delete [] fpi.otf;
  delete [] fpi.ft;
  delete [] fpi.ylp;
  delete [] fpi.xlp;
  delete [] fpi.c1;
  delete [] fpi.c2;
  delete [] fpi.imean;
}
