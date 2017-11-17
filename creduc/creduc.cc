/*
  C routines for the CRISPRED pipeline
  -
  AUTHORS: Jaime de la Cruz Rodriguez (IFA-UU)
           Michiel van Noort (MPS)
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include <omp.h>
#include <algorithm>
//#include <complex.h>
#include <iostream>
#include "types.h"
#include "fftw3.h"
#include "mymath.h"
#include "mmem.h"
#include "mpfit.h"
#include "physc.h"
#include "fpi.h"

using namespace std;
/*

#include "Genetic.h"
#include "Tools.h"
using namespace ROYAC;
*/
double sqr(double var){
  return var * var;
}
double sqr(double &var){
  return var * var;
}

inline int32_t ind(int x, int y, int nx){
  return (y*nx + x);
}
//
// descatter
//
void descatter(int nx, int ny, float *timg, float *tfgain, float *tfpsf, float *res, int nthreads, int verbose){
  //
  // Init variables
  //
  char inam[] = "descatter :";
  float **img = var2dim<float>(timg,ny,nx);
  float **fgain = var2dim<float>(tfgain,ny,nx);
  float **fpsf = var2dim<float>(tfpsf,ny,nx);
  //
  double **psf=mat2d<double>(2*ny, 2*nx);
  memset(psf[0], 0, 4 * nx * ny * sizeof(double));
  double **gain = mat2d<double>(ny, nx);
  //
  int nx2 = nx/2;
  int ny2 = ny/2;
  double suma = 0.0;
  for(int y=0;y<ny;++y) for(int x=0;x<nx;++x) gain[y][x]=fgain[y][x];
  for(int y=0;y<ny;++y) for(int x=0;x<nx;++x) suma+=fpsf[y][x];
  for(int y=0;y<ny;++y) for(int x=0;x<nx;++x) psf[y+ny2][x+nx2]=fpsf[y][x]/suma;
  psf_reorder<double>(psf, 2*ny);
  if(verbose){
    fprintf(stderr, "%s nx=%d, ny=%d, nthreads=%d, psf_integral=%g \n",inam,nx, ny, nthreads, suma );
  }
  //
  // Init FFTW threads
  //
  fftw_init_threads();
  fftw_plan_with_nthreads(nthreads);
  //
  // Compute otf
  //
  complex_t **otf = mat2d<complex_t>(2*ny, nx+1);
  memset(otf[0],0,2*ny*(nx+1)*sizeof(complex_t));
  fftw_plan otfplan=fftw_plan_dft_r2c_2d(2*ny,2*nx,psf[0],(fftw_complex*)(otf[0]),FFTW_ESTIMATE);
  //fftw_execute_dft_r2c(otfplan,psf[0],(fftw_complex*)(otf[0]));
  fftw_execute(otfplan);
  //
  // Temporary arrays for FFTW
  //
  double **im = mat2d<double>(2*ny, 2*nx);
  double **cim = mat2d<double>(2*ny, 2*nx);
  memset(im[0],0,4*nx*ny*sizeof(double));
  complex_t **ft=mat2d<complex_t>(2*ny,nx+1);
  //
  // FFTW Plans
  //
  fftw_plan fplan=fftw_plan_dft_r2c_2d(2*ny,2*nx,im[0],(fftw_complex*)(ft[0]),FFTW_ESTIMATE);
  fftw_plan bplan=fftw_plan_dft_c2r_2d(2*ny,2*nx,(fftw_complex*)(ft[0]),cim[0],FFTW_ESTIMATE);
  //
  double np = nx * ny;
  double np4 = 4.0 * np;
  for(int y=0;y<ny;++y) for(int x=0;x<nx;++x) {
      im[y+ny2][x+nx2]=img[y][x]/(1.0+gain[y][x] * gain[y][x]);
    }
  
  //
  // Iterate solution
  //
  double cs = 1E10, ocs;
  int i = 0, iter_max = 50;
  //
  do{
    memcpy(cim[0],im[0],4*nx*ny*sizeof(double));
    for(int y=0;y<ny;++y) for(int x=0;x<nx;++x) cim[y+ny2][x+nx2]*=gain[y][x];
    //
    // Forward transform
    //
    fftw_execute_dft_r2c(fplan,cim[0],(fftw_complex*)(ft[0]));
    //
    // convolve with OTF
    //
    for(int y=0;y<2*ny;++y) 
      for(int x=0;x<nx+1;++x){
	double tmp=ft[y][x].re;
	ft[y][x].re=otf[y][x].re*tmp-otf[y][x].im*ft[y][x].im;
        ft[y][x].im=otf[y][x].im*tmp+otf[y][x].re*ft[y][x].im;
      }
    //
    // Convert back to image space
    //
    fftw_execute_dft_c2r(bplan,(fftw_complex*)(ft[0]),cim[0]);
    //
    // Normalize result (FFTW does not include the 1/N factor)
    //
    for(int y=0;y<2*ny;++y) for(int x=0;x<2*nx;++x) cim[y][x]/=np4;
    //
    for(int y=0;y<ny;++y) for(int x=0;x<nx;++x) cim[y+ny2][x+nx2]*=gain[y][x];
    //
    // Compute Chisq
    //
    cs=0.0;
    for(int y=0;y<ny;++y) 
      for(int x=0;x<nx;++x){
	double bla=img[y][x]-cim[y+ny2][x+nx2];
        cs+=sqr(im[y+ny2][x+nx2]-bla);
        im[y+ny2][x+nx2]=bla;
      }
    cs/= (double)(np);
    if(verbose) fprintf(stderr, "%s iter = %d, ChiSq = %e \n", inam, i,cs);
  }while((cs>1E-8)&&(++i<iter_max));
  //
  // Destroy FFTW plans
  //
  fftw_destroy_plan(bplan);
  fftw_destroy_plan(fplan);
  fftw_destroy_plan(otfplan);
  //
  //  Write result to res array
  //
  for(int y=0;y<ny;++y) for(int x=0;x<nx;++x) res[(y*nx+x)] = im[y+ny2][x+nx2];
  //
  // Clean-up
  //
  //del_mat2d<float>(img);
  //del_mat2d<float>(fgain);
  //del_mat2d<float>(fpsf);
  del_mat2d<double>(im);
  del_mat2d<double>(cim);
  del_mat2d<double>(psf);
  del_mat2d<double>(gain);
  del_mat2d<complex_t>(otf);
  del_mat2d<complex_t>(ft);
  //
  delete [] img;
  delete [] fpsf;
  delete [] fgain;
  //
}
void convolve(int inx, int iny, int pnx, int pny, float *img, float *fpsf, float *res1, int nthreads, int verbose){
  char inam[] = "convolve :";
  //
  // reshape 2D arrays
  //
  //float **img = var2dim<float>(timg,iny,inx);
  //float **fpsf = var2dim<float>(tfpsf,pny,pnx);
  //
  // Allocate arrays with extra-space for padding
  //
  int npadx = inx + pnx;
  int npady = iny + pny;
  //
  if(verbose){
    fprintf(stderr, "%s Image -> nx=%d, ny=%d\n",inam,inx, iny);
    fprintf(stderr, "%s PSF -> nx=%d, ny=%d\n",inam,pnx, pny);
    fprintf(stderr, "%s FFT -> npad_x=%d, npad_y=%d, nthreads=%d\n",inam, npadx, npady, nthreads);
  }
  //
  double **pimg = mat2d<double>(npady, npadx);
  double **ppsf= mat2d<double>(npady, npadx);
  memset(ppsf[0], 0, npadx * npady * sizeof(double));
  //
  // Pad with the mean of the image
  //
  double suma = 0.0;
  for(int y=0;y<iny;++y) for(int x=0;x<inx;++x) suma += img[ind(x,y,inx)];
  suma/= (double)(inx*iny);
  for(int y=0;y<npady;++y) for(int x=0;x<npadx;++x) pimg[y][x] = suma;
  //
  // Copy image onto padded array
  //
  for(int y=0;y<iny;++y) for(int x=0;x<inx;++x) pimg[y][x] = img[ind(x,y,inx)];
  //
  // Copy PSF to new array
  //
  int pny2 = pny/2;
  int pnx2 = pnx/2;
  //
  suma = 0;
  for(int y=0;y<pny;++y) for(int x=0;x<pnx;++x) suma += fpsf[ind(x,y,pnx)];

  for(int y=0;y<pny2;++y) for(int x=0;x<pnx2;++x) ppsf[y][x] = fpsf[ind(x+pnx2,y+pny2,pnx)]/suma; // first
  //
  int x0 = npadx - pnx2;
  int y0 = npady - pny2;
  for(int y=0;y<pny2;++y) for(int x=0;x<pnx2;++x) ppsf[y+y0][x+x0] = fpsf[ind(x,y,pnx)]/suma;  // Fourth
  for(int y=0;y<pny2;++y) for(int x=0;x<pnx2;++x) ppsf[y][x+x0] = fpsf[ind(x,y+pny2, pnx)]/suma; // third
  for(int y=0;y<pny2;++y) for(int x=0;x<pnx2;++x) ppsf[y+y0][x] = fpsf[ind(x+pnx2,y,pnx)]/suma; // second
  //
  // FFTW plans, otf and other vars
  //
  fftw_init_threads();
  fftw_plan_with_nthreads(nthreads);
  //
  complex_t **ft = mat2d<complex_t>(npady, npadx/2+1);
  complex_t **otf = mat2d<complex_t>(npady, npadx/2+1);
  memset(otf[0],0,npady*(npadx/2+1)*sizeof(complex_t));
  //
  fftw_plan otfplan=fftw_plan_dft_r2c_2d(npady,npadx,ppsf[0],(fftw_complex*)(otf[0]),FFTW_ESTIMATE);
  fftw_plan fplan=fftw_plan_dft_r2c_2d(npady,npadx,pimg[0],(fftw_complex*)(ft[0]),FFTW_ESTIMATE);
  fftw_plan bplan=fftw_plan_dft_c2r_2d(npady,npadx,(fftw_complex*)(ft[0]),pimg[0],FFTW_ESTIMATE);
  //
  // Forward transforms
  //
  if(verbose) fprintf(stderr, "%s convolving ... ", inam);
  //
  fftw_execute(otfplan);
  fftw_execute(fplan);
  //
  // Convolve in Fourier space
  //
  for(int y=0;y<npady;++y) for(int x=0;x<npadx/2+1;++x){
      double tmp=ft[y][x].re;
      ft[y][x].re=otf[y][x].re*tmp-otf[y][x].im*ft[y][x].im;
      ft[y][x].im=otf[y][x].im*tmp+otf[y][x].re*ft[y][x].im;
    }
  fftw_execute(bplan);
  //
  // Destroy FFTW plans
  //
  fftw_destroy_plan(bplan);
  fftw_destroy_plan(fplan);
  fftw_destroy_plan(otfplan);
  //
  // Copy to output array
  //
  double np4 = npadx * npady;
  for(int y=0;y<iny;++y) for(int x=0;x<inx;++x) res1[ind(x,y,inx)] = pimg[y][x] / np4;
  //
  // Clean-up
  //
  del_mat2d<double>(pimg);
  del_mat2d<double>(ppsf);
  del_mat2d<complex_t>(otf);
  del_mat2d<complex_t>(ft);
  //
  if(verbose) fprintf(stderr, "done \n");
  //
}

double get_polcomp2(int nwav, int npar, float wav, double *pars){
  //
  double res = 1.0;
  double twav = 1.0;
  //
  if(npar <= 3) return res;
  //
  int ii0 = 3;
  for(int ii=ii0;ii<npar;ii++){
    twav *= wav;
    res += pars[ii] * twav;
  }
  //
 
  return res;
}
double get_polcomp(int &nwav, int &npar, float wav, double *pars){
  //
  double res = 1.0;
  double twav = 1.0;
  //
  if(npar <= 2) return res;
  //
  int ii = 2;
  while(ii<npar){
    twav *= wav;
    res += pars[ii] * twav;
    ii += 1;
  }
  //
  return res;
}
int fitgain_model(int nwav, int npar, double *pars, double *dev, double **derivs, void *tmp1){
  //
  fgd_t *tmp = (fgd_t*) tmp1;
  //
  // interpolate to shifted grid (cavity error)
  // 
  double *twav = new double [nwav];
  for(int ii = 0;ii<nwav;++ii) twav[ii] = tmp->wav[ii] - pars[1];
  //
  double *spec = intep<double>(tmp->nmean, tmp->xl, tmp->yl, nwav, twav);
  //
  // Compute difference between model and observations
  //
  for(int ii=0;ii<nwav;++ii) dev[ii] = (pars[0] * spec[ii] * get_polcomp(nwav, npar, tmp->wav[ii], pars)) - tmp->idat[ii]; 
  //
  // Clean-up
  //
  delete [] spec;
  delete [] twav;
  //
  return 0;
}
void get_ratio(int &nwav, int &npar, int &nmean, float *wav, float *dat, double *pars, float *xl, float *yl, float *ratio){
  double *twav = new double [nwav];
  for(int ii = 0;ii<nwav;++ii) twav[ii] = wav[ii] - pars[1];
  //
  double *spec = intep<double>(nmean, xl, yl, nwav, twav);
  //
  // Compute difference between model and observations
  //
  for(int ii=0;ii<nwav;++ii) ratio[ii] = dat[ii] / (pars[0] * spec[ii] * get_polcomp(nwav, npar, wav[ii], pars));
  //
  delete [] twav;
  delete [] spec;
}
void get_ratio2(int &nwav, int &npar, int &nmean, float *wav, float *dat, double *pars, float *xl, float *yl, float *ratio){
  double *twav = new double [nwav];
  for(int ii = 0;ii<nwav;++ii) twav[ii] = wav[ii] - pars[1];
  //
  double *spec = intep<double>(nmean, xl, yl, nwav, twav);
  //
  // Compute difference between model and observations
  //
  for(int ii=0;ii<nwav;++ii) ratio[ii] = dat[ii] / (pars[0] * spec[ii] * get_polcomp2(nwav, npar, wav[ii], pars));
  //
  delete [] twav;
  delete [] spec;
}
//
void fitgain(int nwav, int nmean, int npar, int npix, float *xl, float *yl, float *wav, float *dat1, double *pars1, float *ratio1, int nt){
  fprintf(stderr, "cfitgain : nwav=%d, npar=%d, nmean=%d, npix=%d \n", nwav, npar, nmean, npix);
  // return;
  //
  // Reshape vars
  //
  float **dat = var2dim<float>(dat1, npix, nwav);
  double **pars = var2dim<double>(pars1, npix, npar);
  float **ratio = var2dim<float>(ratio1, npix, nwav);


  

  float ma = 0.0;
  size_t ndat = (size_t)npix * (size_t)nwav;
  for(size_t ii = 0; ii<ndat; ii++) ma = std::max(ma, dat1[ii]);
  
  //
  // Init MPFIT struct and loop
  //
  int status;
  mp_result result;
  memset(&result,0,sizeof(result));
  //
  mp_par fitpars[npar];
  memset(&fitpars[0], 0, sizeof(fitpars));
  fitpars[0].limited[0] = 1;
  fitpars[0].limited[1] = 1;
  fitpars[0].limits[0] = 0;
  fitpars[0].limits[1] = 3.0*ma;
  fitpars[1].limited[0] = 1;
  fitpars[1].limited[1] = 1;
  fitpars[1].limits[0] = -0.4;
  fitpars[1].limits[1] = 0.4;
  for(int ii=0;ii<=npar-1;++ii) fitpars[ii].side = 0;
  //
  fgd_t tmp;
  tmp.xl = xl;
  tmp.yl = yl;
  tmp.wav = wav;
  tmp.nmean = nmean;
  //
  double ntot = 100. / (npix - 1.0);
  //
  int tid = 0;
  bool master = false;
  //
  int nskip = 0;
  int oper = -1, per=0, ix=  0, ww=0;
  //
#pragma omp parallel default(shared) firstprivate(ww,ix,tid,master,per,oper,status,result,tmp) num_threads(nt)  
  {
    tid = omp_get_thread_num();
    if(tid == 0) {
      master = true;
      int ntp = omp_get_num_threads();
      fprintf(stderr,"cfitgain : nthreads -> %d\n", ntp);
    }
#pragma omp for schedule(dynamic, 1)
    for(ix = 0;ix<npix;++ix){
      //
      // Avoid masked pixels (pixel = 0.0)
      //
      tmp.idat = dat[ix];
      double sum = 0.0;
      for(ww=0;ww<nwav;++ww) sum += tmp.idat[ww];
      //
      if(sum >= 1.E-1){
	//
	// Fit pixel
	//
	status = mpfit(fitgain_model, nwav, npar, pars[ix], fitpars, 0, (void*) &tmp, &result);
	//
	// Get ratio dat / model
	//
	get_ratio(nwav, npar, nmean, wav, dat[ix], pars[ix], xl, yl, ratio[ix]);
      } else{nskip += 1;}
      //
      // Progress counter
      // 
      per = (int)(ntot * ix);
      if(oper != per && master){
	oper = per;
	fprintf(stderr, "\rcfitgain : fitting data -> %d%s", per, "%");
      }
    }
  }
  fprintf(stderr, " \n");
  fprintf(stderr, "cfitgain : %d (%g %s) skipped pixels\n", nskip, (float)(nskip)/npix * 100.0, "%");

    //
    // Clean-up
    //
    delete [] dat;
    delete [] pars;
    delete [] ratio;
    //

}

int fitgain_model2(int nwav, int npar, double *pars, double *dev, double **derivs, void *tmp1){

  fpi_t *fpip = (fpi_t*) tmp1;
  

  //
  // If reflectivity changes, recompute CRISP profile and 
  // re-apply ratio. Also init Bezier control points
  //
  bool ref = false;
  if(pars[2] != fpip->orh){
    dual_fpi(fpip,pars[2]); // CRISP profile
    //
    //for(int ii=0;ii<fpip->npad/2+1;ii++) fpip->otf[ii] *= fpip->ft[ii];//Convolution
    for(int ii=0;ii<fpip->npad/2+1;ii++) {
      double tmp = fpip->otf[ii].re;
      fpip->otf[ii].re = (fpip->ft[ii].re * fpip->otf[ii].re) - (fpip->ft[ii].im * fpip->otf[ii].im);
      fpip->otf[ii].im = (fpip->ft[ii].re * fpip->otf[ii].im) + (fpip->ft[ii].im * tmp);
    }

    //
    fftw_execute(fpip->bplan);// FFT back
    hermite_control(fpip->npad, fpip->xlp, fpip->ylp, fpip->lp1, fpip->lp2, fpip->fp1, fpip->fp2);//Compute interpolation coeff again
    ref = true;
  }

  //
  // Shift line profile (interpolating)
  //
  if((fpip->ech != pars[1]) || ref){
    double *wav1 = new double [fpip->nwav];
    //
    for(int ww=0;ww<fpip->nwav;ww++) wav1[ww] = fpip->wav[ww] - pars[1];
    //
    hermite(fpip->npad, fpip->xlp, fpip->ylp, fpip->nwav, wav1, fpip->imean, 
	    fpip->lp1, fpip->lp2, fpip->fp1, fpip->fp2);
    delete [] wav1;
    fpip->ech = pars[1];
  }

  
  //
  // Linear component and scale factor
  //
  for(int ww=0;ww<fpip->nwav;ww++){
    dev[ww] = (pars[0] * get_polcomp2(nwav, npar, fpip->wav[ww], pars) * fpip->imean[ww]) - fpip->idat[ww];
  }

  return 0;
}



void fitgain2(int nwav, int nmean, int npar, int npix, float *xl, float *yl, float *wav, float *dat1, double *pars1, float *ratio1, int nt, int line){

  bool master = false;
  int tid = 0;
  int ix = 0;
  int status=0;
  int ww=0, per, oper;
  double sum = 0;
  mp_result result;
  mp_par fitpars[npar];

  
  
  float ma = 0.0;
  size_t ndat = (size_t)npix * (size_t)nwav;
  for(size_t ii = 0; ii<ndat; ii++) ma = std::max(ma, dat1[ii]);
  
  
  memset(&fitpars[0], 0, sizeof(fitpars));

  
  fprintf(stderr, "cfitgain2 : nwav=%d, npar=%d, nmean=%d, npix=%d \n", nwav, npar, nmean, npix);


  //fitpars[0].side = 3;
  //
 
  //
  int nskip = 0;
  double *ipar = NULL;


  //
  //#pragma omp parallel default(shared) firstprivate(ww,ix,tid,master,per,oper,status,result,tmp) num_threads(nt)  
  fpi_t fpi[nt];
  //shared(nskip, fpi, nwav, nmean, npar, npix, xl, yl, wav, dat1, pars1, ratio1, nt, line, stderr, fitpars)
#pragma omp parallel default(shared) firstprivate(master, per, oper) private(status, result, ipar, sum, ww, ix, tid) num_threads(nt)  
  {
    memset(&result,0,sizeof(result));
    double ntot = 100. / (npix - 1.0);
    per = 0;
    oper = -1;

    tid = omp_get_thread_num();

    if(tid == 0) {
      master = true;
      int ntp = omp_get_num_threads();
      fprintf(stderr,"cfitgain2 : nthreads -> %d\n", ntp);
      
      for(int kk=0;kk<ntp;kk++){
	fpi[kk].tid = kk;
	init_fpi(fpi[kk], line);
	init_fftw(fpi[kk],nmean,xl,yl);
	dual_fpi(&fpi[kk],fpi[kk].rhr);
      }
 
      //
      // Init MPFIT struct and loop
      //
      fitpars[0].limited[0] = 1;
      fitpars[0].limited[1] = 1;
      fitpars[0].limits[0] = 0;
      fitpars[0].limits[1] = 3.0*ma;
      fitpars[1].limited[0] = 1;
      fitpars[1].limited[1] = 1;
      fitpars[1].limits[0] = -0.4;
      fitpars[1].limits[1] = 0.4;
      if(npar > 2){
	fitpars[2].limits[0] = fpi[0].rhr -0.1;
	fitpars[2].limits[1] = 0.99;
	fitpars[2].limited[0] = 1;
	fitpars[2].limited[1] = 1;
      }
      for(int ii=0;ii<=npar-1;++ii) fitpars[ii].side = 0;
      
    }
#pragma omp barrier
    
    //fprintf(stderr,"cfitgain2 : thread -> %d, hello!\n", tid);
    
    
    
  //
  fpi[tid].wav = wav;
  fpi[tid].xl = xl;
  fpi[tid].yl = yl;
  fpi[tid].nmean = nmean;
  fpi[tid].nwav = nwav;

  //
  fpi[tid].imean = new double [nwav];
  fpi[tid].orh = fpi[tid].rhr - 0.001; // Just set it to any value so it is recalculated in the loop;
  fpi[tid].ech = -1.0;
  
  //
  // Deconvolve spectrum (in reality we multiply it with the ratio of the two 
  // CRISP profiles Tcrisp/Tcrisp_0).
  //
  for(int ii=0;ii<fpi[tid].npad/2+1;ii++) {
    //fpi[tid].ft[ii] /= fpi[tid].otf[ii];
    double r = fpi[tid].otf[ii].im*fpi[tid].otf[ii].im + fpi[tid].otf[ii].re*fpi[tid].otf[ii].re;
    double tmp = fpi[tid].ft[ii].re;
    fpi[tid].ft[ii].re = (fpi[tid].ft[ii].re * fpi[tid].otf[ii].re + fpi[tid].ft[ii].im * fpi[tid].otf[ii].im) / r;
    fpi[tid].ft[ii].im = (fpi[tid].ft[ii].im * fpi[tid].otf[ii].re - tmp * fpi[tid].otf[ii].im) / r;
  }




#pragma omp for schedule(dynamic, 2) reduction(+:nskip)
    for(ix = 0;ix<npix;++ix){
      //
      // Avoid masked pixels (pixel = 0.0)
      //
      fpi[tid].idat = &dat1[ix*nwav];
      ipar = &pars1[ix*npar];
      sum = 0.0;
      for(ww=0;ww<nwav;ww++) sum += fpi[tid].idat[ww];
      //
      if(sum >= 1.){
	//
	// Fit pixel
	//
	if(ipar[2] < 0.5) ipar[2] = fpi[tid].rhr;

	status = mpfit(fitgain_model2, nwav, npar, ipar, fitpars, 0, (void*) &fpi[tid], &result);

	//
	// Get ratio dat / model
	//
	get_ratio2(nwav, npar, nmean, wav, fpi[tid].idat, ipar, xl, yl, &ratio1[ix*nwav]);
	ipar = NULL;
      } else{nskip += 1;}
      //
      // Progress counter
      // 
      per = (int)(ntot * ix);
      if(oper != per && master){
	oper = per;
	fprintf(stderr, "\rcfitgain2 : fitting data -> %d%s", per, "%");
      }
    }

    //
    // Clean-up
    //
    if(master){
      for(int ii=0;ii<nt;ii++) clean_fpi(fpi[ii]);
    }
#pragma omp barrier
  }
  fprintf(stderr, "\rcfitgain2 : fitting data -> %d%s", (int)100, "%");
  fprintf(stderr, "cfitgain2 : %d (%g %s) skipped pixels\n", nskip, (float)(nskip)/npix * 100.0, "%");
}

void sim_polcal(double *p, pol_t *pol){
  //
  char inam[] = "sim_polcal :";
  //
  //pol_t *pol = (pol_t *) polv;
  //
  double cd = cos(p[16] * dtor);
  double sd = sin(p[16] * dtor);
  //
  // Invert matrix
  //
  double *mi = invert4x4<double>(p);

  //
  // Normalize inverse matrix
  //
  double scl = 0.0;
  for(uint8_t jj = 0; jj<pol->nlc;jj++) scl += mi[jj];

  scl = 1.0 / scl;
  for(uint8_t ii = 0; ii<16;ii++) mi[ii] *= scl;
  //
  // compute observed curves
  //
  double Q, U, a, ca, sa, A, B, C, tmp, C0, C1, C2, C3, ca2, sa2;
  //
  int Nsum = 0, ele = 0;
  //


  for(int o=0;o<pol->nlp;++o){
    //
    Q = cos(dtor2 * pol->lp[o]);
    U = sin(dtor2 * pol->lp[o]);
    //
    for(int r=0;r<pol->nqwp;++r){
      //
      a = 2.0 * (pol->qwp[r] + p[17]) * dtor;
      ca = cos(a);
      sa = sin(a);
      //
      A = ca * ca + sa * sa * cd;
      B = sa * sa + ca * ca * cd;
      C = ca * sa * (1.0 - cd);
      //
      scl = mi[0] * pol->dat[ele];
      for(int ii=1;ii<pol->nlc;++ii) scl += mi[ii] * pol->dat[ele+ii];
      //
      scl = fabs(scl);
      ele += 4;
      //
      C1 = scl * (A * Q + C * U);
      C2 = scl * (C * Q + B * U);
      C3 = scl * (sa  * Q - ca * U) * sd;
      //
      for(int l=0;l<pol->nlc;++l){
	pol->cur[Nsum] = (p[4*l ] * scl + 
			  p[4*l+1] * C1 +
			  p[4*l+2] * C2 +
			  p[4*l+3] * C3);
	Nsum += 1;
      }
      
    }
  }
  pol->ncall += 1;
  delete [] mi;
}
//
int compute_chisq_lm(int ndata, int npar, double *p, double *dev, double **derivs, void *polv){
  //
  pol_t *pol = (pol_t *) polv;
  //
  // Get simulated polcal curves
  //
  sim_polcal(p, pol);
  //
  // Compute model - obs
  //
  pol->Chisq = 0.0;
  //
  for(int ii=0;ii<ndata;ii++) {
    dev[ii] = pol->cur[ii] -  pol->dat[ii];
    pol->Chisq += dev[ii] * dev[ii];
  }
  //
  pol->Chisq /= (double)(pol->ndata);
  //
  return 0;
}

void polcal_2d(pol_t pol, int nthreads, float *dat2d, double *res,float *cs, float *qwp, float *lp){
  //
  //float **bla = var2dim<float>(dat2d, pol.npix, pol.ndata);
  //
  char inam[] = "polcal_2d :";
  fprintf(stderr, "%s ndata=%d, npix=%d, nlc=%d, nqwp=%d, nlp=%d, nele=%d\n",inam,pol.ndata, pol.npix, pol.nlc, pol.nqwp, pol.nlp,pol.npix*pol.ndata);
  //
  int i, tid;
  int npix = pol.npix, ndata = pol.ndata;
  //
  // Prepare L-M settings
  //
  int status, chunk = 10;
  int per=0, oper=-1;
  float ntot = 100.0 / (pol.npix - 1.0);
  //
  //shared(pol, bla,nthreads,chunk, ndone, pars, npar, ndata, npix)
#pragma omp parallel shared(cs, res, dat2d, qwp, lp, npix, ndata, inam, ntot, nthreads, per, oper, chunk) private(i, tid, status) firstprivate(pol) num_threads(nthreads)
  {
    int npar = 18;
    mp_par pars[npar];
    memset(&pars[0], 0, sizeof(pars));
    // Variable limits
    for(i=0;i<16;++i){
      pars[i].limited[0] = 1;
      pars[i].limited[1] = 1;
      pars[i].limits[0] = -1.5;
      pars[i].limits[1] = 1.5;
    }
    // Fix the retardance and offset angle
    pars[16].fixed = 1;
    pars[17].fixed = 1;
    //
    for(i=0;i<npar-2;++i){ 
      pars[i].step = 1.0E-10;
      //pars[ii].side = 0;
    }
    mp_result result;
    memset(&result,0,sizeof(result));       /* Zero results structure */
    //
    mp_config config;
    memset(&config,0,sizeof(config));
    config.maxiter = 40;
    //
    tid = omp_get_thread_num();
    pol.tid = tid;
    if (tid == 0){
	fprintf(stderr,"%s %d threads\n", inam, nthreads);
      }
    pol.cur = new double [ndata];;
    pol.qwp = qwp;
    pol.lp = lp;
    //
#pragma omp for schedule(dynamic, 1)
      for(i=0;i<npix;i++){
	//
	//
	pol.dat = &dat2d[i*ndata];
	//
	status = mpfit(compute_chisq_lm, ndata, npar, 
		       &res[i*18], pars, &config, (void*) &pol, &result);
	//
	cs[i] = pol.Chisq;     
	//
	if(pol.tid==0){
	  per = (int)(ntot*i);
	  if(per != oper){
	    oper = per;
	    fprintf(stderr, "\r%s fitting data -> %d%s",inam, per, "%");
	  }

	}
      }
    
  delete [] pol.cur;
  pol.lp=NULL;
  pol.qwp=NULL;
  }
  fprintf(stderr,"\r%s fitting data -> %d%s", inam, 100, "%\n");
  //
  // Clean-up
  //
  //delete [] res;
}


namespace {

    const int cacheSize = 128;
    const float deltasqr = 4;
    const float beta = 4;

    float distanceCache[cacheSize*cacheSize];

    bool calculateCache( void ) {
        double i2,j2;
        for( int i=0; i<cacheSize; ++i ) {
            i2 = i*i;
            for( int j=0; j<cacheSize; ++j ) {
                j2 = j*j;
                distanceCache[i*cacheSize + j] = pow( sqrt( i2+j2+deltasqr ), -beta );
            }
        }

        return true;
    }

    float inverseDistanceWeight( float **a, int xl, int xh, int yl, int yh, int xb, int yb ) {

        static bool dummy = calculateCache();

        float sum=0.0, weighted_values=0.0, tmp;
        int xabs, yabs;
        for( int y=yl; y<=yh; ++y ) {
            yabs = abs(y-yb);
            for( int x=xl; x<=xh; ++x ) {
                xabs = abs(x-xb);
                if( a[y][x] ) {
                    tmp = distanceCache[yabs*cacheSize + xabs];
                    weighted_values += tmp*a[y][x];
                    sum += tmp;
                }
            }
        }

        return weighted_values/sum;
    }

}

void fillpix( int nx, int ny, float *img, uint8_t *mask, int nt ) {

    int npix = nx * ny;

    float **img2d = var2dim<float>( img, ny, nx );
    uint8_t **mask2d = var2dim<uint8_t>( mask, ny, nx );

    int nbad = nx * ny;
    for( int yy=0; yy<ny; yy++ ) for( int xx=0; xx<nx; xx++ ) nbad -= mask2d[yy][xx];
    // fprintf( stderr, "cfillpix2 : There are %d bad pixels\n", nbad );
    if( nbad<=0 ) return;

    int *idx = new int[nbad];
    int *idy = new int[nbad];
    float *bp = new float[nbad];
    memset( idx, 0, nbad*sizeof( int ) );
    memset( idy, 0, nbad*sizeof( int ) );
    memset( bp, 0, nbad*sizeof( float ) );

    
    int k = 0;
    for( int yy=0; yy<ny; yy++ ) {
        for( int xx=0; xx<nx; xx++ ) {
            if( mask2d[yy][xx] == 0 ) {
                idx[k] = xx;
                idy[k] = yy;
                k++;
            }
        }
    }

    int dl = cacheSize / 2;
    float ntot = nbad>1 ? 100. / (nbad-1.0): 100.0;
    int ch=1;
    int tid=0, per=0, oper=-1;
    uint8_t master=0;
    k = 0;
    #pragma omp parallel default(shared) firstprivate(k, tid, master, per, oper) num_threads(nt)
    {
        tid = omp_get_thread_num();
        master = 0;
        if( tid == 0 ) {
          //  fprintf( stderr, "cfillpix2 : nthreads -> %d\n", nt );
            master = 1;
        }

        #pragma omp for schedule(dynamic, ch)
        for( k=0; k<nbad; k++ ) {
            bp[k] = inverseDistanceWeight( img2d, std::max( idx[k]-dl, 0 ), std::min( idx[k]+dl, nx-1 ),
                                                  std::max( idy[k]-dl, 0 ), std::min( idy[k]+dl, ny-1 ), idx[k], idy[k] );

            if( master ) {
                per = ntot*k;
                if( per > oper ) {
                  //  fprintf( stderr, "\rcfillpix2 : %d %s", per, "%" );
                    oper = per;
                }
            }
        }

    }//end parallel block

    // fprintf( stderr, "\rcfillpix : %d %s\n", 100, "%" );

    for( int ii=0; ii<nbad; ii++ ) {
        img2d[idy[ii]][idx[ii]] = bp[ii];
    }

    delete[] bp;
    delete[] idx;
    delete[] idy;
    delete[] img2d;
    delete[] mask2d;

    return;
}




