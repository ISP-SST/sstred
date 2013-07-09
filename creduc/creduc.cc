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
#include "types.h"
#include "fftw3.h"
#include "mymath.h"
#include "mmem.h"
#include "mpfit.h"
#include "physc.h"

/*
#include "Genetic.h"
#include "Tools.h"
using namespace ROYAC;
*/
float64_t sqr(float64_t var){
  return var * var;
}
float64_t sqr(float64_t &var){
  return var * var;
}

inline int32_t ind(int32_t x, int32_t y, int32_t nx){
  return (y*nx + x);
}
//
// descatter
//
void descatter(int nx, int ny, float32_t *timg, float32_t *tfgain, float32_t *tfpsf, float32_t *res, int32_t nthreads, int32_t verbose){
  //
  // Init variables
  //
  char inam[] = "descatter :";
  float32_t **img = var2dim<float32_t>(timg,ny,nx);
  float32_t **fgain = var2dim<float32_t>(tfgain,ny,nx);
  float32_t **fpsf = var2dim<float32_t>(tfpsf,ny,nx);
  //
  float64_t **psf=mat2d<float64_t>(2*ny, 2*nx);
  memset(psf[0], 0, 4 * nx * ny * sizeof(float64_t));
  float64_t **gain = mat2d<float64_t>(ny, nx);
  //
  int32_t nx2 = nx/2;
  int32_t ny2 = ny/2;
  float64_t suma = 0.0;
  for(int32_t y=0;y<ny;++y) for(int32_t x=0;x<nx;++x) gain[y][x]=fgain[y][x];
  for(int32_t y=0;y<ny;++y) for(int32_t x=0;x<nx;++x) suma+=fpsf[y][x];
  for(int32_t y=0;y<ny;++y) for(int32_t x=0;x<nx;++x) psf[y+ny2][x+nx2]=fpsf[y][x]/suma;
  psf_reorder<float64_t>(psf, 2*ny);
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
  float64_t **im = mat2d<float64_t>(2*ny, 2*nx);
  float64_t **cim = mat2d<float64_t>(2*ny, 2*nx);
  memset(im[0],0,4*nx*ny*sizeof(float64_t));
  complex_t **ft=mat2d<complex_t>(2*ny,nx+1);
  //
  // FFTW Plans
  //
  fftw_plan fplan=fftw_plan_dft_r2c_2d(2*ny,2*nx,im[0],(fftw_complex*)(ft[0]),FFTW_ESTIMATE);
  fftw_plan bplan=fftw_plan_dft_c2r_2d(2*ny,2*nx,(fftw_complex*)(ft[0]),cim[0],FFTW_ESTIMATE);
  //
  float64_t np = nx * ny;
  float64_t np4 = 4.0 * np;
  for(int32_t y=0;y<ny;++y) for(int32_t x=0;x<nx;++x) {
      im[y+ny2][x+nx2]=img[y][x]/(1.0+gain[y][x] * gain[y][x]);
    }
  
  //
  // Iterate solution
  //
  float64_t cs = 1E11, ocs;
  int32_t i = 0, iter_max = 50;
  //
  do{
    memcpy(cim[0],im[0],4*nx*ny*sizeof(float64_t));
    for(int32_t y=0;y<ny;++y) for(int32_t x=0;x<nx;++x) cim[y+ny2][x+nx2]*=gain[y][x];
    //
    // Forward transform
    //
    fftw_execute_dft_r2c(fplan,cim[0],(fftw_complex*)(ft[0]));
    //
    // convolve with OTF
    //
    for(int32_t y=0;y<2*ny;++y) 
      for(int32_t x=0;x<nx+1;++x){
	float64_t tmp=ft[y][x].re;
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
    for(int32_t y=0;y<2*ny;++y) for(int32_t x=0;x<2*nx;++x) cim[y][x]/=np4;
    //
    for(int32_t y=0;y<ny;++y) for(int32_t x=0;x<nx;++x) cim[y+ny2][x+nx2]*=gain[y][x];
    //
    // Compute Chisq
    //
    cs=0.0;
    for(int32_t y=0;y<ny;++y) 
      for(int32_t x=0;x<nx;++x){
	float64_t bla=img[y][x]-cim[y+ny2][x+nx2];
        cs+=sqr(im[y+ny2][x+nx2]-bla);
        im[y+ny2][x+nx2]=bla;
      }
    cs/= (float64_t)(np);
    if(verbose) fprintf(stderr, "%s iter = %d, ChiSq = %e \n", inam, i,cs);
  }while((cs>1E-10)&&(++i<iter_max));
  //
  // Destroy FFTW plans
  //
  fftw_destroy_plan(bplan);
  fftw_destroy_plan(fplan);
  fftw_destroy_plan(otfplan);
  //
  //  Write result to res array
  //
  for(int32_t y=0;y<ny;++y) for(int32_t x=0;x<nx;++x) res[(y*nx+x)] = im[y+ny2][x+nx2];
  //
  // Clean-up
  //
  //del_mat2d<float32_t>(img);
  //del_mat2d<float32_t>(fgain);
  //del_mat2d<float32_t>(fpsf);
  del_mat2d<float64_t>(im);
  del_mat2d<float64_t>(cim);
  del_mat2d<float64_t>(psf);
  del_mat2d<float64_t>(gain);
  del_mat2d<complex_t>(otf);
  del_mat2d<complex_t>(ft);
  //
  delete [] img;
  delete [] fpsf;
  delete [] fgain;
  //
}
void convolve(int32_t inx, int32_t iny, int32_t pnx, int32_t pny, float32_t *img, float32_t *fpsf, float32_t *res1, int32_t nthreads, int32_t verbose){
  char inam[] = "convolve :";
  //
  // reshape 2D arrays
  //
  //float32_t **img = var2dim<float32_t>(timg,iny,inx);
  //float32_t **fpsf = var2dim<float32_t>(tfpsf,pny,pnx);
  //
  // Allocate arrays with extra-space for padding
  //
  int32_t npadx = inx + pnx;
  int32_t npady = iny + pny;
  //
  if(verbose){
    fprintf(stderr, "%s Image -> nx=%d, ny=%d\n",inam,inx, iny);
    fprintf(stderr, "%s PSF -> nx=%d, ny=%d\n",inam,pnx, pny);
    fprintf(stderr, "%s FFT -> npad_x=%d, npad_y=%d, nthreads=%d\n",inam, npadx, npady, nthreads);
  }
  //
  float64_t **pimg = mat2d<float64_t>(npady, npadx);
  float64_t **ppsf= mat2d<float64_t>(npady, npadx);
  memset(ppsf[0], 0, npadx * npady * sizeof(float64_t));
  //
  // Pad with the mean of the image
  //
  float64_t suma = 0.0;
  for(int32_t y=0;y<iny;++y) for(int32_t x=0;x<inx;++x) suma += img[ind(x,y,inx)];
  suma/= (float64_t)(inx*iny);
  for(int32_t y=0;y<npady;++y) for(int32_t x=0;x<npadx;++x) pimg[y][x] = suma;
  //
  // Copy image onto padded array
  //
  for(int32_t y=0;y<iny;++y) for(int32_t x=0;x<inx;++x) pimg[y][x] = img[ind(x,y,inx)];
  //
  // Copy PSF to new array
  //
  int32_t pny2 = pny/2;
  int32_t pnx2 = pnx/2;
  //
  suma = 0;
  for(int32_t y=0;y<pny;++y) for(int32_t x=0;x<pnx;++x) suma += fpsf[ind(x,y,pnx)];

  for(int32_t y=0;y<pny2;++y) for(int32_t x=0;x<pnx2;++x) ppsf[y][x] = fpsf[ind(x+pnx2,y+pny2,pnx)]/suma; // first
  //
  int32_t x0 = npadx - pnx2;
  int32_t y0 = npady - pny2;
  for(int32_t y=0;y<pny2;++y) for(int32_t x=0;x<pnx2;++x) ppsf[y+y0][x+x0] = fpsf[ind(x,y,pnx)]/suma;  // Fourth
  for(int32_t y=0;y<pny2;++y) for(int32_t x=0;x<pnx2;++x) ppsf[y][x+x0] = fpsf[ind(x,y+pny2, pnx)]/suma; // third
  for(int32_t y=0;y<pny2;++y) for(int32_t x=0;x<pnx2;++x) ppsf[y+y0][x] = fpsf[ind(x+pnx2,y,pnx)]/suma; // second
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
  for(int32_t y=0;y<npady;++y) for(int32_t x=0;x<npadx/2+1;++x){
      float64_t tmp=ft[y][x].re;
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
  float64_t np4 = npadx * npady;
  for(int32_t y=0;y<iny;++y) for(int32_t x=0;x<inx;++x) res1[ind(x,y,inx)] = pimg[y][x] / np4;
  //
  // Clean-up
  //
  del_mat2d<float64_t>(pimg);
  del_mat2d<float64_t>(ppsf);
  del_mat2d<complex_t>(otf);
  del_mat2d<complex_t>(ft);
  //
  if(verbose) fprintf(stderr, "done \n");
  //
}


float64_t get_polcomp(int32_t &nwav, int32_t &npar, float32_t wav, float64_t *&pars){
  //
  float64_t res = 1.0;
  float64_t twav = 1.0;
  //
  if(npar <= 2) return res;
  //
  int32_t ii = 2;
  while(ii<npar){
    twav *= wav;
    res += pars[ii] * twav;
    ii += 1;
  }
  //
  return res;
}
int32_t fitgain_model(int32_t nwav, int32_t npar, float64_t *pars, float64_t *dev, float64_t **derivs, void *tmp1){
  //
  fgd_t *tmp = (fgd_t*) tmp1;
  //
  // interpolate to shifted grid (cavity error)
  // 
  float64_t *twav = new float64_t [nwav];
  for(int32_t ii = 0;ii<nwav;++ii) twav[ii] = tmp->wav[ii] - pars[1];
  //
  float64_t *spec = intep<float64_t>(tmp->nmean, tmp->xl, tmp->yl, nwav, twav);
  //
  // Compute difference between model and observations
  //
  for(int32_t ii=0;ii<nwav;++ii) dev[ii] = (pars[0] * spec[ii] * get_polcomp(nwav, npar, tmp->wav[ii], pars)) - tmp->idat[ii]; 
  //
  // Clean-up
  //
  delete [] spec;
  delete [] twav;
  //
  return 0;
}
void get_ratio(int32_t &nwav, int32_t &npar, int32_t &nmean, float32_t *wav, float32_t *dat, float64_t *pars, float32_t *xl, float32_t *yl, float32_t *ratio){
  float64_t *twav = new float64_t [nwav];
  for(int32_t ii = 0;ii<nwav;++ii) twav[ii] = wav[ii] - pars[1];
  //
  float64_t *spec = intep<float64_t>(nmean, xl, yl, nwav, twav);
  //
  // Compute difference between model and observations
  //
  for(int32_t ii=0;ii<nwav;++ii) ratio[ii] = dat[ii] / (pars[0] * spec[ii] * get_polcomp(nwav, npar, wav[ii], pars));
  //
  delete [] twav;
  delete [] spec;
}
//
void fitgain(int32_t nwav, int32_t nmean, int32_t npar, int32_t npix, float32_t *xl, float32_t *yl, float32_t *wav, float32_t *dat1, float64_t *pars1, float32_t *ratio1){
  fprintf(stderr, "cfitgain : nwav=%d, npar=%d, nmean=%d, npix=%d \n", nwav, npar, nmean, npix);
  // return;
  //
  // Reshape vars
  //
  float32_t **dat = var2dim<float32_t>(dat1, npix, nwav);
  float64_t **pars = var2dim<float64_t>(pars1, npix, npar);
  float32_t **ratio = var2dim<float32_t>(ratio1, npix, nwav);
  //
  // Init MPFIT struct and loop
  //
  int32_t status;
  mp_result result;
  memset(&result,0,sizeof(result));
  //
  mp_par fitpars[npar];
  memset(&fitpars[0], 0, sizeof(fitpars));
  fitpars[0].limited[0] = 1;
  fitpars[0].limited[1] = 1;
  fitpars[0].limits[0] = 0;
  fitpars[0].limits[1] = 4096;
  fitpars[1].limited[0] = 1;
  fitpars[1].limited[1] = 1;
  fitpars[1].limits[0] = -0.1;
  fitpars[1].limits[1] = 0.1;
  for(int32_t ii=0;ii<=npar-1;++ii) fitpars[ii].side = 0;
  //
  fgd_t tmp;
  tmp.xl = xl;
  tmp.yl = yl;
  tmp.wav = wav;
  tmp.nmean = nmean;
  //
  float64_t ntot = 100. / (npix - 1.0);
  //
  int32_t nskip = 0;
  int32_t oper = -1, per=0;
  for(int32_t ix = 0;ix<npix;++ix){
    //
    // Avoid masked pixels (pixel = 0.0)
    //
    tmp.idat = dat[ix];
    float64_t sum = 0.0;
    for(int32_t ww=0;ww<nwav;++ww) sum += tmp.idat[ww];
    //
    if(sum >= 1.E-3){
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
    per = (int32_t)(ntot * ix);
    if(oper != per){
      oper = per;
      fprintf(stderr, "\rcfitgain : fitting data -> %d%s", per, "%");
    }
  }
  fprintf(stderr, " \n");
  fprintf(stderr, "cfitgain : %d (%g %s) skipped pixels\n", nskip, (float32_t)(nskip)/npix * 100.0, "%");
  //
  // Clean-up
  //
  delete [] dat;
  delete [] pars;
  delete [] ratio;
  //
}

void sim_polcal(float64_t *p, pol_t *pol){
  //
  char inam[] = "sim_polcal :";
  //
  //pol_t *pol = (pol_t *) polv;
  //
  float64_t cd = cos(p[16] * dtor);
  float64_t sd = sin(p[16] * dtor);
  //
  // Invert matrix
  //
  float64_t *mi = invert4x4<float64_t>(p);

  //
  // Normalize inverse matrix
  //
  float64_t scl = 0.0;
  for(uint8_t jj = 0; jj<pol->nlc;jj++) scl += mi[jj];

  scl = 1.0 / scl;
  for(uint8_t ii = 0; ii<16;ii++) mi[ii] *= scl;
  //
  // compute observed curves
  //
  float64_t Q, U, a, ca, sa, A, B, C, tmp, C0, C1, C2, C3, ca2, sa2;
  //
  int32_t Nsum = 0, ele = 0;
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
int32_t compute_chisq_lm(int32_t ndata, int32_t npar, float64_t *p, float64_t *dev, float64_t **derivs, void *polv){
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
  for(int32_t ii=0;ii<ndata;ii++) {
    dev[ii] = pol->cur[ii] -  pol->dat[ii];
    pol->Chisq += dev[ii] * dev[ii];
  }
  //
  pol->Chisq /= (float64_t)(pol->ndata);
  //
  return 0;
}

void polcal_2d(pol_t pol, int32_t nthreads, float32_t *dat2d, float64_t *res,float32_t *cs, float32_t *qwp, float32_t *lp){
  //
  //float32_t **bla = var2dim<float32_t>(dat2d, pol.npix, pol.ndata);
  //
  char inam[] = "polcal_2d :";
  fprintf(stderr, "%s ndata=%d, npix=%d, nlc=%d, nqwp=%d, nlp=%d, nele=%d\n",inam,pol.ndata, pol.npix, pol.nlc, pol.nqwp, pol.nlp,pol.npix*pol.ndata);
  //
  int32_t i, tid;
  int32_t npix = pol.npix, ndata = pol.ndata;
  //
  // Prepare L-M settings
  //
  int32_t status, chunk = 10;
  int32_t per=0, oper=-1;
  float32_t ntot = 100.0 / (pol.npix - 1.0);
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
      pars[i].limits[0] = -1.2;
      pars[i].limits[1] = 1.2;
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
    config.maxiter = 15;
    //
    tid = omp_get_thread_num();
    pol.tid = tid;
    if (tid == 0){
	fprintf(stderr,"%s %d threads\n", inam, nthreads);
      }
    pol.cur = new float64_t [ndata];;
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
	  per = (int32_t)(ntot*i);
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


#define beta 4.0
#define deltasqr 4.0

template <class T> T inv_dist_wght(T **a,int xl,int xh,int yl,int yh,int xb,int yb)
{
  T weight=0.0,res=0.0;
  for(int y=yl;y<=yh;++y)
    for(int x=xl;x<=xh;++x)
      if(a[y][x]){
        T c = pow(sqrt((T)sqr(x-xb)+(T)sqr(y-yb)+deltasqr),-beta);
        res+= c * a[y][x];
        weight += c;
      }
  return(res/weight);
}

#undef beta
#undef deltasqr

void fillpix(int nx, int ny, float *img, uint8_t *mask, int nt){
  int maxsize = 128;
  int npix = nx*ny;

  float **img2d = var2dim<float>(img, ny, nx);
  uint8_t **mask2d = var2dim<uint8_t>(mask, ny, nx);

  int nbad = nx*ny;
  for(int yy=0;yy<ny;yy++) for(int xx=0;xx<nx;xx++) nbad -= mask2d[yy][xx];
  fprintf(stderr,"cfillpix : There are %d bad pixels\n", nbad);
  if(nbad <= 0) return;  
  
  int *idx = new int [nbad];
  int *idy = new int [nbad];
  float *bp = new float [nbad];
  memset(idx,0, nbad*sizeof(int));
  memset(idy,0, nbad*sizeof(int));
  memset(bp,0, nbad*sizeof(float));

  int k = 0;
  for(int yy=0;yy<ny;yy++){
    for(int xx=0;xx<nx;xx++){
      if(mask2d[yy][xx] == 0){
	idx[k] = xx;
	idy[k] = yy;
	k++;
      }
    }
  }
  
  int dl = maxsize / 2;
  float ntot = 100. / (nbad - 1.0);
  int ch = 1;
  int tid=0,per=0,oper=-1;
  uint8_t master=0;
  k = 0;
#pragma omp parallel default(shared) firstprivate(k, tid, master, per, oper) num_threads(nt)
  {
    tid = omp_get_thread_num();
    master = 0;
    if(tid == 0){
      fprintf(stderr,"cfillpix : nthreads -> %d\n",nt);
      master = 1;
    }
    
#pragma omp for schedule(dynamic, ch)
    for(k = 0; k<nbad; k++){
      bp[k] = inv_dist_wght<float>(img2d,std::max(idx[k]-dl,0),std::min(idx[k]+dl,nx-1),std::max(idy[k]-dl,0),std::min(idy[k]+dl,ny-1), idx[k],idy[k]);
      
      if(master) {
	per = ntot * k;
	if(per > oper){
	  fprintf(stderr,"\rcfillpix : %d %s",per,"%");
	  oper = per;
	}
      }
    }
    
  }//end parallel block
  
  fprintf(stderr,"\rcfillpix : %d %s\n",100,"%");
  
  for(int ii = 0; ii<nbad; ii++){
    img2d[idy[ii]][idx[ii]] = bp[ii];
  }

  delete [] bp;
  delete [] idx;
  delete [] idy;
  delete [] img2d;
  delete [] mask2d;
  //
  return;
}
