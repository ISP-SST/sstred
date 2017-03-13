#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include "types.h"
#include "creduc.h"
//
// creduc.cc wrapper routines in C
//
// Jaime de la Cruz Rodriguez (IFA-UU)

extern "C" {
  short cdescatter(int argc,void *argv[]){
    int nx = ARG_INT(argv, 0);
    int ny = ARG_INT(argv, 1);
    int nthreads = ARG_INT(argv, 6);
    int verbose = ARG_INT(argv, 7);
    float *timg = ARG_FLOAT_ARRAY(argv,2);
    float *tfgain = ARG_FLOAT_ARRAY(argv,3);
    float *tfpsf = ARG_FLOAT_ARRAY(argv,4);
    float *res = ARG_FLOAT_ARRAY(argv,5);
    //
    descatter(nx, ny, timg, tfgain, tfpsf, res, nthreads, verbose);
    //
    return 0;
  }
  short cconvolve(int argc, void *argv[]){
    //
    int nx = ARG_INT(argv, 0);
    int ny = ARG_INT(argv, 1);
    int nx1 = ARG_INT(argv, 2);
    int ny1 = ARG_INT(argv, 3);
    int nthreads = ARG_INT(argv, 7);
    int verbose = ARG_INT(argv, 8);
    //
    float *timg = ARG_FLOAT_ARRAY(argv,4);
    float *tfpsf = ARG_FLOAT_ARRAY(argv,5);
    float *res1 = ARG_FLOAT_ARRAY(argv,6);
    //
    convolve(nx, ny, nx1, ny1, timg, tfpsf, res1, nthreads, verbose);
    //
    return 0;
  }
  short cfitgain(int argc, void *argv[]){
    //
    int nwav = ARG_INT(argv, 0);
    int nmean = ARG_INT(argv, 1);
    int npar = ARG_INT(argv, 2);
    int npix = ARG_INT(argv, 3);
    int nt = ARG_INT(argv, 10);

    //
    float *xl = ARG_FLOAT_ARRAY(argv,4);
    float *yl = ARG_FLOAT_ARRAY(argv,5);
    float *wav = ARG_FLOAT_ARRAY(argv,6);
    float *dat1 = ARG_FLOAT_ARRAY(argv,7);
    //
    double *pars1 = ARG_FLOAT64_ARRAY(argv,8); // At input contains the guess parameters to init L-M
    float *ratio1 = ARG_FLOAT_ARRAY(argv,9); // At input contains the guess parameters to init L-M

    //
    fitgain(nwav, nmean, npar, npix, xl, yl, wav, dat1, pars1, ratio1, nt); 
    //
    return 0;
  }
  short cfitgain2(int argc, void *argv[]){
    //
    int nwav = ARG_INT(argv, 0);
    int nmean = ARG_INT(argv, 1);
    int npar = ARG_INT(argv, 2);
    int npix = ARG_INT(argv, 3);
    int nt = ARG_INT(argv, 10);
    int line = ARG_INT(argv, 11);

    //
    float *xl = ARG_FLOAT_ARRAY(argv,4);
    float *yl = ARG_FLOAT_ARRAY(argv,5);
    float *wav = ARG_FLOAT_ARRAY(argv,6);
    float *dat1 = ARG_FLOAT_ARRAY(argv,7);
    //
    double *pars1 = ARG_FLOAT64_ARRAY(argv,8); // At input contains the guess parameters to init L-M
    float *ratio1 = ARG_FLOAT_ARRAY(argv,9);

    //
    fitgain2(nwav, nmean, npar, npix, xl, yl, wav, dat1, pars1, ratio1, nt, line); 
    //
    return 0;
  }

  /*
  short cpolcal_1d(int argc, void *argv[]){
    //
    pol_t pol;
    //
    pol.nlc = ARG_INT(argv, 0);
    pol.nqwp = ARG_INT(argv, 1);
    pol.nlp = ARG_INT(argv, 2);
    //
    pol.dat = ARG_FLOAT_ARRAY(argv,3);
    pol.qwp =  ARG_FLOAT_ARRAY(argv,4);
    pol.lp =  ARG_FLOAT_ARRAY(argv,5);
    //
    pol.res = ARG_FLOAT64_ARRAY(argv,6);
    pol.cur = ARG_FLOAT_ARRAY(argv,7);
    pol.ndata = pol.nlc * pol.nqwp * pol.nlp;
    //
    polcal_1d(pol);
    //
    return 0;
  }
  
  short cpolcal_fov(int argc, void *argv[]){
    //
    pol_t pol;
    //
    pol.nlc = ARG_INT(argv, 0);
    pol.nqwp = ARG_INT(argv, 1);
    pol.nlp = ARG_INT(argv, 2);
    pol.npix = ARG_INT(argv, 3);
    //
    pol.dat2d = ARG_FLOAT_ARRAY(argv,4);
    pol.qwp =  ARG_FLOAT_ARRAY(argv,5);
    pol.lp =  ARG_FLOAT_ARRAY(argv,6);
    //
    pol.res =  ARG_FLOAT64_ARRAY(argv,7);
    pol.cs =  ARG_FLOAT_ARRAY(argv,8);
    //
    pol.ndata = pol.nlc * pol.nqwp * pol.nlp;
    //
    polcal_fov(pol);
    //
  }
  */
  short cpolcal_2d(int argc, void *argv[]){
    //
    pol_t pol;
    //
    pol.nlc = ARG_INT(argv, 0);
    pol.nqwp = ARG_INT(argv, 1);
    pol.nlp = ARG_INT(argv, 2);
    pol.npix = ARG_INT(argv, 3);
    //
    float *dat2d = ARG_FLOAT_ARRAY(argv,4);
    float *qwp =  ARG_FLOAT_ARRAY(argv,5);
    float *lp =  ARG_FLOAT_ARRAY(argv,6);
    //
    double *res =  ARG_FLOAT64_ARRAY(argv,7);
    float *cs =  ARG_FLOAT_ARRAY(argv,8);
    //
    int nthreads = ARG_INT(argv, 9);
    //
    pol.ndata = pol.nlc * pol.nqwp * pol.nlp;
    //
    polcal_2d(pol, nthreads, dat2d, res, cs, qwp, lp);
    //
    return 0;
  }
  
  short cfillpix(int argc, void *argv[]){
    //
    int nx = ARG_INT(argv, 0);
    int ny = ARG_INT(argv, 1);
    // float *img = ARG_FLOAT_ARRAY(argv,2);
    // uint8_t *mask = ARG_BYTE_ARRAY(argv,3);
    int nt = ARG_INT(argv, 4);
    //
    //fillpix(nx, ny, img, mask, nt);

    fillpix(nx, ny, (float*)(argv[2]), (uint8_t*)(argv[3]), nt);
    //
    return 0;
  }

};


