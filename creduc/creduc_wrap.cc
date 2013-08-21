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
    int32_t nx = ARG_INT(argv, 0);
    int32_t ny = ARG_INT(argv, 1);
    int32_t nthreads = ARG_INT(argv, 6);
    int32_t verbose = ARG_INT(argv, 7);
    float32_t *timg = ARG_FLOAT_ARRAY(argv,2);
    float32_t *tfgain = ARG_FLOAT_ARRAY(argv,3);
    float32_t *tfpsf = ARG_FLOAT_ARRAY(argv,4);
    float32_t *res = ARG_FLOAT_ARRAY(argv,5);
    //
    descatter(nx, ny, timg, tfgain, tfpsf, res, nthreads, verbose);
    //
    return 0;
  }
  short cconvolve(int argc, void *argv[]){
    //
    int32_t nx = ARG_INT(argv, 0);
    int32_t ny = ARG_INT(argv, 1);
    int32_t nx1 = ARG_INT(argv, 2);
    int32_t ny1 = ARG_INT(argv, 3);
    int32_t nthreads = ARG_INT(argv, 7);
    int32_t verbose = ARG_INT(argv, 8);
    //
    float32_t *timg = ARG_FLOAT_ARRAY(argv,4);
    float32_t *tfpsf = ARG_FLOAT_ARRAY(argv,5);
    float32_t *res1 = ARG_FLOAT_ARRAY(argv,6);
    //
    convolve(nx, ny, nx1, ny1, timg, tfpsf, res1, nthreads, verbose);
    //
    return 0;
  }
  short cfitgain(int argc, void *argv[]){
    //
    int32_t nwav = ARG_INT(argv, 0);
    int32_t nmean = ARG_INT(argv, 1);
    int32_t npar = ARG_INT(argv, 2);
    int32_t npix = ARG_INT(argv, 3);
    int32_t nt = ARG_INT(argv, 10);

    //
    float32_t *xl = ARG_FLOAT_ARRAY(argv,4);
    float32_t *yl = ARG_FLOAT_ARRAY(argv,5);
    float32_t *wav = ARG_FLOAT_ARRAY(argv,6);
    float32_t *dat1 = ARG_FLOAT_ARRAY(argv,7);
    //
    float64_t *pars1 = ARG_FLOAT64_ARRAY(argv,8); // At input contains the guess parameters to init L-M
    float32_t *ratio1 = ARG_FLOAT_ARRAY(argv,9); // At input contains the guess parameters to init L-M

    //
    fitgain(nwav, nmean, npar, npix, xl, yl, wav, dat1, pars1, ratio1, nt); 
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
    float32_t *dat2d = ARG_FLOAT_ARRAY(argv,4);
    float32_t *qwp =  ARG_FLOAT_ARRAY(argv,5);
    float32_t *lp =  ARG_FLOAT_ARRAY(argv,6);
    //
    float64_t *res =  ARG_FLOAT64_ARRAY(argv,7);
    float32_t *cs =  ARG_FLOAT_ARRAY(argv,8);
    //
    int32_t nthreads = ARG_INT(argv, 9);
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


