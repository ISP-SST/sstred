#include <iostream>
#include <complex.h>
#include "fftw3.h"
#include "types.h"

struct cfpi{
  double erh, erl, ech, ecl, orh;
  double rhr, rlr, slr, shr, dhr, dlr, w0, dw; // FPI vars
  double *calp;
  int line, npad, nmean, nwav;
  fftw_plan otfplan, bplan, fplan;
  double complex *otf, *ft,*otfm;
  double *xlp, *ylp, *imean;
  double *tr, *tw;
  float *xl, *yl, *wav, *idat;
  double *c1,*c2;
  bool scl;
};
typedef struct cfpi fpi_t;


template <class T> T sqr(T val){
  return val * val;
}
template <class T> T sqr2(T val){
  return val * val * val * val;
}

void init_fpi(fpi_t &fpi, int line);
void clean_fpi(fpi_t &fpi);
void init_fftw(fpi_t &fpi, int nmean, float *xl, float *yl);
void dual_fpi(fpi_t *fpi, double erh);
complex_t cmul(complex_t a, complex_t b);
complex_t cdiv(complex_t a, complex_t b);
void bezier3_control(int n, double *x, double *y, double *&c1, double *&c2);
void bezier3(int n, double *x, double *y, int np, double *xp, double *yp, double *c1, double *c2);
