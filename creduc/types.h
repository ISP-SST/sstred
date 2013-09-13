#ifndef __TYPES_H__  // __TYPES_H__
#define __TYPES_H__
typedef float float32_t;
typedef double float64_t;
typedef float fp_t;

struct mcomplex{
  double re;
  double im;
};
typedef struct mcomplex complex_t;

struct fgd{
  float32_t *xl;
  float32_t *yl;
  float32_t *idat;
  float32_t *wav;
  int32_t nmean;
  complex_t *oftm;
  complex_t *otft;
  float64_t *psf;
  int32_t npsf;
};

struct pol{
  int32_t nqwp, nlp, nlc, ndata, npix, nthread, tid;
  float64_t *res, Chisq, *cur;
  float32_t *dat, *qwp, *lp, *cs, *dat2d;
  int32_t ncall;
};
typedef struct fgd fgd_t;
typedef struct pol pol_t;
#endif               // __TYPES_H__
