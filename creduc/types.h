typedef float float32_t;
typedef double float64_t;
typedef float fp_t;

struct complex{
  float64_t re;
  float64_t im;
};
struct fgd{
  float32_t *xl;
  float32_t *yl;
  float32_t *idat;
  float32_t *wav;
  int32_t nmean;
};

struct pol{
  int32_t nqwp, nlp, nlc, ndata, npix, nthread, tid;
  float64_t *res, Chisq, *cur;
  float32_t *dat, *qwp, *lp, *cs, *dat2d;
  int32_t ncall;
};
typedef struct complex complex_t;
typedef struct fgd fgd_t;
typedef struct pol pol_t;
