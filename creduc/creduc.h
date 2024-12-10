
//
#define ARG_BYTE(av, i) (*((uint8_t *)((av)[i])))
#define ARG_SHORT(av, i) (*((short *)((av)[i])))
#define ARG_INT(av, i) (*((int *)((av)[i])))
#define ARG_FLOAT(av, i) (*((float *)((av)[i])))
#define ARG_FLOAT64(av, i) (*((double *)((av)[i])))
//
#define ARG_BYTE_ARRAY(av, i) ((uint8_t *)((av)[i]))
#define ARG_SHORT_ARRAY(av, i) ((short *)((av)[i]))
#define ARG_INT_ARRAY(av, i) ((int *)((av)[i]))
#define ARG_FLOAT_ARRAY(av, i) ((float *)((av)[i]))
#define ARG_FLOAT64_ARRAY(av, i) ((double *)((av)[i]))
//
// Functions
//
extern void descatter(int nx, int ny, float *img, float *fgain, float *fpsf, float *res, int nthreads, int verbose);
//
extern void convolve(int inx, int iny, int pnx, int pny, float *timg, float *tfpsf, float *res1,int nthreads, int verbose);
extern void fitgain(int nwav, int nmean, int nlc, int npar, int npix, float *xl, float *yl, float *wav, float *dat1, double *pars1, float *ratio1, double* sig, int nt);
extern void fitgain2(int nwav, int nmean, int npar, int npix, float *xl, float *yl, float *wav, float *dat1, double *pars1, float *ratio1, int nt, int line);
extern void polcal_1d(pol_t &pol);
extern void polcal_fov(pol_t &pol);
extern void polcal_2d(pol_t pol, int nthreads, float *dat2d, double *res,float *cs, float *qwp, float *lp);
extern void polcal_2d_mask(pol_t pol, int nthreads, float *dat2d, double *res,float *cs, float *qwp, float *lp, uint8_t *mask);
extern void fillpix(int nx, int ny, float *img,  uint8_t *mask, int nt);
