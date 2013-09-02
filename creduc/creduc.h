
//
#define ARG_BYTE(av, i) (*((uint8_t *)((av)[i])))
#define ARG_SHORT(av, i) (*((short *)((av)[i])))
#define ARG_INT(av, i) (*((int32_t *)((av)[i])))
#define ARG_FLOAT(av, i) (*((float32_t *)((av)[i])))
#define ARG_FLOAT64(av, i) (*((float64_t *)((av)[i])))
//
#define ARG_BYTE_ARRAY(av, i) ((uint8_t *)((av)[i]))
#define ARG_SHORT_ARRAY(av, i) ((short *)((av)[i]))
#define ARG_INT_ARRAY(av, i) ((int32_t *)((av)[i]))
#define ARG_FLOAT_ARRAY(av, i) ((float32_t *)((av)[i]))
#define ARG_FLOAT64_ARRAY(av, i) ((float64_t *)((av)[i]))
//
// Functions
//
extern void descatter(int nx, int ny, float32_t *img, float32_t *fgain, float32_t *fpsf, float32_t *res, int nthreads, int verbose);
//
extern void convolve(int32_t inx, int32_t iny, int32_t pnx, int32_t pny, float32_t *timg, float32_t *tfpsf, float32_t *res1,int32_t nthreads, int32_t verbose);
extern void fitgain(int32_t nwav, int32_t nmean, int32_t npar, int32_t npix, float32_t *xl, float32_t *yl, float32_t *wav, float32_t *dat1, float64_t *pars1, float32_t *ratio1, int32_t nt);
extern void fitgain2(int32_t nwav, int32_t nmean, int32_t npar, int32_t npix, float32_t *xl, float32_t *yl, float32_t *wav, float32_t *dat1, float64_t *pars1, float32_t *ratio1, int32_t nt);
extern void polcal_1d(pol_t &pol);
extern void polcal_fov(pol_t &pol);
extern void polcal_2d(pol_t pol, int32_t nthreads, float32_t *dat2d, float64_t *res,float32_t *cs, float32_t *qwp, float32_t *lp);
extern void fillpix(int nx, int ny, float *img,  uint8_t *mask, int nt);
