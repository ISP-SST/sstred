#ifndef LGRDMATCH
#define LGRDMATCH

void ana_stretch(const int ny, const int nx,  double *im, int npy, int npx, const double *gr,  double *out);
void ana_stretch_full_matrix(const int ny, const int nx,int npy, int npx, const double *gr,  double *outx, double *outy);
void ana_gridmatch(int const ny, int const nx, double* __restrict__ p1, double* __restrict__ p2,
		   int const nyg, int const nxg, int* __restrict__ gy, int* __restrict__ gx, int const dy, int const dx,
		   double const gwid, int stretch_clip, double* __restrict__ out);

#endif
