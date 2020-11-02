#ifndef IMTOOLS
#define IMTOOLS

void nearest2D(int const ny, int const nx, double* __restrict__ y, double* __restrict__ x, double*  d_in, int const ny1, int const nx1, double*  yy_in, double*  xx_in, double*  res_in, int const nthreads, double const missing);

void bilint2D(int const ny, int const nx, double* __restrict__ y, double* __restrict__ x, double*  d_in, int const ny1, int const nx1, double*  yy_in, double*  xx_in, double*  res_in, int const nthreads, double const missing);

#endif
