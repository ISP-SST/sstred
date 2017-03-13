#include "mmath.h"
#include <stdlib.h>
#include <stdint.h>

double *bezier2(int n, unsigned char *x, unsigned char *y, int np, double *xp, int type) {
    double *res;

    switch(type) {
    case 4:
        res=ccbezier2<float>(n, (float*)x, (float*)y , np, xp);
        break;
    case 5:
        res=ccbezier2<double>(n, (double*)x, (double*)y , np, xp);
        break;
    }

    return res;
}
double *bezier3(int n, unsigned char *x, unsigned char *y, int np, double *xp, int type) {
    double *res;

    switch(type) {
    case 4:
        res=ccbezier3<float>(n, (float*)x, (float*)y , np, xp);
        break;
    case 5:
        res=ccbezier3<double>(n, (double*)x, (double*)y , np, xp);
        break;
    }

    return res;
}
