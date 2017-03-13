#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <sys/types.h>
#include "mymath.h"
#include "export.h"

//# ifdef __cplusplus
extern "C" {
    //#endif // __cplusplus

    typedef unsigned char byte;
    static void free_array(unsigned char *arg)
    {
        // fprintf(stderr,"freeing memory!\n");
        free(arg);
    }

    IDL_VPTR cbezier2(int argc,IDL_VPTR argv[]) {
        IDL_VPTR x = argv[0];
        IDL_VPTR y = argv[1];
        IDL_VPTR xp_in = argv[2];

        int n =      x->value.s.arr->dim[0];
        int np = xp_in->value.s.arr->dim[0];
        double *res;

        //
        // Check input data types to be float or double(s)
        //
        if(x->type < 4 || y->type < 4 || xp_in->type < 4 ||
                x->type > 5 || y->type > 5 || xp_in->type > 5) {

            fprintf(stderr,"cbezier2: ERROR, input data must be double or float!\n");
            short *dum = (short*)malloc(sizeof(short));
            IDL_MEMINT dims[]= {1};
            IDL_VPTR result=IDL_ImportArray(1,dims,IDL_TYP_INT,(byte*)dum,free_array,0);
            return result;

        }

        //
        // Check xp type
        //
        double *xp;
        if(xp_in->type == 4) {
            xp = (double*)malloc(np*sizeof(double));
            for(int k=0; k<np; k++) xp[k] = ((float*)xp_in->value.s.arr->data)[k];
        } else {
            xp = (double*)xp_in->value.s.arr->data;
        }

        //
        // Check individual types
        //
        if((x->type + y->type) == 10 || (x->type + y->type) == 8) {
            int type = x->type;
            res = bezier2(n, (byte*)x->value.s.arr->data,(byte*)y->value.s.arr->data, np, xp, type);
        } else {
            double *xx,*yy;
            if(x->type == 5) {
                xx = (double*)x->value.s.arr->data;
            } else {
                xx = (double*)malloc(n*sizeof(double));
                for(int k=0; k<n; k++) xx[k] = ((float*)x->value.s.arr->data)[k];
            }
            if(y->type == 5) {
                yy = (double*)y->value.s.arr->data;
            } else {
                yy = (double*)malloc(n*sizeof(double));
                for(int k=0; k<n; k++) yy[k] = ((float*)y->value.s.arr->data)[k];
            }

            int type = 5;
            res = bezier2(n, (byte*)xx,(byte*)yy, np, xp, type);


            if(x->type != 5) free(xx);
            if(y->type != 5) free(yy);
        }

        IDL_MEMINT dims[]= {np};
        IDL_VPTR result=IDL_ImportArray(1,dims,IDL_TYP_DOUBLE, (byte*)res,free_array,0);


        if(xp_in->type != 5) free(xp);
        return result;
    }


    IDL_VPTR cbezier3(int argc,IDL_VPTR argv[]) {
        IDL_VPTR x = argv[0];
        IDL_VPTR y = argv[1];
        IDL_VPTR xp_in = argv[2];

        int n =      x->value.s.arr->dim[0];
        int np = xp_in->value.s.arr->dim[0];
        double *res;

        //
        // Check input data types to be float or double(s)
        //
        if(x->type < 4 || y->type < 4 || xp_in->type < 4 ||
                x->type > 5 || y->type > 5 || xp_in->type > 5) {

            fprintf(stderr,"cbezier3: ERROR, input data must be double or float!\n");
            short *dum = (short*)malloc(sizeof(short));
            IDL_MEMINT dims[]= {1};
            IDL_VPTR result=IDL_ImportArray(1,dims,IDL_TYP_INT,(byte*)dum,free_array,0);
            return result;

        }

        //
        // Check xp type
        //
        double *xp;
        if(xp_in->type == 4) {
            xp = (double*)malloc(np*sizeof(double));
            for(int k=0; k<np; k++) xp[k] = ((float*)xp_in->value.s.arr->data)[k];
        } else {
            xp = (double*)xp_in->value.s.arr->data;
        }

        //
        // Check individual types
        //
        if((x->type + y->type) == 10 || (x->type + y->type) == 8) {
            int type = x->type;
            res = bezier3(n, (byte*)x->value.s.arr->data,(byte*)y->value.s.arr->data, np, xp, type);
        } else {
            double *xx,*yy;
            if(x->type == 5) {
                xx = (double*)x->value.s.arr->data;
            } else {
                xx = (double*)malloc(n*sizeof(double));
                for(int k=0; k<n; k++) xx[k] = ((float*)x->value.s.arr->data)[k];
            }
            if(y->type == 5) {
                yy = (double*)y->value.s.arr->data;
            } else {
                yy = (double*)malloc(n*sizeof(double));
                for(int k=0; k<n; k++) yy[k] = ((float*)y->value.s.arr->data)[k];
            }

            int type = 5;
            res = bezier3(n, (byte*)xx,(byte*)yy, np, xp, type);


            if(x->type != 5) free(xx);
            if(y->type != 5) free(yy);
        }

        IDL_MEMINT dims[]= {np};
        IDL_VPTR result=IDL_ImportArray(1,dims,IDL_TYP_DOUBLE, (byte*)res,free_array,0);


        if(xp_in->type != 5) free(xp);
        return result;
    }


    int IDL_Load(void)
    {
        static IDL_SYSFUN_DEF2 function_addr[] = {{ cbezier2, "CBEZIER2", 3, 3, 0, 0 },{ cbezier3, "CBEZIER3", 3, 3, 0, 0 }
        };

        return IDL_SysRtnAdd(function_addr, TRUE, IDL_CARRAY_ELTS(function_addr));
    }
    //# ifdef __cplusplus
};
//#endif // __cplusplus
