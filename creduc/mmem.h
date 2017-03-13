#include <sys/types.h>
#include <stdint.h>
#include <algorithm>
/*
  Memory allocation functions -> templates overview
 */
template <class T> T **mat2d(int32_t nx1, int32_t nx2);
template <class T> T ***mat3d(int32_t nx1, int32_t nx2, int32_t nx3);
template <class T> T ****mat4d(int32_t nx1, int32_t nx2, int32_t nx3, int32_t nx4);
template <class T> T **var2dim(T *data, int32_t nx1, int32_t nx2);
template <class T> T ***var3dim(T *data, int32_t nx1, int32_t nx2, int32_t nx3);
template <class T> T ****var4dim(T *data, int32_t nx1, int32_t nx2, int32_t nx3, int32_t nx4);
//
template <class T> void del_mat2d(T **p);
template <class T> void del_mat3d(T ***p);
template <class T> void del_mat4d(T ****p);
//

//
/*
  Definitions
 */
template <class T> T **mat2d(int32_t nx1, int32_t nx2) {
    T **p;
    p = new T* [nx1];
    p[0] = new T [nx1 * nx2];
    for(int32_t x1=1; x1<=nx1-1; ++x1) p[x1] = p[x1-1] + nx2;
    return p;
}
//
template <class T> T ***mat3d(int32_t nx1, int32_t nx2, int32_t nx3) {
    //
    T ***p;
    p=new T** [nx1];
    p[0]=new T* [nx1*nx2];
    p[0][0]=new T [nx1*nx2*nx3];
    for(int x2=1; x2<=(nx2-1); ++x2) p[0][x2]=p[0][x2-1]+nx3;
    for(int x1=1; x1<=(nx1-1); ++x1) {
        p[x1]=p[x1-1]+nx2;
        p[x1][0]=p[x1-1][0]+nx2*nx3;
        for(int x2=1; x2<=(nx2-1); ++x2) p[x1][x2]=p[x1][x2-1]+nx3;
    }
    return p;
}
//
template <class T> T ***mat3d2(int32_t nx1, int32_t nx2, int32_t nx3) {

}
//
template <class T> T ****mat4d(int32_t nx1, int32_t nx2, int32_t nx3, int32_t nx4) {
    T ****p;
    p = new T*** [nx1];
    p[0] = new T** [nx1 * nx2];
    p[0][0] = new T* [nx1 * nx2 * nx3];
    p[0][0][0] = new T [nx1 * nx2 * nx3 * nx4];
    //
    for(int32_t x3=1; x3<=nx3-1; ++x3) p[0][0][x3]=p[0][0][x3-1] + nx4;
    for(int32_t x2=1; x2<=nx2-1; ++x2) {
        p[0][x2] = p[0][x2-1] + nx3;
        p[0][x2][0] = p[0][x2-1][0] + (nx3 * nx4);
        for(int32_t x3=1; x3<=nx3-1; ++x3) p[0][x2][x3] = p[0][x2][x3-1] + nx4;
    }
    for(int32_t x1=1; x1<=nx1-1; ++x1) {
        p[x1] = p[x1-1] + nx2;
        p[x1][0]=p[x1-1][0]+nx2*nx3;
        p[x1][0][0]=p[x1-1][0][0] + (nx2 * nx3 * nx4);
        for(int x3=1; x3<=nx3-1; ++x3) p[x1][0][x3]=p[x1][0][x3-1] + nx4;
        for(int x2=1; x2<=nx2-1; ++x2) {
            p[x1][x2]=p[x1][x2-1] + nx3;
            p[x1][x2][0]=p[x1][x2-1][0] + (nx3*nx4);
            for(int x3=1; x3<=nx3-1; ++x3) p[x1][x2][x3]=p[x1][x2][x3-1] + nx4;
        }
    }
    return p;
}
//
template <class T> T **var2dim(T *data, int32_t nx1, int32_t nx2) {
    T **p;
    p=new T* [nx1];
    p[0] = data;
    for(int32_t x1=1; x1<=nx1-1; ++x1) p[x1]=p[x1-1] + nx2;
    return p;
}
//
template <class T> T ***var3dim(T *data, int32_t nx1, int32_t nx2, int32_t nx3) {
    T ***p;
    p=new T** [nx1];
    p[0]=new T* [nx1 * nx2];
    p[0][0]=data;
    for(int32_t x2=1; x2<=nx2-1; ++x2) p[0][x2]=p[0][x2-1] + nx3;
    for(int32_t x1=1; x1<=nx1-1; ++x1) {
        p[x1]=p[x1-1] + nx2;
        p[x1][0]=p[x1-1][0] + nx2 * nx3;
        for(int32_t x2=1; x2<=nx2-1; ++x2) p[x1][x2]=p[x1][x2-1] + nx3;
    }
    return p;
}
//
template <class T> T ****var4dim(T *data, int32_t nx1, int32_t nx2, int32_t nx3, int32_t nx4) {
    T ****p;
    p=new T*** [nx1];
    p[0]=new T** [nx1*nx2];
    p[0][0]=new T* [nx1*nx2*nx3];
    p[0][0][0]=data;
    for(int x3=1; x3<=nx3-1; ++x3) p[0][0][x3]=p[0][0][x3-1] + nx4;
    for(int x2=1; x2<=nx2-1; ++x2) {
        p[0][x2]=p[0][x2-1] + nx3;
        p[0][x2][0]=p[0][x2-1][0] + (nx3 * nx4);
        for(int x3=1; x3<=nx3-1; ++x3) p[0][x2][x3]=p[0][x2][x3-1] + nx4;
    }
    for(int x1=1; x1<=nx1-1; ++x1) {
        p[x1]=p[x1-1] + nx2;
        p[x1][0]=p[x1-1][0] + (nx2 * nx3);
        p[x1][0][0]=p[x1-1][0][0] + (nx2 * nx3 * nx4);
        for(int x3=1; x3<=nx3-1; ++x3) p[x1][0][x3]=p[x1][0][x3-1] + nx4;
        for(int x2=1; x2<=nx2-1; ++x2) {
            p[x1][x2]=p[x1][x2-1] + nx3;
            p[x1][x2][0]=p[x1][x2-1][0] + (nx3 * nx4);
            for(int x3=1; x3<=nx3-1; ++x3) p[x1][x2][x3]=p[x1][x2][x3-1] + nx4;
        }
    }
    return p;
}
//
template <class T> void del_mat2d(T **p) {
    delete[] (p[0]);
    delete[] (p);
}
//
template <class T> void del_mat3d(T ***p) {
    delete[] (p[0][0]);
    delete[] (p[0]);
    delete[] (p);
}
//
template <class T> void del_mat4d(T ****p) {
    delete[] (p[0][0][0]);
    delete[] (p[0][0]);
    delete[] (p[0]);
    delete[] (p);
}
//
template <class T> void psf_reorder(T **f, int32_t np)
{
    int32_t nh = np/2;
    T *buf=new T [nh];
    for(int32_t x=0; x<nh; ++x) {
        memcpy(buf,f[x],nh*sizeof(T));
        memcpy(f[x],f[x+nh]+nh,nh*sizeof(T));
        memcpy(f[x+nh]+nh,buf,nh*sizeof(T));
        //
        memcpy(buf,f[x+nh],nh*sizeof(T));
        memcpy(f[x+nh],f[x]+nh,nh*sizeof(T));
        memcpy(f[x]+nh,buf,nh*sizeof(T));
    }
    delete[] buf;
}

