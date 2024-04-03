#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "sep_mask.h"

int32_t find_connected(uint8_t *mask, int32_t sx, int32_t sy, int32_t *idx,
		    int32_t i0)
{
    int32_t nx1=0,nx2=0,j,l,k;
    int32_t sr[]={sx, -1, 1, -sx, sx-1, sx+1, -sx-1, -sx+1};

    mask[i0] = 2;
    idx[nx1++] = i0;
    // now check the surrounding
    do {
	nx2 = 0;
	for(l=0;l<nx1;l++) {
	    for(j=0;j<8;j++) {
		if(mask[(k=idx[l]+sr[j])] == 1) {
		    mask[k] = 2;
		    idx[nx1++] = k;
		    nx2++;
		}
	    }
	}
    } while(nx2 > 0);
    return(nx1);
}

short sep_mask(int argc, void *argv[])
{
    int32_t i=0;
    uint8_t *mask    = ARG_BYTE_ARRAY(argv, i++);
    int32_t *newmask = ARG_INT_ARRAY(argv, i++);
    int32_t sx       = ARG_INT(argv, i++);
    int32_t sy       = ARG_INT(argv, i++);
    int32_t alim     = ARG_INT(argv, i++);
    int32_t ner      = ARG_INT(argv, i++);

    int32_t j,k,l,ntot,nmax=0,*idx,nx1,nx2,nfound=0;
    int32_t sr[]={sx, -1, 1, -sx, sx-1, sx+1, -sx-1, -sx+1};

    ntot = sx*sy;
    for(i=sx;i<ntot-sx;i++) {
	if(mask[i] != 0) nmax++;
    }
    idx = (int32_t *) malloc(sizeof(int32_t)*nmax);

    if(!ner) {   // Ignore chunks conected to borders
	for(i=sx+1;i<(2*sx-1);i++) {              // Bottom
	    if(mask[i] == 0) continue;
	    nx1 = find_connected(mask,sx,sy,idx,i);
	    for(j=0;j<nx1;j++) mask[idx[j]] = 0;
	}
	for(i=(ntot-2*sx+1);i<(ntot-sx-1);i++) {  // Top
	    if(mask[i] == 0) continue;
	    nx1 = find_connected(mask,sx,sy,idx,i);
	    for(j=0;j<nx1;j++) mask[idx[j]] = 0;
	}
	for(i=sx+1;i<(ntot-sx+1);i += sx) {       // Left
	    if(mask[i] == 0) continue;
	    nx1 = find_connected(mask,sx,sy,idx,i);
	    for(j=0;j<nx1;j++) mask[idx[j]] = 0;
	}
	for(i=2*sx-2;i<ntot-sx-1;i += sx) {       // Right
	    if(mask[i] == 0) continue;
	    nx1 = find_connected(mask,sx,sy,idx,i);
	    for(j=0;j<nx1;j++) mask[idx[j]] = 0;
	}
    }
	
    for(i=sx+1;i<ntot-sx-1;i++) {
	if(mask[i] == 0) continue;
	nx1 = find_connected(mask,sx,sy,idx,i);
	if(nx1 >= alim) {
	    nfound++;
	    for(j=0;j<nx1;j++) newmask[idx[j]] = nfound;
	}
	for(j=0;j<nx1;j++) mask[idx[j]] = 0;
    }
			     
    free(idx);
    return(0);
}
