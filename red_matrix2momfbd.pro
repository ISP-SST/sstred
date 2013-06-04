function red_matrix2momfbd, lc0, lc1, lc2, lc3, imodmat
                                ;
  inam = 'red_matrix2momfbd : '
  chan = 0
  nx = lc0.clip[chan,0,1] - lc0.clip[chan,0,0] + 1
  ny = lc0.clip[chan,1,1] - lc0.clip[chan,1,0] + 1
  print,' '
  print,inam +'Clipping modulation matrix to '+red_stri(nx)+' x '+red_stri(ny)

                                ;
  dim=size(imodmat, /dimension)
                                ;
  mny = dim[2] - 1
  mnx = dim[1] - 1
                                ;
  imm = replicate({patch:lc0.patch, clip:lc0.clip},[4,4])
  dim = size(lc0.patch, /dimension)
                                ;
  for x = 0, dim[0]-1 do begin
     for y = 0, dim[1]-1 do begin
                                ;
        pnx = lc0.patch[x,y].xh - lc0.patch[x,y].xl + 1
        pny = lc0.patch[x,y].yh - lc0.patch[x,y].yl + 1
        cpmm = dindgen(4,4,pnx,pny)      
        pmm = dindgen(pnx,pny)
                                ;
        clip = red_getclips(lc0, x, y)
                                ;
        tpmm=imodmat[*,(clip[0]>0)<mnx:(clip[1]<mnx),(clip[2]>0):(clip[3]<mny)]
        for l=0,3 do begin
           pmm[*,*]=median(tpmm[l,*,*])                       ; fill extra bits with the mean...
           pmm[(-clip[0])>0,(-clip[2])>0]=reform(tpmm[l,*,*]) ; copy what we do have
                                ; print, (-clip[0])>0 , (-clip[2])>0 
                                ;imm[l,0].patch[x,y].img =
                                ;reform(red_convl2d(fillnan(pmm),lc0.patch[x,y].psf[*,*,0]))
                                ;
           dip = size(pmm, /dim)
           res = fltarr(dip)
                                ;f77_convlim, dip[0], dip[1], fillnan(pmm), dip[0], dip[1], reform(lc0.patch[x,y].psf[*,*,0]), res
           
           imm[l,0].patch[x,y].img = red_convolve(red_fillnan(pmm), reform(lc0.patch[x,y].psf[*,*,0]))
                                ; imm[l,0].patch[x,y].img = reform(red_convl2d(fillnan(pmm),lc0.patch[x,y].psf[*,*,0]))
        end
                                ;
        clip = red_getclips(lc1, x,y )
                                ;
                                ; if not[clip[1]-clip[0]+1 eq pnx] or not[clip[3]-clip[2]+1 eq pny] then stop
        tpmm=imodmat[*,(clip[0]>0)<mnx:(clip[1]<mnx),(clip[2]>0):(clip[3]<mny)]
        for l=0,3 do begin
           pmm[*,*]=median(tpmm[4+l,*,*])                       ; fill extra bits with the mean...
           pmm[(-clip[0])>0,(-clip[2])>0]=reform(tpmm[4+l,*,*]) ; copy what we do have
                                ;
           dip = size(pmm, /dim)
           res = fltarr(dip)
                                ;f77_convlim, dip[0], dip[1], fillnan(pmm), dip[0], dip[1], reform(lc1.patch[x,y].psf[*,*,0]), res
           imm[l,1].patch[x,y].img = red_convolve(red_fillnan(pmm), reform(lc1.patch[x,y].psf[*,*,0]))
                                ;
                                ;imm[l,1].patch[x,y].img = reform(red_convl2d(fillnan(pmm),lc1.patch[x,y].psf[*,*,0]))
        end
                                ;
        clip = red_getclips(lc2, x, y)
                                ;
                                ;  if not[clip[1]-clip[0]+1 eq pnx] or not[clip[3]-clip[2]+1 eq pny] then stop
        tpmm=imodmat[*,(clip[0]>0<mnx):(clip[1]<mnx),(clip[2]>0):(clip[3]<mny)]
        for l=0,3 do begin
           pmm[*,*]=median(tpmm[8+l,*,*])                       ; fill extra bits with the mean...
           pmm[(-clip[0])>0,(-clip[2])>0]=reform(tpmm[8+l,*,*]) ; copy what we do have
                                ;
           dip = size(pmm, /dim)
           res = fltarr(dip)
                                ;f77_convlim, dip[0], dip[1], fillnan(pmm), dip[0], dip[1], reform(lc2.patch[x,y].psf[*,*,0]), res
           imm[l,2].patch[x,y].img = red_convolve(red_fillnan(pmm), reform(lc2.patch[x,y].psf[*,*,0]))
                                ;
                                ;imm[l,2].patch[x,y].img = reform(red_convl2d(fillnan(pmm),lc2.patch[x,y].psf[*,*,0]))
        end
;
        clip = red_getclips(lc3, x, y)
                                ;
        if not[clip[1]-clip[0]+1 eq pnx] or not[clip[3]-clip[2]+1 eq pny] then stop
                                ;
        tpmm=imodmat[*,(clip[0]>0<mnx):(clip[1]<mnx),(clip[2]>0):(clip[3]<mny)]
                                ;
        for l=0,3 do begin
           pmm[*,*]=median(tpmm[12+l,*,*])                       ; fill extra bits with the mean...
           pmm[(-clip[0])>0,(-clip[2])>0]=reform(tpmm[12+l,*,*]) ; copy what we do have
                                ;
           dip = size(pmm, /dim)
           res = fltarr(dip)
                                ;f77_convlim, dip[0], dip[1], fillnan(pmm), dip[0], dip[1], reform(lc3.patch[x,y].psf[*,*,0]), res
           imm[l,3].patch[x,y].img = red_convolve(red_fillnan(pmm), reform(lc3.patch[x,y].psf[*,*,0]))
                                ;
                                ;imm[l,3].patch[x,y].img = reform(red_convl2d(fillnan(pmm),reform(lc3.patch[x,y].psf[*,*,0])))
        end
     end
  end
                                ;
  return, imm
end
