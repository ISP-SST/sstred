; docformat = 'rst'

;+
; Smooth an inverse modulation matrix map with PSFs from the image
; restoration process.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
; 
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    lc : 
;   
;       The momfbd data for the different LC states.
;   
;    imodmat : 
;   
;       The inverse modulation matrix.
;   
; 
; :Keywords:
; 
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2014-01-10 : PS. Cleanup (remove unused code).
;
;-
function red_matrix2momfbd, lc, imodmat

  ;; Name of this function
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  Nstokes = 4
  Nlc = n_elements(lc)
  
  chan = 0
  Nx = lc[0].clip[chan,0,1] - lc[0].clip[chan,0,0] + 1
  Ny = lc[0].clip[chan,1,1] - lc[0].clip[chan,1,0] + 1
  print,' '
  print, inam +' : Clipping modulation matrix to '+red_stri(Nx)+' x '+red_stri(Ny)

  dim=size(imodmat, /dimension)

  mny = dim[2] - 1
  mnx = dim[1] - 1

  imm = replicate({patch:lc[0].patch, clip:lc[0].clip},[4,4])
  dim = size(lc[0].patch, /dimension)

  for ix = 0, dim[0]-1 do begin
    for iy = 0, dim[1]-1 do begin
      
      pnx = lc[0].patch[ix,iy].xh - lc[0].patch[ix,iy].xl + 1
      pny = lc[0].patch[ix,iy].yh - lc[0].patch[ix,iy].yl + 1
      cpmm = dindgen(4,4,pnx,pny)      
      pmm = dindgen(pnx,pny)

      for ilc = 0, Nlc-1 do begin
      
        clip = red_getclips(lc[ilc], ix, iy)
        tpmm=imodmat[*,(clip[0]>0)<mnx:(clip[1]<mnx),(clip[2]>0):(clip[3]<mny)]

        for istokes=0,Nstokes-1 do begin
          pmm[*,*]=median(tpmm[ilc*4 + istokes,*,*])                       ; fill extra bits with the median
          pmm[(-clip[0])>0,(-clip[2])>0]=reform(tpmm[ilc*4 + istokes,*,*]) ; copy what we do have
          dip = size(pmm, /dim)
          res = fltarr(dip)
          imm[istokes,ilc].patch[ix,iy].img = red_convolve(red_fillnan(pmm) $
                                                           , reform(lc[ilc].patch[ix,iy].psf[*,*,0]))
        endfor                  ; istokes

      endfor                    ; ilc
      
;      clip = red_getclips(lc1, ix,iy )
;      tpmm=imodmat[*,(clip[0]>0)<mnx:(clip[1]<mnx),(clip[2]>0):(clip[3]<mny)]
;
;      for istokes=0,Nstokes-1 do begin
;        pmm[*,*]=median(tpmm[4+istokes,*,*])                       ; fill extra bits with the mean...
;        pmm[(-clip[0])>0,(-clip[2])>0]=reform(tpmm[4+istokes,*,*]) ; copy what we do have
;
;        dip = size(pmm, /dim)
;        res = fltarr(dip)
;        imm[istokes,1].patch[ix,iy].img = red_convolve(red_fillnan(pmm) $
;                                                       , reform(lc1.patch[ix,iy].psf[*,*,0]))
;      endfor                    ; istokes
;
;      clip = red_getclips(lc2, ix, iy)
;      tpmm=imodmat[*,(clip[0]>0<mnx):(clip[1]<mnx),(clip[2]>0):(clip[3]<mny)]
;
;      for istokes=0,Nstokes-1 do begin
;        pmm[*,*]=median(tpmm[8+istokes,*,*])                       ; fill extra bits with the mean...
;        pmm[(-clip[0])>0,(-clip[2])>0]=reform(tpmm[8+istokes,*,*]) ; copy what we do have
;
;        dip = size(pmm, /dim)
;        res = fltarr(dip)
;        imm[istokes,2].patch[ix,iy].img = red_convolve(red_fillnan(pmm) $
;                                                       , reform(lc2.patch[ix,iy].psf[*,*,0]))
;      endfor                    ; istokes
;
;      clip = red_getclips(lc3, ix, iy)
;      if not[clip[1]-clip[0]+1 eq pnx] or not[clip[3]-clip[2]+1 eq pny] then stop
;      tpmm=imodmat[*,(clip[0]>0<mnx):(clip[1]<mnx),(clip[2]>0):(clip[3]<mny)]
;
;      for istokes=0,Nstokes-1 do begin
;        pmm[*,*]=median(tpmm[12+istokes,*,*])                       ; fill extra bits with the mean...
;        pmm[(-clip[0])>0,(-clip[2])>0]=reform(tpmm[12+istokes,*,*]) ; copy what we do have
;
;        dip = size(pmm, /dim)
;        res = fltarr(dip)
;        imm[istokes,3].patch[ix,iy].img = red_convolve(red_fillnan(pmm) $
;                                                       , reform(lc3.patch[ix,iy].psf[*,*,0]))
;      endfor                    ; istokes

    endfor                      ; iy
  endfor                        ; ix

  return, imm
  
end
