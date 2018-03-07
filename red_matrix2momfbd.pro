; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    lc0 : 
;   
;   
;   
;    lc1 : 
;   
;   
;   
;    lc2 : 
;   
;   
;   
;    lc3 : 
;   
;   
;   
;    imodmat : 
;   
;   
;   
; 
; :Keywords:
; 
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2014-01-10 : PS cleanup (remove unused code)
; 
;   2016-08-09 : Re-write to also handle channels with reversed align_clip
;
;   2018-03-07 : THI:  Re-writing the CRISP polarimetry to properly
;                apply clips/offsets to the demodulation matrices.
; 
;-
function red_matrix2momfbd, lc0, lc1, lc2, lc3, imodmat, xoff=xoff, yoff=yoff

  imm = replicate({patch:lc0.patch, clip:lc0.clip},[4,4])

  dim = size(imodmat, /dimension)
  imodmat2 = make_array( dim, type=size(imodmat, /type) )
  do_copy = 1
  if keyword_set(xoff) && keyword_set(yoff) then begin
    xodim = size(xoff, /dimension)
    yodim = size(yoff, /dimension)
    if max(dim[1:2]-xodim) eq 0 && max(dim[1:2]-yodim) eq 0 then begin
      for i=0, dim[0]-1 do imodmat2[i,*,*] = red_applyoffsets( reform(imodmat[i,*,*]), xoff, yoff )
      do_copy = 0
    endif
  endif
  if do_copy gt 0 then begin
    for i=0, dim[0]-1 do imodmat2[i,*,*] = reform(imodmat[i,*,*])
  endif
  
  ; The momfbd results are transposed so we need to also transpose the demodulation matrix
  imodmat2 = transpose( temporary(imodmat2), [0,2,1] )

  mny = dim[2] - 1
  mnx = dim[1] - 1
  dim = size(lc0.patch, /dimension)

  for x = 0, dim[0]-1 do begin
     for y = 0, dim[1]-1 do begin

        pnx = lc0.patch[x,y].xh - lc0.patch[x,y].xl + 1
        pny = lc0.patch[x,y].yh - lc0.patch[x,y].yl + 1
        cpmm = dindgen(4,4,pnx,pny)      
        pmm = dindgen(pnx,pny)

        clip = red_getclips(lc0, x, y, /aligned)
        tpmm = imodmat2[*,(clip[0]>0)<mnx:(clip[1]<mnx),(clip[2]>0):(clip[3]<mny)]
        for l=0,3 do begin
           pmm[*,*] = median(tpmm[l,*,*])                       ; fill extra bits with the mean...
           pmm[(-clip[0])>0,(-clip[2])>0]=reform(tpmm[l,*,*])   ; copy what we do have
           imm[l,0].patch[x,y].img = red_convolve(red_fillnan(pmm), reform(lc0.patch[x,y].psf[*,*,0]))
        end

        clip = red_getclips(lc1, x, y, /aligned )
        tpmm = imodmat2[*,(clip[0]>0)<mnx:(clip[1]<mnx),(clip[2]>0):(clip[3]<mny)]
        for l=0,3 do begin
           pmm[*,*]=median(tpmm[4+l,*,*])                       ; fill extra bits with the mean...
           pmm[(-clip[0])>0,(-clip[2])>0]=reform(tpmm[4+l,*,*]) ; copy what we do have
           imm[l,1].patch[x,y].img = red_convolve(red_fillnan(pmm), reform(lc1.patch[x,y].psf[*,*,0]))
        end

        clip = red_getclips(lc2, x, y, /aligned)
        tpmm = imodmat2[*,(clip[0]>0)<mnx:(clip[1]<mnx),(clip[2]>0):(clip[3]<mny)]
        for l=0,3 do begin
           pmm[*,*]=median(tpmm[8+l,*,*])                       ; fill extra bits with the mean...
           pmm[(-clip[0])>0,(-clip[2])>0]=reform(tpmm[8+l,*,*]) ; copy what we do have
           imm[l,2].patch[x,y].img = red_convolve(red_fillnan(pmm), reform(lc2.patch[x,y].psf[*,*,0]))
        end

        clip = red_getclips(lc3, x, y, /aligned)
        tpmm = imodmat2[*,(clip[0]>0)<mnx:(clip[1]<mnx),(clip[2]>0):(clip[3]<mny)]
        for l=0,3 do begin
           pmm[*,*]=median(tpmm[12+l,*,*])                       ; fill extra bits with the mean...
           pmm[(-clip[0])>0,(-clip[2])>0]=reform(tpmm[12+l,*,*]) ; copy what we do have
           imm[l,3].patch[x,y].img = red_convolve(red_fillnan(pmm), reform(lc3.patch[x,y].psf[*,*,0]))
        end
     end
  end

  return, imm
  
end
