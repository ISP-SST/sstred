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
;-
function red_matrix2momfbd, lc0, lc1, lc2, lc3, imodmat

  inam = 'red_matrix2momfbd : '
  chan = 0
  nx = abs(lc0.clip[chan,0,1] - lc0.clip[chan,0,0]) + 1
  ny = abs(lc0.clip[chan,1,1] - lc0.clip[chan,1,0]) + 1
  print,' '
  print,'Clipping modulation matrix to '+red_stri(nx)+' x '+red_stri(ny)

  
  imm = replicate({patch:lc0.patch, clip:lc0.clip},[4,4])

  align_clip = [ lc0.clip[chan,0,0], lc0.clip[chan,0,1], $
                 lc0.clip[chan,1,0], lc0.clip[chan,1,1] ]

  dim = size(imodmat, /dimension)
  imodmat2 = make_array( dim[0], nx, ny, type=size(imodmat, /type) )
  for i=0, dim[0]-1 do begin
      imodmat2[i,*,*] = red_clipim( reform(imodmat[i,*,*],dim[1],dim[2]), align_clip )
  endfor
  
  dim=size(imodmat2, /dimension)
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
           dip = size(pmm, /dim)
           res = fltarr(dip)
           imm[l,0].patch[x,y].img = red_convolve(red_fillnan(pmm), reform(lc0.patch[x,y].psf[*,*,0]))
        end

        clip = red_getclips(lc1, x, y, /aligned )
        tpmm = imodmat2[*,(clip[0]>0)<mnx:(clip[1]<mnx),(clip[2]>0):(clip[3]<mny)]
        for l=0,3 do begin
           pmm[*,*]=median(tpmm[4+l,*,*])                       ; fill extra bits with the mean...
           pmm[(-clip[0])>0,(-clip[2])>0]=reform(tpmm[4+l,*,*]) ; copy what we do have

           dip = size(pmm, /dim)
           res = fltarr(dip)
           imm[l,1].patch[x,y].img = red_convolve(red_fillnan(pmm), reform(lc1.patch[x,y].psf[*,*,0]))
        end

        clip = red_getclips(lc2, x, y, /aligned)
        tpmm = imodmat2[*,(clip[0]>0<mnx):(clip[1]<mnx),(clip[2]>0):(clip[3]<mny)]
        for l=0,3 do begin
           pmm[*,*]=median(tpmm[8+l,*,*])                       ; fill extra bits with the mean...
           pmm[(-clip[0])>0,(-clip[2])>0]=reform(tpmm[8+l,*,*]) ; copy what we do have

           dip = size(pmm, /dim)
           res = fltarr(dip)
           imm[l,2].patch[x,y].img = red_convolve(red_fillnan(pmm), reform(lc2.patch[x,y].psf[*,*,0]))
        end

        clip = red_getclips(lc3, x, y, /aligned)
        tpmm = imodmat2[*,(clip[0]>0<mnx):(clip[1]<mnx),(clip[2]>0):(clip[3]<mny)]
        for l=0,3 do begin
           pmm[*,*]=median(tpmm[12+l,*,*])                       ; fill extra bits with the mean...
           pmm[(-clip[0])>0,(-clip[2])>0]=reform(tpmm[12+l,*,*]) ; copy what we do have

           dip = size(pmm, /dim)
           res = fltarr(dip)
           imm[l,3].patch[x,y].img = red_convolve(red_fillnan(pmm), reform(lc3.patch[x,y].psf[*,*,0]))
        end
     end
  end

  return, imm
  
end
