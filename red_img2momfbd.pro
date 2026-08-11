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
;    img : 
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
; 
;-
function red_img2momfbd, lc0, img
  inam = 'red_img2momfbd : '
  chan = 0
  nx = lc0.clip[chan,0,1] - lc0.clip[chan,0,0] + 1
  ny = lc0.clip[chan,1,1] - lc0.clip[chan,1,0] + 1
  print,' '
  print,inam +'Clipping modulation matrix to '+red_stri(nx)+' x '+red_stri(ny)

                                ;
  dim=size(img, /dimension)
                                ;
                                ;
  mny = dim[1] - 1
  mnx = dim[0] - 1
                                ;
  imm = {patch:lc0.patch, clip:lc0.clip}
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
        tpmm=img[(clip[0]>0)<mnx:(clip[1]<mnx),(clip[2]>0):(clip[3]<mny)]
        for l=0,0 do begin
           pmm[*,*]=median(tpmm[*,*])                       ; fill extra bits with the mean...
           pmm[(-clip[0])>0,(-clip[2])>0]=reform(tpmm[*,*]) ; copy what we do have
           
           dip = size(pmm, /dim)
           res = fltarr(dip)
           imm.patch[x,y].img = red_convolve(red_fillnan(pmm), reform(lc0.patch[x,y].psf[*,*,0]))
           
        end
     end
  end
  return, imm
end
