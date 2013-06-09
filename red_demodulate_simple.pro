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
;   immt1 : 
;   
;   
;   
;    l0 : 
;   
;   
;   
;    l1 : 
;   
;   
;   
;    l2 : 
;   
;   
;   
;    l3 : 
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
function red_demodulate_simple,immt1, l0, l1, l2, l3
  inam = 'red_demodulate_simple'
  dim = size(l0 ,/dimension)
  res = fltarr(dim[0], dim[1], 4)
                                ;
  immt = reform(immt1, [4, 4, dim[0], dim[1]])
                                ;
  for z = 0, 3 do begin
     res[*,*,z]=reform(immt[0,z,*,*]) * l0 +$
                reform(immt[1,z,*,*]) * l1 +$
                reform(immt[2,z,*,*]) * l2 +$ 
                reform(immt[3,z,*,*]) * l3
  endfor
                                ;
  return,res
end
