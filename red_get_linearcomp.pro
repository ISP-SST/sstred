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
;    iwav : 
;   
;   
;   
;    pp : 
;   
;   
;   
;    npar : 
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
function red_get_linearcomp, iwav, pp, npar
                                ;
  dim = size(pp, /dimension)
  res = fltarr(dim[1], dim[2]) + 1.0
  if(npar le 2) then return, res
                                ;
  wavf = 1.0
                                ;
  for ii = 2L, npar - 1 do begin
     wavf *= iwav
     res += pp[ii,*,*] * wavf
  endfor
                                ;
  return, res
end
