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
function red_get_linearcomp, iwav, pp, npar, reflect = reflect
                                ;
  dim = size(pp, /dimension)
  res = fltarr(dim[1], dim[2]) + 1.0
  if(npar le 2 AND ~keyword_set(reflect)) then return, res
  if(npar le 3 AND keyword_set(reflect)) then return, res

                                ;
  wavf = 1.0
                                ;

  ii0 = 2
  if(keyword_set(reflect)) then ii0 += 1 ; Are we also fitting reflectivities?

  for ii = ii0, npar - 1 do begin
     wavf *= iwav
     res += reform(pp[ii,*,*]) * wavf
  endfor
                                ;
  return, res
end
