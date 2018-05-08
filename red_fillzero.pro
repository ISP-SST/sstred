; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
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
;    var : 
;   
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
function red_fillzero, var

  idx = where(var lt 1.e-5, count, complement = idx1)
  if(count gt 1) then var[idx] = median(var[idx1])

  return, var

end
