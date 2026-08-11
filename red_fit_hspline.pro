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
;    pp : 
;   
;   
;   
; 
; :Keywords:
; 
;    iwav  : 
;   
;   
;   
;    fl  : 
;   
;   
;   
;    wl  : 
;   
;   
;   
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
;   
;   2013-07-11 : MGL. Use red_intepf, not intepf.
; 
; 
;-
function red_fit_hspline, pp, iwav = iwav, fl = fl, wl = wl, bezier = bezier

  if keyword_set(bezier) then begin
    res = float(red_bezier3(iwav, pp, wl, /linear)) ; linear extrapolation outside bound
  endif else begin
    res = float(red_intepf(iwav, pp, wl, /linear)) ; linear extrapolation outside bounds
  endelse

  res = temporary(res) - fl
  return, temporary(res)

end
