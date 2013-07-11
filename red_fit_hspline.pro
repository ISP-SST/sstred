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
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
;   
;   2013-07-11 : MGL. Use red_intepf, not intepf.
; 
; 
;-
function red_fit_hspline, pp, iwav = iwav, fl = fl, wl = wl
  res = float(red_intepf(iwav, pp, wl, /linear)) ; linear extrapolation outside bounds
  return, (res - fl)
end
