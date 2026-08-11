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
;    hh : 
;   
;   
;   
;    mi : 
;   
;   
;   
;    ss : 
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
function red_reform_frac_time, hh, mi, ss
  return, hh * 3600.d + mi * 60.d0 + ss
end
