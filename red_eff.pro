; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    SST pipeline
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
;   dmm : in
; 
; 
; 
; :History:
; 
;   2018-04-11 : Split from red__polcal.pro.
; 
;-
function red_eff, dmm

  return, 1./sqrt(4*total(dmm^2, 1))

end
