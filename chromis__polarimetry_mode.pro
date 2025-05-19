; docformat = 'rst'

;+
; Boolean, are we running Chromis in polarimetry mode?
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
; 
;   True if we have T and R cameras.
;
; :History:
; 
;    2025-05-19 : MGL. First version.
; 
;-
function chromis::polarimetry_mode

  return, total(strmatch(*self.cameras,'*-[TR]')) gt 0

end
