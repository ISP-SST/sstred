; docformat = 'rst'

;+
; Boolean, do we have polarimetric data?
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
;   True if we have T and R cameras, False otherwise.
;
; :History:
; 
;    2025-05-19 : MGL. First version.
; 
;-
function red::polarimetric_data

  return, total(strmatch(*self.cameras,'*-[TR]')) gt 0

end
