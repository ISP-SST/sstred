; docformat = 'rst'

;+
; Strings that match raw data directories of CHROMIS.
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
;    A string array with wildcard-bracketed raw data directory names
;    for CHROMIS.
; 
; :History:
; 
;     2022-07-29 : MGL. First version.
; 
;-
function chromis::raw_data_dirs, dummy

  return, '*CHROMIS-' + ['darks', 'flats', 'pinholes', 'polcal', 'data'] + '*'

end


