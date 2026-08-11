; docformat = 'rst'

;+
; Strings that match raw data directories of CRISP2 as well as CRISP
; after the 2022 camera upgrade.
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
;    for CRISP2 as well as CRISP after the 2022 camera upgrade.
; 
; :History:
; 
;     2022-07-29 : MGL. First version.
; 
;-
function crisp2::raw_data_dirs, dummy

  ;; Should match directories with 'CRISP' tags as well as with
  ;; 'CRISP2' tags. Assume both kinds of data do not exist for the
  ;; same observing day.
  
  return, '*CRISP*-' + ['darks', 'flats', 'pinholes', 'polcal', 'data'] + '*'

end


