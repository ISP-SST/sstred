; docformat = 'rst'

;+
; Strings that match raw data directories of CRISP with old Sarnoff
; cameras.
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
;    for CRISP with the old Sarnoff cameras.
; 
; :History:
; 
;     2022-07-29 : MGL. First version.
; 
;-
function crisp::raw_data_dirs, dummy

  ;; CRISP with old cameras. We presume that we will not populate the
  ;; database with *really* old CRISP data with different directory
  ;; names
  return, '*' + ['Darks','Flats','Pinholes','Polcal','Science'] + '*'

end

