; docformat = 'rst'

;+
; Return the grid pitch in arcsec of the pinhole array for a specific
; date.
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
;    The grid pitch in arcsec.
; 
; :Params:
; 
;    date : in, type=string
; 
;      The date in ISO format YYYY-MM-DD.
; 
; :History:
; 
;   2024-06-19 : MGL. First version.
; 
;-
function red_pinhole_pitch_arcsec, isodate

  ;; The pinhole array installed just before the 2023 season.
  ;; Calculated in 2024 from the pitch of the old array and the ratio
  ;; of pitches in pixels of the old and new array.
  if isodate lt 2023-01-01 then return, 4.13 ; ["]

  ;; The pinhole array used 2004-2012. Based on calibration with
  ;; SDO/HMI images in 2013, data from 2012.
  return, 5.12                  ; ["]

end
