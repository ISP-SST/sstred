; docformat = 'rst'

;+
; Get WAVELNTH keyword from FITS header, taking WAVEUNIT into account.
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
;    The wavelength in meters.
; 
; :Params:
; 
; 
; 
; 
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;    2017-06-19 : MGL. First version.
; 
; 
; 
; 
;-
pro red_fitspar_getwavelnth, hdr $
                             , wavelnth = wavelnth, haswav = haswav $
                             , waveunit = waveunit, hasunit = hasunit

  wavelnth = fxpar(hdr, 'WAVELNTH', count = haswav)
  waveunit = fxpar(hdr, 'WAVEUNIT', count = hasunit)
  
  if haswav and hasunit then wavelnth *= 10d^long(waveunit)
  
end
