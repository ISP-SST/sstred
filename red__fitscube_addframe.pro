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
;    2016-03-24 : MGL. First version.
; 
;    2017-11-01 : MGL. Get the cube dimensions from the file instead
;                 of from keywords.
; 
;-
pro red::fitscube_addframe, fileassoc, frame $
                              , ituning = ituning $
                              , istokes = istokes $
                              , iscan = iscan

  if n_elements(ituning) eq 0 then ituning = 0L
  if n_elements(istokes) eq 0 then istokes = 0L
  if n_elements(iscan)   eq 0 then iscan   = 0L
  
  ;; Get dimensions from the file
  lun = (size(fileassoc,/struc)).file_lun
  fs = fstat(lun)
  dimensions = long(fxpar(headfits(fs.name), 'NAXIS*'))
  Nx      = dimensions[0]
  Ny      = dimensions[1]
  Ntuning = dimensions[2]
  Nstokes = dimensions[3]
  Nscans  = dimensions[4]

  ;; Calculate the frame number
  iframe = long(ituning) + long(istokes)*Ntuning $
           + long(iscan)*Ntuning*Nstokes

  ;; Write the frame
  fileassoc[iframe] = frame
  
end
