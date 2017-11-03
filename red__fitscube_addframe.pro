; docformat = 'rst'

;+
; Add a frame to a fitscube file.
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
; 
; :Params:
;
;    fileassoc : in, type=integer
;
;       The assoc variable set up to access the main data part of the
;       file in which to add the frame.
;
;    frame : in, type=array
;
;       The image frame to add.
; 
; :Keywords:
; 
;    iframe : in, optional, type=integer, default="based on ituning, istokes, and iscan"
;
;       The frame index in the data cube seen as a 3D cube. 
;
;    ituning  : in, optional, type=integer, default=0
;
;       The tuning index, used to calculate iframe.
;
;    istokes  : in, optional, type=integer, default=0
;
;       The stokes index, used to calculate iframe.
;
;    iscan : in, optional, type=integer, default=0
;
;       The scan index, used to calculate iframe.
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
;    2017-11-03 : MGL. New keyword iframe.
; 
;-
pro red::fitscube_addframe, fileassoc, frame $
                            , iframe = iframe $
                            , ituning = ituning $
                            , istokes = istokes $
                            , iscan = iscan

  if n_elements(iframe) eq 0 then begin
    
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
    
  endif
  
  ;; Write the frame
  fileassoc[iframe] = frame
  
end
