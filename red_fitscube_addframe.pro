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
;    fitscube_info : in/out, optional, type=structure
;
;       Structure with dimensions, lun and header.
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
;    2020-03-23 : OA. Added 'fitscube_info' keyword.
; 
;-
pro red_fitscube_addframe, filename_or_fileassoc, frame $
                           , iframe = iframe $
                           , ituning = ituning $
                           , istokes = istokes $
                           , iscan = iscan $
                           , fitscube_info = fitscube_info
  
  if size(filename_or_fileassoc, /tname) eq 'STRING' then open_and_close = 1 else open_and_close = 0

  if open_and_close then begin
    ;; We have the file name, open the file and set up an assoc
    ;; variable.
    filename = filename_or_fileassoc
    red_fitscube_open, filename, fileassoc, fitscube_info, /update
    Nx      = fitscube_info.dimensions[0]
    Ny      = fitscube_info.dimensions[1]
    Ntuning = fitscube_info.dimensions[2]
    Nstokes = fitscube_info.dimensions[3]
    Nscans  = fitscube_info.dimensions[4]
  endif else begin
    ;; We have an assoc variable, get array dimensions from the file.
    fileassoc = filename_or_fileassoc
    lun = (size(fileassoc,/struc)).file_lun
    if ~keyword_set(fitscube_info) then begin
      fs = fstat(lun)
      filename = fs.name
      hdr = headfits(filename)
      dimensions = long(fxpar(hdr, 'NAXIS*'))
      Nx      = dimensions[0]
      Ny      = dimensions[1]
      Ntuning = dimensions[2]
      Nstokes = dimensions[3]
      Nscans  = dimensions[4]
      fitscube_info = {dimensions:  [Nx, Ny, Ntuning, Nstokes, Nscans] $
                   , lun:       lun $
                   , header:    hdr $
                  }
    endif else begin
      Nx = fitscube_info.dimensions[0]
      Ny = fitscube_info.dimensions[1]
      Ntuning = fitscube_info.dimensions[2]
      Nstokes = fitscube_info.dimensions[3]
      Nscans  = fitscube_info.dimensions[4]
    endelse
  endelse

  sz = size(frame,/dimensions)
  if sz[0] ne Nx or sz[1] ne Ny then begin  ; shouldn't happen
    print,"Size of frame doesn't correspond to the size of the cube. Stop in red_fitscube_addframe."
    stop
  endif
  
    
  if n_elements(iframe) eq 0 then begin
    
    if n_elements(ituning) eq 0 then ituning = 0L
    if n_elements(istokes) eq 0 then istokes = 0L
    if n_elements(iscan)   eq 0 then iscan   = 0L   
    
    ;; Calculate the frame number
    iframe = long(ituning) + long(istokes)*Ntuning $
             + long(iscan)*Ntuning*Nstokes
    
  endif
  
  ;; Write the frame
  fileassoc[iframe] = frame
  
  ;; Close if we opened.
  if open_and_close then red_fitscube_close, fileassoc, fitscube_info
  
end
