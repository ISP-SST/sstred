; docformat = 'rst'

;+
; Get a frame from an SST fitscube file.
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
; 
; :Params:
;
;    filename_or_fileassoc : in, type="integer or string"
;
;       If a string, the name of the file from which to get the frame.
;       Otherwise the assoc variable set up to access the main data
;       part of the already opened file.
;
;    frame : out, type=array
;
;       The read image frame.
; 
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
;     2017-11-03 : MGL. First version.
; 
;     2018-01-12 : MGL. Doesn't need to be a method, rename
;                  red::fitscube_getframe to red_fitscube_getframe. 
; 
; 
;-
pro red_fitscube_getframe, filename_or_fileassoc, frame $
                           , iframe = iframe $
                           , ituning = ituning $
                           , istokes = istokes $
                           , iscan = iscan
  
  if size(filename_or_fileassoc, /tname) eq 'STRING' then open_and_close = 1 else open_and_close = 0
  
  if open_and_close then begin
    ;; We have the file name, open the file and set up an assoc
    ;; variable.
    filename = filename_or_fileassoc
    red_fitscube_open, filename, fileassoc, fitscube_info
    Nx      = fitscube_info.dimensions[0]
    Ny      = fitscube_info.dimensions[1]
    Ntuning = fitscube_info.dimensions[2]
    Nstokes = fitscube_info.dimensions[3]
    Nscans  = fitscube_info.dimensions[4]
    hdr = fitscube_info.header
  endif else begin
    ;; We have an assoc variable, get array dimensions from the file.
    fileassoc = filename_or_fileassoc
    lun = (size(fileassoc,/struc)).file_lun
    fs = fstat(lun)
    filename = fs.name
    hdr = headfits(filename)
    dimensions = long(fxpar(hdr, 'NAXIS*'))
    Nx      = dimensions[0]
    Ny      = dimensions[1]
    Ntuning = dimensions[2]
    Nstokes = dimensions[3]
    Nscans  = dimensions[4]
  endelse

  if n_elements(iframe) eq 0 then begin

    if n_elements(ituning) eq 0 then ituning = 0L
    if n_elements(istokes) eq 0 then istokes = 0L
    if n_elements(iscan)   eq 0 then iscan   = 0L
  
    ;; Calculate the frame number
    iframe = long(ituning) + long(istokes)*Ntuning $
             + long(iscan)*Ntuning*Nstokes
  endif

  if (size(fileassoc,/struc)).type_name eq 'FLOAT' then begin
    ;; Get the frame
    frame = fileassoc[iframe]
  endif else begin
    ;; Transform integer data to float
    bzero  = fxpar(hdr, 'BZERO',  count=nzero )
    bscale = fxpar(hdr, 'BSCALE', count=nscale)
    if nzero eq 1 and nscale eq 1 then begin
      ;; Get the frame and rescale
      frame = float(fileassoc[iframe]*bscale + bzero)
    endif else begin
      ;; This should not happen but just get the integer frame if it
      ;; does.
      frame = fileassoc[iframe]
    endelse
  endelse
  
  ;; Close if we opened.
  if open_and_close then red_fitscube_close, fileassoc, fitscube_info

end
