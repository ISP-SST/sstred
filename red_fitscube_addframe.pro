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
    hdr = headfits(filename)
    bitpix = fxpar(hdr, 'BITPIX')
    Nx = fxpar(hdr, 'NAXIS1')
    Ny = fxpar(hdr, 'NAXIS2')
    case bitpix of
      16 : array_structure = intarr(Nx, Ny)
      -32 : array_structure = fltarr(Nx, Ny)
      else : stop
    endcase
    Nlines = where(strmatch(hdr, 'END *'), Nmatch)
    Npad = 2880 - (80L*Nlines mod 2880)
    Nblock = (Nlines-1)*80/2880+1 ; Number of 2880-byte blocks
    offset = Nblock*2880          ; Offset to start of data
    ;; Must be an existing file!
    openu, lun, filename, /get_lun, /swap_if_little_endian
    fileassoc = assoc(lun, array_structure, offset)
  endif else begin
    ;; We have an assoc variable, get array dimensions from the file.
    fileassoc = filename_or_fileassoc
    lun = (size(fileassoc,/struc)).file_lun
    if ~keyword_set(fitscube_info) then begin
      fs = fstat(lun)
      filename = fs.name
      hdr = headfits(filename)
      Nx = fxpar(hdr, 'NAXIS1')
      Ny = fxpar(hdr, 'NAXIS2')
    endif else begin
      Nx = fitscube_info.dimensions[0]
      Ny = fitscube_info.dimensions[1]
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
    
    ;; Get dimensions from the file
    if ~keyword_set(fitscube_info) then begin
      Ntuning = fxpar(hdr,'NAXIS3')
      Nstokes = fxpar(hdr,'NAXIS4')
      Nscans = fxpar(hdr,'NAXIS5')
    endif else begin
      if open_and_close then begin ;shouldn't happen
        print, "One should use 'fitscube_info' only with fileassoc. Stop in red_fitscube_addframe."
        stop
      endif
      Ntuning = fitscube_info.dimensions[2]
      Nstokes = fitscube_info.dimensions[3]
      Nscans  = fitscube_info.dimensions[4]
    endelse
    
    ;; Calculate the frame number
    iframe = long(ituning) + long(istokes)*Ntuning $
             + long(iscan)*Ntuning*Nstokes
    
  endif
  
  ;; Write the frame
  fileassoc[iframe] = frame
  
  ;; Close if we opened.
  if open_and_close then free_lun, lun
  
end
