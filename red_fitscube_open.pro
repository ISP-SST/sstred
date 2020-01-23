; docformat = 'rst'

;+
; Open a fitscube file.
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
; :Params:
; 
;   filename : in, type=string
; 
;      The name of the fitscube file to open.
; 
;   fileassoc : out, optional, type="assoc variable"
;   
;      An assoc variable set up to access the frames in the file.
;
;   fitscube_info : out, optional, type=struct
;
;      Info about the opened fitscube file, including the header, the
;      cube dimensions, and the logical unit.
;
; :Keywords:
; 
;   header : out, optional, type=strarr
; 
;      The header of the file.
;
;   lun : out, optional, type=integer
;
;      The logical unit of the opened file.
;
;   update : in, optional, type=boolean
;
;      Open with openu rather than openr.
; 
; :History:
; 
;      2019-04-03 : MGL. First version.
; 
;      2019-08-23 : MGL. New keyword update.
;
;-
pro red_fitscube_open, filename, fileassoc, fitscube_info $
                       , header = header $
                       , lun = lun $
                       , update = update

  ;; Open the file
  if keyword_set(update) then begin
    openu, lun, filename, /get_lun, /swap_if_little_endian
  endif else begin
    openr, lun, filename, /get_lun, /swap_if_little_endian
  endelse
  
  ;; Was opening all we wanted?
  if ~arg_present(header) $
     and ~arg_present(fileassoc) $
     and ~arg_present(fitscube_info) then return
  
  ;; Header
  header = headfits(lun)

  if ~arg_present(fileassoc) $
     and ~arg_present(fitscube_info) then return

  ;; Header info
  naxis = fxpar(header, 'NAXIS*')
  Nx      = naxis[0]
  Ny      = naxis[1]
  Ntuning = naxis[2]
  Nstokes = naxis[3]
  Nscans  = naxis[4]
  bitpix = fxpar(header, 'BITPIX')

  ;; Find starting point of data
  Nlines = where(strmatch(header, 'END *'), Nmatch) + 1 ; Starts at 0, so add one
;  Npad = 2880 - (80L*Nlines mod 2880)
  Nblock = (Nlines-1)*80/2880+1 ; Number of 2880-byte blocks
  offset = Nblock*2880          ; Offset to start of data

  ;; Set up an assoc variable.
  case bitpix of
    16  : array_structure = intarr(Nx, Ny)
    -32 : array_structure = fltarr(Nx, Ny)
    else : stop
  endcase

  fileassoc = assoc(lun, array_structure, offset)
  ;; Some information is available in the struct returned after
  ;; calling x=size(fileassoc,/struct). Among them the assoc variable
  ;; type in x.type and x.type_name, as well as its dimensions in
  ;; x.n_elements, x.n_dimensions, and x.dimensions. Also the file
  ;; offset (space reserved for the header) is available in
  ;; x.file_offset.

  if ~arg_present(fitscube_info) then return

  ;; Set up a struct with info about the fitscube.
  fitscube_info = {dimensions:  [Nx, Ny, Ntuning, Nstokes, Nscans] $
                   , lun:       lun $
                   , header:    header $
                  }
  ;; Some information is available in the struct returned after
  ;; calling x=fstat(lun). Among them the filename x.name and the file
  ;; size in bytes x.size.
  

  
end
