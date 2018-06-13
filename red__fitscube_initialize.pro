; docformat = 'rst'

;+
; Create a fits file to write data cube slices into when you don't
; want to have the entire cube in memory at the same time.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl
; 
; 
; 
; :Params:
; 
;    filename : in, type=string
;
;       The name of the file to create.
;
;
;    hdr : in, type=strarr
;
;       The complete fits header of the file. The NAXIS* and BITPIX
;       keywords are used to set the file up.
;
;
;    lun : out, type=integer
;
;       The logical unit of the opened file.
;
;
;    fileassoc : out, type="file assoc variable"
;
;       An assoc variable that can be used to write (and read!) data
;       slices. 
;
;    dimensions : in, type=intarr
;
;       The dimensions of the data cube.
;
; :History:
; 
;    2016-03-24 : MGL. First version.
; 
;    2017-06-05 : MGL. Remove WCS keywords.
;
;    2017-09-07 : MGL. Changed red_fitsaddpar --> red_fitsaddkeyword. 
;
;    2017-10-06 : MGL. Headers for tabulating spatial coordinates.
;
;    2017-10-18 : MGL. Move adding of WCS related file header keywords
;                 to fitscube_addwcs method.
; 
;-
pro red::fitscube_initialize, filename, hdr, lun, fileassoc, dimensions
  
  ;; This kind of FITS cube is always five-dimensional, but dimensions
  ;; can be degenerate.
  if n_elements(dimensions) ge 1 then Nx      = long(dimensions[0]) else Nx = 1L
  if n_elements(dimensions) ge 2 then Ny      = long(dimensions[1]) else Ny = 1L
  if n_elements(dimensions) ge 3 then Ntuning = long(dimensions[2]) else Ntuning = 1L
  if n_elements(dimensions) ge 4 then Nstokes = long(dimensions[3]) else Nstokes = 1L
  if n_elements(dimensions) ge 5 then Nscans  = long(dimensions[4]) else Nscans = 1L

  ;; Make sure the headers describing dimensions and coordinates are
  ;; correct.
  dimensions = [Nx, Ny, Ntuning, Nstokes, Nscans]
  fxaddpar, hdr, 'NAXIS', 5, 'Number of data axes'
  for i = 0, 4 do fxaddpar, hdr, 'NAXIS'+strtrim(i+1, 2), dimensions[i]

  ;; Some more standard edits of the header
  red_fitsaddkeyword, anchor = anchor, AFTER = 'SIMPLE' $
                      , hdr, 'DATE', red_timestamp(/iso) $                      ; DATE with time
                      , 'Creation UTC date of FITS header'                      ;
  red_fitsaddkeyword, anchor = anchor, hdr, 'FILENAME', file_basename(filename) ; New file name

  ;; VAR_KEYS keyword will have to be added again if variable keywords
  ;; are added or copied.
  red_fitsdelkeyword, hdr, 'VAR_KEYS' 

  ;; Open the new fits file
  openw, lun, filename, /get_lun, /swap_if_little_endian

  ;; Strip trailing empty lines from hdr
  hdr = hdr[0:where(strmid(hdr,0,3) eq 'END')]
  
  ;; Make byte-version of header and write it to the file.
  Nlines = n_elements(hdr)     
  bhdr=reform(byte(hdr),80L*Nlines)
  ;bhdr = replicate(32B, 80L*Nlines)
  ;for n = 0L, Nlines-1 do bhdr[80*n] = byte(hdr[n])
  writeu, lun, bhdr

  ;; Pad to mod 2880 bytes
  Npad = 2880 - (80L*Nlines mod 2880) ; Number of bytes of padding
  if Npad eq 2880 then Npad = 0
  if Npad gt 0 then writeu, lun, replicate(32B,Npad)
  StartPos = 80L*Nlines + Npad
  print, 'StartPos=', StartPos

  ;; Now prepare for the data part
  Nslices = Ntuning * Nstokes * Nscans
  Nelements = Nx * Ny * Nslices

  bitpix = fxpar(hdr, 'BITPIX')
  case bitpix of
     16 : array_structure = intarr(Nx, Ny)
    -32 : array_structure = fltarr(Nx, Ny)
    else : stop
  endcase
;  help, array_structure
;  print, 'Nelements', Nelements

  ;; Pad the file, including the data part, to mod 2880 bytes
  EndPos = StartPos + Nelements*abs(bitpix)/8 ; End position of data
  Npad = 2880 - (EndPos mod 2880)
  if Npad gt 0 and Npad lt 2880 then begin
    bpad = bytarr(Npad)
    rec = assoc(lun,bpad,EndPos)
    rec[0] = bpad
  endif
;  print, 'EndPos=', EndPos
;  print, 'Npad', Npad

  ;; Set up an assoc variable for writing the frames/slices
  fileassoc = assoc(lun, array_structure, StartPos)

end
