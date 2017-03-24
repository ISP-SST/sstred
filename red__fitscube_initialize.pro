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
;    head : in, type=strarr
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
; 
; 
;-
pro red::fitscube_initialize, filename, head, lun, fileassoc, dimensions 
  
  ;; This kind of FITS cube is always five-dimensional, but dimensions
  ;; can be degenerate.
  if n_elements(dimensions) ge 1 then Nx = long(dimensions[0])      else Nx = 1L
  if n_elements(dimensions) ge 2 then Ny = long(dimensions[1])      else Ny = 1L
  if n_elements(dimensions) ge 3 then Nlambda = long(dimensions[2]) else Nlambda = 1L
  if n_elements(dimensions) ge 4 then Nstokes = long(dimensions[3]) else Nstokes = 1L
  if n_elements(dimensions) ge 5 then Nscans  = long(dimensions[4]) else Nscans = 1L

  dimensions = [Nx, Ny, Nlambda, Nstokes, Nscans]


  ;; Make sure the headers describing dimensions and coordinates are
  ;; correct. 

  fxaddpar, head, 'NAXIS', 5, 'Number of data axes'
  for i = 0, 4 do fxaddpar, head, 'NAXIS'+strtrim(i, 2), dimensions[i]

  ;; The CTYPE{1,2} keywords should change to HPLN-TAN/HPLT-TAN if we
  ;; know the position and orientation? /MGL

  fxaddpar, head, 'CTYPE1', 'x', '[arcsec]', after='DATE'
  fxaddpar, head, 'CUNIT1', 'arcsec', 'Unit along axix 1'
  fxaddpar, head, 'CDELT1', float(self.image_scale), '[arcsec] x-coordinate increment'

  fxaddpar, head, 'CTYPE2', 'y', '[arcsec]'
  fxaddpar, head, 'CUNIT2', 'arcsec', 'Unit along axix 2'
  fxaddpar, head, 'CDELT2', float(self.image_scale), '[arcsec] y-coordinate increment'

  fxaddpar, head, 'CNAME3', 'Wavelength, increases with dim. 3'
  fxaddpar, head, 'CTYPE3', 'WAVE-TAB',   'Wavelength, tabulated'
  fxaddpar, head, 'CUNIT3', 'm', 'Wavelength unit'
  fxaddpar, head, 'CRPIX3', 1, 'Coord. 3 lower left pixel'
  fxaddpar, head, 'CRVAL3', 1, 'Coord. 3 indexes table starting at 1'
  ;; Should have similar keywords as 5 if tabulated?
  fxaddpar, head, 'PS3_0', 'EXTNAME-WCS-TABLES', 'Extension w/tabulated coordinates'
  fxaddpar, head, 'PS3_1', 'WAVE-TABULATION', 'TTYPE column name for WAVE coords'

  fxaddpar, head, 'CTYPE4', 'Stokes',     '[I,Q,U,V]'

  
  fxaddpar, head, 'CNAME5', 'Time since DATEREF, increases with dim. 3 and 5'
  fxaddpar, head, 'CTYPE5', 'TIME-TAB', 'Time since DATEREF, tabulated'
  fxaddpar, head, 'CUNIT5', 's'
  fxaddpar, head, 'CRPIX5', 1, 'Coord. 5 lower left pixel'
  fxaddpar, head, 'CRVAL5', 1, 'Coord. 5 indexes table starting at 1'
  fxaddpar, head, 'CDELT5', 1, 'Stride is built into PC5_5, not CDELT5'
  fxaddpar, head, 'PS5_0', 'EXTNAME-WCS-TABLES', 'Extension w/tabulated coordinates'
  fxaddpar, head, 'PS5_1', 'TIME-TABULATION', 'TTYPE column name for TIME coords'
  fxaddpar, head, 'PC5_5', Nlambda, 'STRIDE for repetition = '+strtrim(Nlambda, 2)
  fxaddpar, head, 'PC5_3', 1, '5th coord varies with 3rd coord, STRIDE=1'


 
  
  
  ;; Open the new fits file
  openw, lun, filename, /get_lun, /swap_if_little_endian

  hdr = head                    ; Protect input
  if n_elements(tabhdu) gt 0 then begin
    ;; Binary extension header 
    extname = tag_names(tabhdu)
    for i = 0, n_elements(extname)-1 do red_create_extheaders, tabhdu, extname[i], hdr
  endif

  ;; Make byte-version of header and write it to the file.
  Nlines = n_elements(hdr)     ; Number of header lines
  bhdr = replicate(32B, 80L*Nlines)
  for n = 0L, Nlines-1 do bhdr[80*n] = byte(hdr[n])
  writeu, lun, bhdr

  ;; Pad to mod 2880 bytes
  Npad = 2880 - (80L*Nlines mod 2880) ; Number of bytes of padding
  if npad GT 0 then writeu, lun, replicate(32B,Npad)
  Nblock = (Nlines-1)*80/2880+1 ; Number of 2880-byte blocks
  offset = Nblock*2880          ; Offset to start of data

  ;; Now prepare for the data part

  dims = [Nx, Ny]
  Nelements = Nx*Ny
  Nslices = Ny
 
  bitpix = fxpar(head, 'BITPIX')
  case bitpix of
      16 : array_structure = intarr(dims)
     -32 : array_structure = fltarr(dims)
     else : stop
  endcase

  ;; Pad the file, including the data part, to mod 2880 bytes
  endpos = offset + Nelements*abs(bitpix)/8 ; End position of data
  Npad = 2880 - (endpos mod 2880)
  if npad GT 0 and npad LT 2880 then begin
     bpad = bytarr(Npad)
     rec = assoc(lun,bpad,endpos)
     rec[0] = bpad
  endif

  ;; Write the WCS extension. (This was inspired by Stein's
  ;; wcs_crosstabulation.pro)
  if n_elements(wcs_time) gt 0 or n_elements(wcs_wavelength) gt 0 then begin
    free_lun, lun               ; The file is expected to be closed for this.
    fxbhmake,bdr,1,'EXTNAME-WCS-TABLES','For storing tabulated WCS coordinates'  
    ;; Assume the time and wavelength arrays already have the correct
    ;; dimensions, and that the proper C* keywords are set.
    colno = 1
    if n_elements(wcs_time_coordinate) then begin
      fxbaddcol, colno, bdr, wcs_time_coordinate, 'TIME-TABULATION' $
                 , 'Table of times', TUNIT = 's'
      t_added = 1
      colno++
    endif else t_added = 0
    if n_elements(wcs_wavelength_coordinate) then begin
      fxbaddcol, colno, bdr, wcs_wavelength_coordinate, 'WAVE-TABULATION' $
                 , 'Table of wavelengths', TUNIT = 'm'
      w_added = 1
      colno++
    endif else w_added = 0
    ;;
    ;; Now, for the writing:
    ;;
    ;; 1. Create binary table extension w/the header;; bdr, file_unit
    ;; & extension_no are outputs
    ;;
    fxbcreate, lun, filename, bdr, extension_no
    ;;
    ;; 2. Actually write the data - column 1,  row 1
    ;;
    colno = 1
    if n_elements(wcs_time_coordinate) then begin
      fxbwrite, lun, wcs_time_coordinate, colno, 1
      colno++
    endif
    if n_elements(wcs_wavelength_coordinate) then begin
      fxbwrite, lun, wcs_wavelength_coordinate, colno, 1
      colno++
    endif
    ;;
    ;; Finished - crossing our fingers here!!
    ;;
    fxbfinish, lun
    ;; Open the fits file again, so the image data can be written
    ;; outside of this subroutine.
    openu, lun, filename, /get_lun, /swap_if_little_endian
  endif

  ;; Write other binary extension if any.
  if n_elements(tabhdu) gt 0 then begin
    free_lun, lun               ; The file is expected to be closed for this.
    red_write_tabhdu, tabhdu, filename
    ;; Open the fits file again, so the image data can be written
    ;; outside of this subroutine.
    openu, lun, filename, /get_lun, /swap_if_little_endian
  endif

  ;; Set up an assoc variable for writing the slices
  fileassoc = assoc(lun, array_structure, offset)

end
