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
pro red::fitscube_initialize, filename, hdr, lun, fileassoc, dimensions $
                              , wcs_time_coordinate = wcs_time_coordinate $
                              , wcs_wave_coordinate = wcs_wave_coordinate $
                              , scannumber = scannumber
  
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

  ;; If the main header didn't have the EXTEND keyword, add it now.
  red_fitsaddpar, hdr, 'EXTEND', !true, 'The file has extension(s).'

  ;; The CTYPE{1,2} keywords should change to HPLN-TAN/HPLT-TAN if we
  ;; know the position and orientation. The position coordinates then
  ;; go into CRVAL{1,2}. /MGL


  ;; No rotations
  red_fitsaddpar, hdr, 'PC1_1', 1.0, 'No rotations', anchor = anchor, after = 'FILENAME', /force
  red_fitsaddpar, hdr, 'PC2_2', 1.0, 'No rotations', anchor = anchor
  red_fitsaddpar, hdr, 'PC3_3', 1.0, 'No rotations', anchor = anchor
  red_fitsaddpar, hdr, 'PC4_4', 1.0, 'No rotations', anchor = anchor
  red_fitsaddpar, hdr, 'PC5_5', 1.0, 'No rotations', anchor = anchor 

  ;; Spatial coordinate 1
  red_fitsaddpar, hdr, 'CTYPE1', 'INSX-TAN', 'Instrument X', anchor = anchor
;  red_fitsaddpar, hdr, 'CTYPE1', 'HPLN-TAN', 'SOLAR X'; , anchor = anchor , after = 'FILENAME'
  red_fitsaddpar, hdr, 'CUNIT1', 'arcsec', 'Unit along axis 1', anchor = anchor 
  red_fitsaddpar, hdr, 'CDELT1', float(self.image_scale), '[arcsec] x-coordinate resolution', anchor = anchor 
  red_fitsaddpar, hdr, 'CRPIX1', (Nx+1.0)/2.0, 'Center pixel in x direction', anchor = anchor 
  red_fitsaddpar, hdr, 'CRVAL1', 0., 'Coord. 1 coordinates start', anchor = anchor 
;  red_fitsaddpar, hdr, 'CRVAL1', ??, '[arcsec] Center x solar coordinate', anchor = anchor 

  ;; Spatial coordinate 2
  red_fitsaddpar, hdr, 'CTYPE2', 'INSY-TAN', 'Instrument Y', anchor = anchor 
;  red_fitsaddpar, hdr, 'CTYPE2', 'HPLT-TAN', '[arcsec]', 'SOLAR Y', anchor = anchor 
  red_fitsaddpar, hdr, 'CUNIT2', 'arcsec', 'Unit along axis 2', anchor = anchor 
  red_fitsaddpar, hdr, 'CDELT2', float(self.image_scale), '[arcsec] y-coordinate resolution', anchor = anchor 
  red_fitsaddpar, hdr, 'CRPIX2', (Ny+1.0)/2.0, 'Center pixel in y direction', anchor = anchor 
  red_fitsaddpar, hdr, 'CRVAL2', 0., 'Coord. 2 coordinates start', anchor = anchor 
;  red_fitsaddpar, hdr, 'CRVAL2', ??, '[arcsec] Center y solar coordinate', anchor = anchor 

  ;; Tuning, tabulated wavelength
  red_fitsaddpar, hdr, 'CTYPE3', 'WAVE-TAB', 'Wavelength, function of tuning and scan number', anchor = anchor 
  red_fitsaddpar, hdr, 'CNAME3', 'Wavelength', anchor = anchor 
  red_fitsaddpar, hdr, 'CUNIT3', 'nm', 'Wavelength unit, tabulated for dim. 3 and 5', anchor = anchor 
  red_fitsaddpar, hdr, 'PS3_0', 'WCS-TAB', 'EXTNAME; EXTVER=EXTLEVEL=1 is default', anchor = anchor 
  red_fitsaddpar, hdr, 'PS3_1', 'WAVE+TIME', 'TTYPE for column w/coordinates', anchor = anchor 
;  red_fitsaddpar, hdr, 'PS3_2', 'WAVE-INDEX', 'TTYPE for INDEX, anchor = anchor 
  red_fitsaddpar, hdr, 'PV3_3', 1, 'Coord. 3 is tabulated coordinate number 1', anchor = anchor 
  red_fitsaddpar, hdr, 'CRPIX3', 0, 'Unity transform', anchor = anchor 
  red_fitsaddpar, hdr, 'CRVAL3', 0, 'Unity transform', anchor = anchor 
  red_fitsaddpar, hdr, 'CDELT3', 1, 'Unity transform', anchor = anchor 
  
  ;; Stokes
  red_fitsaddpar, hdr, 'CTYPE4', 'STOKES', 'Stokes vector [I,Q,U,V]', anchor = anchor
  red_fitsaddpar, hdr, 'CRPIX4', 1, 'First (and only) quantity is I', anchor = anchor 
  red_fitsaddpar, hdr, 'CRVAL4', 1, 'First (and only) quantity is I', anchor = anchor 
  red_fitsaddpar, hdr, 'CDELT4', 1, '[1,2,3,4] = [I,Q,U,V]', anchor = anchor 
    
  ;; Scan number = repetition = major time dimension (although time
  ;; varies during scans as well).
  red_fitsaddpar, hdr, 'CTYPE5', 'TIME-TAB', 'Time, function of tuning and scan number', anchor = anchor 
  red_fitsaddpar, hdr, 'CNAME5', 'Time since DATEREF, increases with dim. 3 and 5', anchor = anchor 
  red_fitsaddpar, hdr, 'CUNIT5', 's', anchor = anchor 
  red_fitsaddpar, hdr, 'PS5_0', 'WCS-TAB',   'EXTNAME; EXTVER=EXTLEVEL=1 is default', anchor = anchor 
  red_fitsaddpar, hdr, 'PS5_1', 'WAVE+TIME', 'TTYPE for column w/coordinates', anchor = anchor 
;  red_fitsaddpar, hdr,'PS5_2', 'TIME-INDEX','TTYPE for INDEX, anchor = anchor 
  red_fitsaddpar, hdr, 'PV5_3', 2, 'Coord. 5 is tabulated coordinate number 2', anchor = anchor 
  red_fitsaddpar, hdr, 'CRPIX5', 0, 'Unity transform', anchor = anchor 
  red_fitsaddpar, hdr, 'CRVAL5', 0, 'Unity transform', anchor = anchor 
  red_fitsaddpar, hdr, 'CDELT5', 1, 'Unity transform', anchor = anchor 
 
  ;; Open the new fits file
  openw, lun, filename, /get_lun, /swap_if_little_endian

  ;; Make byte-version of header and write it to the file.
  Nlines = n_elements(hdr)     
  bhdr = replicate(32B, 80L*Nlines)
  for n = 0L, Nlines-1 do bhdr[80*n] = byte(hdr[n])
  writeu, lun, bhdr

  ;; Pad to mod 2880 bytes
  Npad = 2880 - (80L*Nlines mod 2880) ; Number of bytes of padding
  if npad GT 0 then writeu, lun, replicate(32B,Npad)
  Nblock = (Nlines-1)*80/2880+1 ; Number of 2880-byte blocks
  offset = Nblock*2880          ; Offset to start of data
print, 'offset=', offset
  
  ;; Now prepare for the data part
  Nslices = Ntuning * Nstokes * Nscans
  Nelements = Nx * Ny * Nslices

;  naxis = sxpar(hdr,'NAXIS')
;  dims = lonarr(2)
;  Nelements = 1L
;  Nslices = 1L
;  for i = 1, 2 do begin
;     dims[i-1] = sxpar(hdr,'NAXIS'+strtrim(i, 2))
;     Nelements *= dims[i-1]
;     if i gt 1 then Nslices *= dims[i-1]
;  endfor

  bitpix = fxpar(hdr, 'BITPIX')
  case bitpix of
    16 : array_structure = intarr(Nx, Ny)
    -32 : array_structure = fltarr(Nx, Ny)
    else : stop
  endcase
help, array_structure
print, 'Nelements', Nelements

  ;; Pad the file, including the data part, to mod 2880 bytes
  endpos = offset + Nelements*abs(bitpix)/8 ; End position of data
  Npad = 2880 - (endpos mod 2880)
  if npad GT 0 and npad LT 2880 then begin
    bpad = bytarr(Npad)
    rec = assoc(lun,bpad,endpos)
    rec[0] = bpad
  endif
print, 'endpos=', endpos
print, 'Npad', Npad


  ;; Write the WCS extension. (This was inspired by Stein's
  ;; wcs_crosstabulation.pro)
  if n_elements(wcs_time_coordinate) gt 0 $
     or n_elements(wcs_wave_coordinate) gt 0 $
     or n_elements(scannumbers) gt 0 $
  then begin
    free_lun, lun               ; The file is expected to be closed for this.
    fxbhmake,bdr,1,'EXTNAME-WCS-TABLES','For storing tabulated WCS coordinates'  
    ;; Assume the time and wavelength arrays already have the correct
    ;; dimensions, and that the proper C* keywords are set.
    colno = 1

    if n_elements(wcs_time_coordinate) ne 0 then begin
      fxbaddcol, colno, bdr, wcs_time_coordinate, 'TIME-TABULATION' $
                 , 'Table of times', TUNIT = 's'
      t_added = 1
      colno++
    endif else t_added = 0

    if n_elements(wcs_wave_coordinate) ne 0 then begin
      fxbaddcol, colno, bdr, wcs_wave_coordinate, 'WAVE-TABULATION' $
                 , 'Table of wavelengths', TUNIT = 'm'
      w_added = 1
      colno++
    endif else w_added = 0

    if n_elements(scannumber) ne 0 then begin
      ;; Not part of wcs system, store in other extension? Alternate
      ;; coordinate of the fifth dimension?
      fxbaddcol, colno, bdr, scannumber, 'SCAN-TABULATION' $
                 , 'Table of scan numbers'
      s_added = 1
      colno++
    endif else s_added = 0
print, t_added, w_added,  s_added
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
    if t_added then begin
      fxbwrite, lun, wcs_time_coordinate, colno, 1
      colno++
    endif
    if w_added then begin
      fxbwrite, lun, wcs_wave_coordinate, colno, 1
      colno++
    endif
    if s_added then begin
      fxbwrite, lun, scannumber, colno, 1
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

;  ;; Write other binary extension if any.
;  if n_elements(tabhdu) gt 0 then begin
;    free_lun, lun               ; The file is expected to be closed for this.
;    red_write_tabhdu, tabhdu, filename
;    ;; Open the fits file again, so the image data can be written
;    ;; outside of this subroutine.
;    openu, lun, filename, /get_lun, /swap_if_little_endian
;  endif

  ;; Set up an assoc variable for writing the frames/slices
  fileassoc = assoc(lun, array_structure, offset)

end
