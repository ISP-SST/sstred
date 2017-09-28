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
;    2017-06-05 : MGL. Remove WCS keywords.
;
;    2017-09-07 : MGL. Changed red_fitsaddpar --> red_fitsaddkeyword. 
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

  ;; If the main header didn't have the EXTEND keyword, add it now.
  red_fitsaddkeyword, hdr, 'EXTEND', !true, 'The file has extension(s).'

  ;; The CTYPE{1,2} keywords should change to HPLN-TAN/HPLT-TAN if we
  ;; know the position and orientation. The position coordinates then
  ;; go into CRVAL{1,2}. /MGL


  ;; No rotations
  red_fitsaddkeyword, hdr, 'PC1_1', 1.0, 'No rotations', anchor = anchor, after = 'FILENAME', /force
  red_fitsaddkeyword, hdr, 'PC2_2', 1.0, 'No rotations', anchor = anchor
  red_fitsaddkeyword, hdr, 'PC3_3', 1.0, 'No rotations', anchor = anchor
  red_fitsaddkeyword, hdr, 'PC4_4', 1.0, 'No rotations', anchor = anchor
  red_fitsaddkeyword, hdr, 'PC5_5', 1.0, 'No rotations', anchor = anchor 

  ;; Spatial coordinate 1
  red_fitsaddkeyword, hdr, 'CTYPE1', 'INSX-TAN', 'Instrument X', anchor = anchor
;  red_fitsaddkeyword, hdr, 'CTYPE1', 'HPLN-TAN', 'SOLAR X'; , anchor = anchor , after = 'FILENAME'
  red_fitsaddkeyword, hdr, 'CUNIT1', 'arcsec', 'Unit along axis 1', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CDELT1', float(self.image_scale), '[arcsec] x-coordinate resolution', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CRPIX1', (Nx+1.0)/2.0, 'Center pixel in x direction', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CRVAL1', 0., 'Coord. 1 coordinates start', anchor = anchor 
;  red_fitsaddkeyword, hdr, 'CRVAL1', ??, '[arcsec] Center x solar coordinate', anchor = anchor 

  ;; Spatial coordinate 2
  red_fitsaddkeyword, hdr, 'CTYPE2', 'INSY-TAN', 'Instrument Y', anchor = anchor 
;  red_fitsaddkeyword, hdr, 'CTYPE2', 'HPLT-TAN', '[arcsec]', 'SOLAR Y', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CUNIT2', 'arcsec', 'Unit along axis 2', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CDELT2', float(self.image_scale), '[arcsec] y-coordinate resolution', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CRPIX2', (Ny+1.0)/2.0, 'Center pixel in y direction', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CRVAL2', 0., 'Coord. 2 coordinates start', anchor = anchor 
;  red_fitsaddkeyword, hdr, 'CRVAL2', ??, '[arcsec] Center y solar coordinate', anchor = anchor 

  ;; Tuning, tabulated wavelength
  red_fitsaddkeyword, hdr, 'CTYPE3', 'WAVE-TAB', 'Wavelength, function of tuning and scan number', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CNAME3', 'Wavelength', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CUNIT3', 'nm', 'Wavelength unit, tabulated for dim. 3 and 5', anchor = anchor 
  red_fitsaddkeyword, hdr, 'PS3_0', 'WCS-TAB', 'EXTNAME; EXTVER=EXTLEVEL=1 is default', anchor = anchor 
  red_fitsaddkeyword, hdr, 'PS3_1', 'WAVE+TIME', 'TTYPE for column w/coordinates', anchor = anchor 
;  red_fitsaddkeyword, hdr, 'PS3_2', 'WAVE-INDEX', 'TTYPE for INDEX, anchor = anchor 
  red_fitsaddkeyword, hdr, 'PV3_3', 1, 'Coord. 3 is tabulated coordinate number 1', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CRPIX3', 0, 'Unity transform', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CRVAL3', 0, 'Unity transform', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CDELT3', 1, 'Unity transform', anchor = anchor 
  
  ;; Stokes
  red_fitsaddkeyword, hdr, 'CTYPE4', 'STOKES', 'Stokes vector [I,Q,U,V]', anchor = anchor
  red_fitsaddkeyword, hdr, 'CRPIX4', 1, 'First (and only) quantity is I', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CRVAL4', 1, 'First (and only) quantity is I', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CDELT4', 1, '[1,2,3,4] = [I,Q,U,V]', anchor = anchor 
    
  ;; Scan number = repetition = major time dimension (although time
  ;; varies during scans as well).
  red_fitsaddkeyword, hdr, 'CTYPE5', 'UTC--TAB', 'Time, function of tuning and scan number', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CNAME5', 'Time since DATEREF, increases with dim. 3 and 5', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CUNIT5', 's', anchor = anchor 
  red_fitsaddkeyword, hdr, 'PS5_0', 'WCS-TAB',   'EXTNAME; EXTVER=EXTLEVEL=1 is default', anchor = anchor 
  red_fitsaddkeyword, hdr, 'PS5_1', 'WAVE+TIME', 'TTYPE for column w/coordinates', anchor = anchor 
;  red_fitsaddkeyword, hdr,'PS5_2', 'TIME-INDEX','TTYPE for INDEX, anchor = anchor 
  red_fitsaddkeyword, hdr, 'PV5_3', 2, 'Coord. 5 is tabulated coordinate number 2', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CRPIX5', 0, 'Unity transform', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CRVAL5', 0, 'Unity transform', anchor = anchor 
  red_fitsaddkeyword, hdr, 'CDELT5', 1, 'Unity transform', anchor = anchor 
 
  ;; Open the new fits file
  openw, lun, filename, /get_lun, /swap_if_little_endian

  ;; Make byte-version of header and write it to the file.
  Nlines = n_elements(hdr)     
  bhdr = replicate(32B, 80L*Nlines)
  for n = 0L, Nlines-1 do bhdr[80*n] = byte(hdr[n])
  writeu, lun, bhdr

  ;; Pad to mod 2880 bytes
  Npad = 2880 - (80L*Nlines mod 2880) ; Number of bytes of padding
  if Npad GT 0 then writeu, lun, replicate(32B,Npad)
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

  ;; Set up an assoc variable for writing the frames/slices
  fileassoc = assoc(lun, array_structure, offset)

end
