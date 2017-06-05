; docformat = 'rst'

;+
; Add extensions to a fits file.
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
;       The name of the file.
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
pro red::fitscube_addwcs, filename, wcs_wave_coordinate, wcs_time_coordinate $
                          , scannumber = scannumber

  ;; Write the WCS extension. (This was inspired by Stein's
  ;; wcs_crosstabulation.pro)

  if max([n_elements(wcs_time_coordinate) $
          , n_elements(wcs_wave_coordinate) $
          , n_elements(scannumbers)]) eq 0 then return

  ;; Get the dimensions. Allow one of them to be scalar, the other has
  ;; to be (Ntunes x Nscans).
  wdims = size(wcs_time_coordinate, /dim)
  tdims = size(wcs_time_coordinate, /dim)
  if n_elements(wdims) eq 1 and n_elements(tdims) eq 1 then begin
    ;; Both are scalars!
    return
  endif else if n_elements(wdims) eq 1 then begin
    Ntunes = tdims[0]
    Nscans = tdims[1]
  endif else begin
    Ntunes = wdims[0]
    Nscans = wdims[1]
  endelse
  
  ;; An array to store both wave and time coordinates
  wave_time_coords = dblarr(2, Ntunes, Nscans)
  wave_time_coords[0, *, *] = wcs_wave_coordinate
  wave_time_coords[1, *, *] = wcs_time_coordinate

  ;; If the main header doesn't have the EXTEND keyword, add it now.
  fxhmodify, filename, 'EXTEND', !true, 'The file has extension(s).'

  ;; Make the binary extension header
  fxbhmake, bdr, 1, 'WCS-TAB', 'For storing tabulated WCS tabulations'
  
  fxbaddcol, 1, bdr, wave_time_coords,  'WAVE+TIME',  'Coordinate array'
  fxbaddcol, 2, bdr, findgen(Ntunes)+1, 'WAVE-INDEX', 'Index for tuning/wavelength'
  fxbaddcol, 3, bdr, findgen(Nscans)+1, 'TIME-INDEX', 'Index for repeat/time'
  
  fxbcreate, unit, filename, bdr, extension_no
  fxbwrite, unit, wave_time_coords, 1, 1 ; Write coordinates as column 1, row 1
  fxbfinish, unit
  

;
;  if n_elements(scannumber) ne 0 then begin
;    ;; Not part of wcs system, store in other extension? Alternate
;    ;; coordinate of the fifth dimension?
;    fxbaddcol, colno, bdr, scannumber, 'SCAN-TABULATION' $
;               , 'Table of scan numbers'
;    colno++
;  endif 

;  ;;
;  ;; Now, for the writing:
;  ;;
;  ;; 1. Create binary table extension w/the header;; bdr, file_unit
;  ;; & extension_no are outputs
;  ;;
;  fxbcreate, lun, filename, bdr, extension_no
;
;  ;;
;  ;; 2. Actually write the data - column 1,  row 1
;  ;;
;  if n_elements(scannumber) gt 0 then begin
;    fxbwrite, lun, scannumber, colno, 1
;    print, 'Wrote scan'
;  endif
;
;  ;;
;  ;; Finished - crossing our fingers here!!
;  ;;
;  fxbfinish, lun

;  ;; Write other binary extension if any.
;  if n_elements(tabhdu) gt 0 then begin
;    free_lun, lun               ; The file is expected to be closed for this.
;    red_write_tabhdu, tabhdu, filename
;    ;; Open the fits file again, so the image data can be written
;    ;; outside of this subroutine.
;    openu, lun, filename, /get_lun, /swap_if_little_endian
;  endif

end
