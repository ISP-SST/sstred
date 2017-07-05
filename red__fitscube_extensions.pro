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
pro red::fitscube_extensions, filename $
                              , wcs_time_coordinate = wcs_time_coordinate $
                              , wcs_wave_coordinate = wcs_wave_coordinate $
                              , scannumber = scannumber

  ;; Write the WCS extension. (This was inspired by Stein's
  ;; wcs_crosstabulation.pro)

  if max([n_elements(wcs_time_coordinate) $
          , n_elements(wcs_wave_coordinate) $
          , n_elements(scannumbers)]) eq 0 then return
  
  fxbhmake,bdr,1,'WCS-TAB','For storing tabulated WCS coordinates'  
  ;; Assume the time and wavelength arrays already have the correct
  ;; dimensions, and that the proper C* keywords are set.
  colno = 1

  if n_elements(wcs_time_coordinate) gt 0 then begin
    fxbaddcol, colno, bdr, wcs_time_coordinate, 'TIME-TABULATION' $
               , 'Table of times', TUNIT = 's'
    colno++
    print, 'Added time'
  endif
  if n_elements(wcs_wave_coordinate) gt 0 then begin
    fxbaddcol, colno, bdr, wcs_wave_coordinate, 'WAVE-TABULATION' $
               , 'Table of wavelengths', TUNIT = 'm'
    colno++
    print, 'Added wave'
  endif
  if n_elements(scannumber) gt 0 then begin
    ;; Not part of wcs system, store in other extension? Alternate
    ;; coordinate of the fifth dimension?
    fxbaddcol, colno, bdr, scannumber, 'SCAN-TABULATION' $
               , 'Table of scan numbers'
    colno++
    print, 'Added scan'
  endif


;  if n_elements(wcs_time_coordinate) ne 0 then begin
;    fxbaddcol, colno, bdr, wcs_time_coordinate, 'TIME-TABULATION' $
;               , 'Table of times', TUNIT = 's'
;    t_added = 1
;    colno++
;  endif else t_added = 0
;
;  if n_elements(wcs_wave_coordinate) ne 0 then begin
;    fxbaddcol, colno, bdr, wcs_wave_coordinate, 'WAVE-TABULATION' $
;               , 'Table of wavelengths', TUNIT = 'm'
;    w_added = 1
;    colno++
;  endif else w_added = 0
;
;  if n_elements(scannumber) ne 0 then begin
;    ;; Not part of wcs system, store in other extension? Alternate
;    ;; coordinate of the fifth dimension?
;    fxbaddcol, colno, bdr, scannumber, 'SCAN-TABULATION' $
;               , 'Table of scan numbers'
;    s_added = 1
;    colno++
;  endif else s_added = 0
;  print, t_added, w_added,  s_added

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
  if n_elements(wcs_time_coordinate) gt 0 then begin
    fxbwrite, lun, wcs_time_coordinate, colno, 1
    colno++
    print, 'Wrote time'
  endif
  if n_elements(wcs_wave_coordinate) gt 0 then begin
    fxbwrite, lun, wcs_wave_coordinate, colno, 1
    colno++
    print, 'Wrote wave'
  endif
  if n_elements(scannumber) gt 0 then begin
    fxbwrite, lun, scannumber, colno, 1
    colno++
    print, 'Wrote scan'
  endif

  ;;
  ;; Finished - crossing our fingers here!!
  ;;
  fxbfinish, lun

  ;; Write other binary extension if any.
  if n_elements(tabhdu) gt 0 then begin
    free_lun, lun               ; The file is expected to be closed for this.
    red_write_tabhdu, tabhdu, filename
    ;; Open the fits file again, so the image data can be written
    ;; outside of this subroutine.
    openu, lun, filename, /get_lun, /swap_if_little_endian
  endif

end
