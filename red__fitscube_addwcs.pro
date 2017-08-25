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
;    2016-08-25 : MGL. Add spatial coordinates.
; 
; 
; 
;-
pro red::fitscube_addwcs, filename $
                          , wcs_hpln_coordinate $
                          , wcs_hplt_coordinate $
                          , wcs_wave_coordinate $
                          , wcs_time_coordinate

  ;; Write the WCS extension. (This was inspired by Stein's
  ;; wcs_crosstabulation.pro)

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Get the dimensions (using x for hpln and y for hplt).
  xdims = size(wcs_hpln_coordinate, /dim)
  ydims = size(wcs_hplt_coordinate, /dim)
  wdims = size(wcs_wave_coordinate, /dim)
  tdims = size(wcs_time_coordinate, /dim)

  Nxdims = n_elements(xdims)
  Nydims = n_elements(ydims)
  Nwdims = n_elements(wdims)
  Ntdims = n_elements(tdims)

  indx0 = where([Nxdims, Nydims, Nwdims, Ntdims] eq 0, N0)
  indx4 = where([Nxdims, Nydims, Nwdims, Ntdims] eq 4, N4)
  
  ;; At least one of the arrays have to be non-scalar. Any non-scalars
  ;; must have the same dimensions: (2 x 2 x Ntunes x Nscans). (Maybe
  ;; we have to allow for degenerate trailing dimensions for scan-only
  ;; cubes?)
  if N0+N4 ne 4 or N4 eq 0 then begin
    print, inam + ' : At least one of the arrays has the wrong dimensions.'
    help, wcs_hpln_coordinate $
          , wcs_hplt_coordinate $
          , wcs_wave_coordinate $
          , wcs_time_coordinate
    stop
    retall
  endif

  ;; Get the dimensions from a non-scalar
  case (where(indx4))[0] of
    0 : dims = xdims
    1 : dims = ydims
    2 : dims = wdims
    3 : dims = tdims
    else: begin
      print, 'This should not happen!'
      stop
    end
  endcase
  Nx     = dims[0]
  Ny     = dims[1]
  Ntunes = dims[2]
  Nscans = dims[3]
  
  ;; An array to store both wave and time coordinates
  wcs_coords = dblarr(4, Nx, Ny, Ntunes, Nscans)
  wcs_coords[0, *, *, *, *] = wcs_hpln_coordinate
  wcs_coords[1, *, *, *, *] = wcs_hplt_coordinate
  wcs_coords[2, *, *, *, *] = wcs_wave_coordinate
  wcs_coords[3, *, *, *, *] = wcs_time_coordinate

  ;; If the main header doesn't have the EXTEND keyword, add it now.
  fxhmodify, filename, 'EXTEND', !true, 'The file has extension(s).'

  ;; Make the binary extension header
  fxbhmake, bdr, 1, 'WCS-TAB', 'For storing tabulated WCS tabulations'
  
  fxbaddcol, 1, bdr, wcs_coords,  'POINT+WAVE+TIME',  'Coordinate array'
  fxbaddcol, 2, bdr, findgen(Ntunes)+1, 'HPLN-INDEX', 'Index for helioprojective longitude'
  fxbaddcol, 3, bdr, findgen(Nscans)+1, 'HPLT-INDEX', 'Index for helioprojective latitude'
  fxbaddcol, 4, bdr, findgen(Ntunes)+1, 'WAVE-INDEX', 'Index for tuning/wavelength'
  fxbaddcol, 5, bdr, findgen(Nscans)+1, 'TIME-INDEX', 'Index for repeat/time'
  
  fxbcreate, unit, filename, bdr, extension_no
  fxbwrite, unit, wcs_coords, 1, 1 ; Write coordinates as column 1, row 1
  fxbfinish, unit

end
