; docformat = 'rst'

;+
; Add a WCS extension to a fitscube file.
;
; This was inspired by Stein's wcs_crosstabulation.pro
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; 
; :Params:
; 
;    filename : in, type=string
;
;       The name of the file.
;
;    wcs : in, type=struct
;
;       The WCS coordinates hpln, hplt, wave, time as fltarr(2, 2,
;       Ntunes, Nscans) arrays.
;
; :History:
; 
;    2016-03-24 : MGL. First version.
; 
;    2016-08-25 : MGL. Add spatial coordinates.
; 
;    2016-09-28 : MGL. Add a new parameter wcs. 
; 
;    2016-10-11 : MGL. Remove keywords
;                 wcs_{hpln,hplt,wave,time}_coordinate.
; 
;-
pro red::fitscube_addwcs, filename, wcs

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Get the dimensions (using x for hpln and y for hplt).
  xdims = size(wcs.hpln, /dim)
  ydims = size(wcs.hplt, /dim)
  wdims = size(wcs.wave, /dim)
  tdims = size(wcs.time, /dim)

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
    help, wcs
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
  wcs_coords[0, *, *, *, *] = wcs.hpln
  wcs_coords[1, *, *, *, *] = wcs.hplt
  wcs_coords[2, *, *, *, *] = wcs.wave
  wcs_coords[3, *, *, *, *] = wcs.time

  ;; If the main header doesn't have the EXTEND keyword, add it now.
  fxhmodify, filename, 'EXTEND', !true, 'The file has extension(s).'

  ;; Make the binary extension header
  fxbhmake, bdr, 1, 'WCS-TAB', 'For storing tabulated WCS tabulations'
  
  fxbaddcol, 1, bdr, wcs_coords, 'POINT+WAVE+TIME', 'Coordinate array'
  fxbaddcol, 2, bdr, findgen(Ntunes)+1, 'HPLN-INDEX', 'Index for helioprojective longitude'
  fxbaddcol, 3, bdr, findgen(Nscans)+1, 'HPLT-INDEX', 'Index for helioprojective latitude'
  fxbaddcol, 4, bdr, findgen(Ntunes)+1, 'WAVE-INDEX', 'Index for tuning/wavelength'
  fxbaddcol, 5, bdr, findgen(Nscans)+1, 'TIME-INDEX', 'Index for repeat/time'
  
  fxbcreate, unit, filename, bdr, extension_no
  fxbwrite, unit, wcs_coords, 1, 1 ; Write coordinates as column 1, row 1
  fxbfinish, unit

end
