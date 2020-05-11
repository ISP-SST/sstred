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
; :Keywords:
;
;    csyer_spatial_value : in, optional, type=float, default=60.
;
;       Upper limit of systematic error in the spatial cordinates (in
;       arcsec), used to set the value of FITS header keywords CSYER1
;       and CSYER2. The default 60 arcsec is for when the rotation is
;       unspecified.
;
;    csyer_spatial_comment : in, optional, type=string
;
;       The comment for FITS header keywords CSYER1 and CSYER2.
;
;    dimensions : in, type=intarr
;
;       The dimensions of the data cube.
;
;    no_extension : in, optional, type=boolean
;
;       Only add the header keywords, not the extension.
;
;    update : in, optional, type=boolean
;
;       Update an existing WCS extension.
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
;    2017-10-18 : MGL. WCS related file header keywords moved here
;                 from fitscube_initialize method. New keyword
;                 dimensions. 
;
;    2018-01-16 : MGL. Allow for degenerate trailing dimensions of
;                 scan-only cubes.
;
;    2018-11-29 : MGL. New keyword no_extension.
;
;    2019-09-13 : MGL. A version that is not a class method.
;
;    2020-03-23 : MGL. New keyword update.
;
;    2020-03-27 : MGL. New keywords csyer_spatial_value and
;                 csyer_spatial_comment.
;
;-
pro red_fitscube_addwcs, filename, wcs $
                         , csyer_spatial_value = csyer_spatial_value $
                         , csyer_spatial_comment = csyer_spatial_comment $
                         , dimensions = dimensions $
                         , no_extension = no_extension $
                         , update = update

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  if n_elements(csyer_spatial_value) eq 0 then begin
    csyer_spatial_value = 120.  ; 2 arc minutes
    csyer_spatial_comment = '[arcsec] Orientation unknown'
  endif 
  
  ;; The code leading up to the definition of WriteTimeIndex,
  ;; WriteWaveIndex, TabulateWave, TabulateTime has to match the code
  ;; in red_fitscube_addwcsheader. 
  
  ;; Get the dimensions of the coordinate vectors (using x for hpln
  ;; and y for hplt). 
  xdims = size(wcs.hpln, /dim)
  ydims = size(wcs.hplt, /dim)
  wdims = size(wcs.wave, /dim)
  tdims = size(wcs.time, /dim)

  Nxdims = n_elements(xdims)
  Nydims = n_elements(ydims)
  Nwdims = n_elements(wdims)
  Ntdims = n_elements(tdims)

  ;; Degenerate trailing dimensions?
  if Nxdims eq 3 then begin
    xdims = [xdims, 1]
    Nxdims = n_elements(xdims)
  endif
  if Nydims eq 3 then begin
    ydims = [ydims, 1]
    Nydims = n_elements(ydims)
  endif
  if Nwdims eq 3 then begin
    wdims = [wdims, 1]
    Nwdims = n_elements(wdims)
  endif
  if Ntdims eq 3 then begin
    tdims = [tdims, 1]
    Ntdims = n_elements(tdims)
  endif
  
  indx0 = where([Nxdims, Nydims, Nwdims, Ntdims] eq 0, N0)
  indx2 = where([Nxdims, Nydims, Nwdims, Ntdims] eq 2, N2)
  indx3 = where([Nxdims, Nydims, Nwdims, Ntdims] eq 3, N3)
  indx4 = where([Nxdims, Nydims, Nwdims, Ntdims] eq 4, N4)
  
  ;; At least one of the arrays have to be non-scalar. Any non-scalars
  ;; must have the same dimensions: (2 x 2 x Ntunes x Nscans).

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

  ;; Finally, the dimensions ("sizes") of the tabulated coordinate
  ;; array.
  Sx = dims[0]
  Sy = dims[1]
  if n_elements(dims) ge 3 then Sw = dims[2] else Sw = 1
  if n_elements(dims) ge 4 then St = dims[3] else St = 1

  ;; Get the dimensions of the data cube. This kind of FITS cube is
  ;; always five-dimensional, but dimensions can be degenerate.
  if n_elements(dimensions) ge 1 then Nx      = long(dimensions[0]) else Nx      = Nxdims 
  if n_elements(dimensions) ge 2 then Ny      = long(dimensions[1]) else Ny      = Nydims 
  if n_elements(dimensions) ge 3 then Ntuning = long(dimensions[2]) else Ntuning = Nwdims 
  if n_elements(dimensions) ge 4 then Nstokes = long(dimensions[3])
  if n_elements(dimensions) ge 5 then Nscans  = long(dimensions[4]) else Nscans  = Ntdims 

  ;; We don't want to tabulate WAVE and/or TIME if the dimension is
  ;; just 1 long, then no tabulation needed.
  TabulateWave = Sw gt 1
  TabulateTime = St gt 1 || Sw gt 1
  ;; We don't want to write the index arrays for WAVE and TIME if the
  ;; coordinate is not tabulated, or there is exactly one coordinate
  ;; tabulated for each pixel because then the default index array
  ;; works.
  WriteWaveIndex = TabulateWave && (Sw ne Ntuning)
  WriteTimeIndex = TabulateTime && (St ne Nscans || Sw ne Ntuning)
  
  ;; An array to store the coordinates: wcs_coordinates. First
  ;; dimension is the number of tabulated coordinates, then follow the
  ;; numbers of each kind of tabulated coordinate.
  arraydims = [2 + TabulateWave + TabulateTime, Sx, Sy]
  if TabulateWave then red_append, arraydims, Sw
  if TabulateTime then red_append, arraydims, St
  wcs_coords = dblarr(arraydims)
  colno = 0

  case 1 of
    TabulateWave && TabulateTime : begin
      ;; Both WAVE and TIME
      wcs_coords[colno++, *, *, *, *] = wcs.hpln
      wcs_coords[colno++, *, *, *, *] = wcs.hplt
      wcs_coords[colno++, *, *, *, *] = wcs.wave
      wcs_coords[colno++, *, *, *, *] = wcs.time
    end
    TabulateWave || TabulateTime : begin
      ;; WAVE or TIME but not both
      wcs_coords[colno++, *, *, *] = wcs.hpln
      wcs_coords[colno++, *, *, *] = wcs.hplt
      if TabulateWave then wcs_coords[colno++, *, *, *] = wcs.wave
      if TabulateTime then wcs_coords[colno++, *, *, *] = wcs.time      
    end
    else : begin
      ;; Neither WAVE nor TIME
      wcs_coords[colno++, *, *] = wcs.hpln
      wcs_coords[colno++, *, *] = wcs.hplt
    end
  endcase

  hdr = headfits(filename)
  
  if keyword_set(update) then begin

    if n_elements(csyer_spatial_value) gt 0 then begin
      ;; Change the header keywords for systematic errors.
      red_fitsaddkeyword, hdr, anchor = 'CDELT1' $
                          , 'CSYER1', csyer_spatial_value, csyer_spatial_comment
      red_fitsaddkeyword, hdr, anchor = 'CDELT2' $
                          , 'CSYER2', csyer_spatial_value, csyer_spatial_comment
    endif
    modfits, filename, 0, hdr

    ;; Open an existing WCS-TAB extension for updating.
    
    fxbopen, bunit, filename, 'WCS-TAB', bhdr, access = 'RW'

  endif else begin

    ;; Construct the TTYPE1 keyword
    ttypes = ['HPLN', 'HPLT']
    if TabulateWave then red_append, ttypes, ['WAVE']
    if TabulateTime then red_append, ttypes, ['TIME']
    ttype = strjoin(ttypes, '+')


    ;; Modify the main header. ---------------------------------------------------------------
    red_fitscube_addwcsheader, hdr, wcs $
                               , dimensions = dimensions $
                               , csyer_spatial_value = csyer_spatial_value $
                               , csyer_spatial_comment = csyer_spatial_comment 

    ;; If the new header needs another block of 80 lines to fit, this
    ;; operation may take a long time because the whole data part of
    ;; the file needs to be moved on disk. Can we change the way the
    ;; wcs info is handled, so the header is complete before the data
    ;; part is written? Can the extension below be added to the file
    ;; before the data part?

    print
    print, inam+' : Writing a new file header with WCS information.'
    print
    print, inam+' : If this header needs another 2880-byte block, the data'
    print, inam+' : parts of the file have to be moved on the disk. This may'
    print, inam+' : take a long time, depending on the size of the file.'
    print
    tic
    modfits, filename, 0, hdr
    toc
    print, inam+' : Done writing the new header!'
    print

    ;; Return now if we don't actually want to add the extension.
    if keyword_set(no_extension) then return

    ;; Make the binary extension. ------------------------------------------------------------

    ;; Make the extension header
    fxbhmake, bdr, 1, 'WCS-TAB', 'For storing tabulated WCS coordinates'
    colno = 1
    fxbaddcol, colno++, bdr, wcs_coords, ttype, 'Coordinate array'
    fxbaddcol, colno++, bdr, [1., Nx], 'HPLN-INDEX', 'Index for helioprojective longitude'
    fxbaddcol, colno++, bdr, [1., Ny], 'HPLT-INDEX', 'Index for helioprojective latitude'
    if WriteWaveIndex then fxbaddcol, colno++, bdr, findgen(Ntuning)+1., 'WAVE-INDEX', 'Index for tuning/wavelength'
    if WriteTimeIndex then fxbaddcol, colno++, bdr, findgen(Nscans)+1.,  'TIME-INDEX', 'Index for repeat/time'

    ;; Create the extension and leave the file open
    fxbcreate, bunit, filename, bdr, extension_no
    
  endelse

  ;; Write the extension
  colno = 1
  rowno = 1
  fxbwrite, bunit, wcs_coords, colno++, rowno                            ; Write coordinates
  fxbwrite, bunit, [1., Nx],   colno++, rowno                            ; Write hpln index array
  fxbwrite, bunit, [1., Ny],   colno++, rowno                            ; Write hplt index array
  if WriteWaveIndex then fxbwrite, bunit, findgen(Sw)+1., colno++, rowno ; Write wave index array
  if WriteTimeIndex then fxbwrite, bunit, findgen(St)+1., colno++, rowno ; Write time index array
  fxbfinish, bunit

end
