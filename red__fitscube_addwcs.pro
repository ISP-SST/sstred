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
;-
pro red::fitscube_addwcs, filename, wcs, dimensions = dimensions

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

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
;  if N0+N4 ne 4 || N4 eq 0 then begin
;    print, inam + ' : At least one of the arrays has the wrong dimensions.'
;    help, wcs
;    stop
;    retall
;  endif

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

  ;; Construct the TTYPE1 keyword
  ttypes = ['HPLN', 'HPLT']
  if TabulateWave then red_append, ttypes, ['WAVE']
  if TabulateTime then red_append, ttypes, ['TIME']
  ttype = strjoin(ttypes, '+')


  ;; TIME reference value, all times are seconds since midnight. One
  ;; should really check for the existence of MJDREF and JDREF but
  ;; within the pipeline we can be sure we don't use them.
  dateref = self.isodate+'T00:00:00.000000' ; Midnight
  
  ;; Modify the main header. ---------------------------------------------------------------
  
  hdr = headfits(filename)
  Naxis = fxpar(hdr,'NAXIS')

  ;; If the main header doesn't have the EXTEND keyword, add it now.
  red_fitsaddkeyword, hdr, 'EXTEND', !true, 'The file has extension(s).'

  ;; DATEREF is added to the WCS time coordinate
  red_fitsaddkeyword, hdr, 'DATEREF', dateref, 'Reference time in ISO-8601', after = 'DATE'

  anchor = 'FILENAME'
  
  red_fitsaddkeyword, anchor = anchor, hdr, 'PC1_1', 1.0, 'No rotations' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'PC2_2', 1.0, 'No rotations' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'PC3_3', 1.0, 'No rotations' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'PC4_4', 1.0, 'No rotations' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'PC5_5', 1.0, 'No rotations'

  ;; The header could be a copy from a file that has WCS in it, so
  ;; better remove old WCS related keywords.
  keywords = strmid(hdr, 0, 8)
  for iax = 0, Naxis-1 do begin
    ckeywords = keywords[where(strmatch(keywords.trim(),'C*'+strtrim(iax+1, 2)), Nc)]
    pkeywords = keywords[where(strmatch(keywords.trim(),'P[SV]*'+strtrim(iax+1, 2)+'_*'), Np)]
    for ikey = 0, Nc-1 do red_fitsdelkeyword, hdr, ckeywords[ikey]
    for ikey = 0, Np-1 do red_fitsdelkeyword, hdr, pkeywords[ikey]
  endfor                        ; iax

  ;; Now add the new keywords
  coordno = 1                   ; Tabulated coordinate number, for PVi_3 keywords
  
  ;; First spatial dimension, corner coordinates always tabulated 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CTYPE1', 'HPLN-TAB', 'SOLAR X'
  red_fitsaddkeyword, anchor = anchor, hdr, 'CUNIT1', 'arcsec', 'Unit along axis 1'
  red_fitsaddkeyword, anchor = anchor, hdr, 'CNAME1', 'Spatial X'
  red_fitsaddkeyword, anchor = anchor, hdr, 'PS1_0', 'WCS-TAB', 'EXTNAME; EXTVER=EXTLEVEL=1 is default'
  red_fitsaddkeyword, anchor = anchor, hdr, 'PS1_1', ttype, 'TTYPE for column w/coordinates'
  red_fitsaddkeyword, anchor = anchor, hdr, 'PS1_2', 'HPLN-INDEX', 'TTYPE for INDEX'
  red_fitsaddkeyword, anchor = anchor, hdr, 'PV1_3', coordno++, 'Coord. 1 tabulated coordinate number' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CRPIX1', 0, 'Unity transform' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CRVAL1', 0, 'Unity transform' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CDELT1', 1, 'Unity transform' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CSYER1', 60, 'Orientation unknown' 

  ;; Second spatial dimension, corner coordinates always tabulated 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CTYPE2', 'HPLT-TAB', 'SOLAR Y' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CUNIT2', 'arcsec', 'Unit along axis 2' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CNAME2', 'Spatial Y' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'PS2_0', 'WCS-TAB', 'EXTNAME; EXTVER=EXTLEVEL=1 is default' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'PS2_1', ttype, 'TTYPE for column w/coordinates' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'PS2_2', 'HPLT-INDEX', 'TTYPE for INDEX' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'PV2_3', coordno++, 'Coord. 2 tabulated coordinate number' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CRPIX2', 0, 'Unity transform' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CRVAL2', 0, 'Unity transform' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CDELT2', 1, 'Unity transform' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CSYER2', 60, 'Orientation unknown' 

  ;; Tuning, tabulated wavelength, tabulated for actual scans but not
  ;; for, e.g., wideband cubes.
  if TabulateWave then begin
    red_fitsaddkeyword, anchor = anchor, hdr, 'CTYPE3', 'WAVE-TAB', 'Wavelength, function of tuning and scan number' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'CNAME3', 'Wavelength' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'CUNIT3', 'nm', 'Wavelength unit, tabulated for dim. 3 and 5' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'PS3_0', 'WCS-TAB', 'EXTNAME; EXTVER=EXTLEVEL=1 is default' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'PS3_1', ttype, 'TTYPE for column w/coordinates' 
    if WriteWaveIndex then red_fitsaddkeyword, anchor = anchor, hdr, 'PS3_2', 'WAVE-INDEX', 'TTYPE for INDEX'
    red_fitsaddkeyword, anchor = anchor, hdr, 'PV3_3', coordno++, 'Coord. 3 tabulated coordinate number' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'CRPIX3', 0, 'Unity transform' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'CRVAL3', 0, 'Unity transform' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'CDELT3', 1, 'Unity transform' 
  endif else begin
    red_fitsaddkeyword, anchor = anchor, hdr, 'CTYPE3', 'WAVE'
    red_fitsaddkeyword, anchor = anchor, hdr, 'CNAME3', 'Wavelength' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'CUNIT3', 'nm', 'Wavelength unit'
    red_fitsaddkeyword, anchor = anchor, hdr, 'CRVAL3', wcs.wave[0], 'Just a single wavelength'  , /force
  endelse
  
  ;; Stokes
  red_fitsaddkeyword, anchor = anchor, hdr, 'CTYPE4', 'STOKES', 'Stokes vector [I,Q,U,V]'
  red_fitsaddkeyword, anchor = anchor, hdr, 'CRPIX4', 1, 'Index of Stokes components in pixel 1' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CRVAL4', 1, 'The first Stokes index is 1' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CDELT4', 1, 'Stokes indices [1,2,3,4] --> [I,Q,U,V]' 
    
  ;; Scan number = repetition = major time dimension. But time varies
  ;; during scans as well, so we can only avoid tabulating if both St
  ;; and Sw are unity.
  if TabulateTime || TabulateWave then begin
    red_fitsaddkeyword, anchor = anchor, hdr, 'CTYPE5', 'UTC--TAB', 'Time, function of tuning and scan number' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'CNAME5', 'Time since DATEREF, increases with dim. 3 and 5' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'CUNIT5', 's' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'PS5_0', 'WCS-TAB',   'EXTNAME; EXTVER=EXTLEVEL=1 is default' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'PS5_1', ttype, 'TTYPE for column w/coordinates' 
    if WriteTimeIndex then red_fitsaddkeyword, anchor = anchor, hdr, 'PS5_2', 'TIME-INDEX','TTYPE for INDEX' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'PV5_3', coordno++, 'Coord. 5 tabulated coordinate number' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'CRPIX5', 0, 'Unity transform' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'CRVAL5', 0, 'Unity transform' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'CDELT5', 1, 'Unity transform' 
  endif else begin
    red_fitsaddkeyword, anchor = anchor, hdr, 'CTYPE5', 'UTC', 'Time'  
    red_fitsaddkeyword, anchor = anchor, hdr, 'CNAME5', 'Time since DATEREF, increases with dim. 3 and 5' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'CUNIT5', 's' 
    red_fitsaddkeyword, anchor = anchor, hdr, 'CRVAL5', wcs.time[0], 'Just a single time'  
  endelse

  ;; If the new header needs another block of 80 lines to fit, this
  ;; operation may take a long time because the whole data part of the
  ;; file needs to be moved on disk. Can we change the way the wcs
  ;; info is handled, so the header is complete before the data part
  ;; is written? Can the extension below be added to the file before
  ;; the data part?

  print
  print, inam+' : Writing a new file header with WCS information.'
  print
  print, inam+' : If this header needs another 2880-byte block, the data'
  print, inam+' : parts of the file have to be moved on the disk. This may'
  print, inam+' : take a long time, depending on the size of the file.'
  print
  modfits, filename, 0, hdr
  print, inam+' : Done writing the new header!'
  print
  
  ;; Make the binary extension. ------------------------------------------------------------

  ;; Make the extension header
  fxbhmake, bdr, 1, 'WCS-TAB', 'For storing tabulated WCS coordinates'
  colno = 1
  fxbaddcol, colno++, bdr, wcs_coords, ttype, 'Coordinate array'
  fxbaddcol, colno++, bdr, [1., Nx], 'HPLN-INDEX', 'Index for helioprojective longitude'
  fxbaddcol, colno++, bdr, [1., Ny], 'HPLT-INDEX', 'Index for helioprojective latitude'
  if WriteWaveIndex then fxbaddcol, colno++, bdr, findgen(Ntuning)+1., 'WAVE-INDEX', 'Index for tuning/wavelength'
  if WriteTimeIndex then fxbaddcol, colno++, bdr, findgen(Nscans)+1.,  'TIME-INDEX', 'Index for repeat/time'

  ;; Write the extension
  fxbcreate, unit, filename, bdr, extension_no
  colno = 1
  rowno = 1
  fxbwrite, unit, wcs_coords, colno++, rowno                     ; Write coordinates
  fxbwrite, unit, [1., Nx],   colno++, rowno                     ; Write hpln index array
  fxbwrite, unit, [1., Ny],   colno++, rowno                     ; Write hplt index array
  if WriteWaveIndex then fxbwrite, unit, findgen(Sw)+1., colno++, rowno ; Write wave index array
  if WriteTimeIndex then fxbwrite, unit, findgen(St)+1., colno++, rowno ; Write time index array
  fxbfinish, unit

end
