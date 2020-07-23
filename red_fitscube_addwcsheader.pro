; docformat = 'rst'

;+
; Add WCS keywords to a FITS header.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    hdr : in, type=strarr
; 
;       The header in which to add the WCS keywords.
;
;    wcs : in, type=struct
;
;       The WCS coordinates hpln, hplt, wave, time as fltarr(2, 2,
;       Ntunes, Nscans) arrays.
; 
; :Keywords:
;
;    dimensions : in, type=array
;   
;       The dimensions of the fitscube.
; 
; 
; :History:
; 
;    2018-11-26 : MGL. First version.
;
;    2020-04-30 : MGL. New keywords csyer_spatial_value,
;                 csyer_spatial_comment.
; 
;-
pro red_fitscube_addwcsheader, hdr, wcs $
                               , csyer_spatial_value = csyer_spatial_value $
                               , csyer_spatial_comment = csyer_spatial_comment $
                               , dimensions = dimensions
  
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  if n_elements(csyer_spatial_value) eq 0 then begin
    csyer_spatial_value = 120.  ; 2 arc minutes
    csyer_spatial_comment = '[arcsec] Orientation unknown'
  endif 

  
  ;; The code leading up to the definition of WriteTimeIndex,
  ;; WriteWaveIndex, TabulateWave, TabulateTime has to match the code
  ;; in red::fitsube_addwcs. 
  
  Naxis = fxpar(hdr,'NAXIS')
  
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
  
;  print, inam+' : WriteWaveIndex, WriteTimeIndex', WriteWaveIndex, WriteTimeIndex
;  print, dimensions
;  print, dims
;  help, Sx, Sy, Sw, St, Nx, Ny, Ntuning, Nstokes, Nscans
;
;
;  stop
  
;  ;; An array to store the coordinates: wcs_coordinates. First
;  ;; dimension is the number of tabulated coordinates, then follow the
;  ;; numbers of each kind of tabulated coordinate.
;  arraydims = [2 + TabulateWave + TabulateTime, Sx, Sy]
;  if TabulateWave then red_append, arraydims, Sw
;  if TabulateTime then red_append, arraydims, St
;  wcs_coords = dblarr(arraydims)
;  colno = 0
;
;  case 1 of
;    TabulateWave && TabulateTime : begin
;      ;; Both WAVE and TIME
;      wcs_coords[colno++, *, *, *, *] = wcs.hpln
;      wcs_coords[colno++, *, *, *, *] = wcs.hplt
;      wcs_coords[colno++, *, *, *, *] = wcs.wave
;      wcs_coords[colno++, *, *, *, *] = wcs.time
;    end
;    TabulateWave || TabulateTime : begin
;      ;; WAVE or TIME but not both
;      wcs_coords[colno++, *, *, *] = wcs.hpln
;      wcs_coords[colno++, *, *, *] = wcs.hplt
;      if TabulateWave then wcs_coords[colno++, *, *, *] = wcs.wave
;      if TabulateTime then wcs_coords[colno++, *, *, *] = wcs.time      
;    end
;    else : begin
;      ;; Neither WAVE nor TIME
;      wcs_coords[colno++, *, *] = wcs.hpln
;      wcs_coords[colno++, *, *] = wcs.hplt
;    end
;  endcase

  ;; Construct the TTYPE1 keyword
  ttypes = ['HPLN', 'HPLT']
  if TabulateWave then red_append, ttypes, ['WAVE']
  if TabulateTime then red_append, ttypes, ['TIME']
  ttype = strjoin(ttypes, '+')


  ;; Modify the main header. ---------------------------------------------------------------
  
  
  ;; If the main header doesn't have the EXTEND keyword, add it now.
  red_fitsaddkeyword, hdr, 'EXTEND', !true, 'The file has extension(s).'

  anchor = 'FILENAME'

  red_fitsaddkeyword, anchor = anchor, hdr, 'PC1_1', 1.0, 'No rotations' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'PC2_2', 1.0, 'No rotations' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'PC3_3', 1.0, 'No rotations' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'PC4_4', 1.0, 'No rotations' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'PC5_5', 1.0, 'No rotations'

  ;; The header could be a copy from a file that has WCS in it, so
  ;; better remove old WCS related keywords. (Add to delkeywords if
  ;; more WCS keywords are going to be used.)
;  keywords = strmid(hdr, 0, 8)
;  for iax = 0, Naxis-1 do begin
;    ckeywords = keywords[where(strmatch(keywords.trim(),'C*'+strtrim(iax+1, 2)), Nc)]
;    pkeywords = keywords[where(strmatch(keywords.trim(),'P[SV]*'+strtrim(iax+1, 2)+'_*'), Np)]
;    for ikey = 0, Nc-1 do red_fitsdelkeyword, hdr, ckeywords[ikey]
;    for ikey = 0, Np-1 do red_fitsdelkeyword, hdr, pkeywords[ikey]
;  endfor                        ; iax
;
  red_append, delkeywords, ['CDELT1', 'CDELT2', 'CDELT3', 'CDELT4', 'CDELT5']
  red_append, delkeywords, ['CNAME1', 'CNAME2', 'CNAME3', 'CNAME3', 'CNAME5']
  red_append, delkeywords, ['CRPIX1', 'CRPIX2', 'CRPIX3', 'CRPIX4', 'CRPIX5']
  red_append, delkeywords, ['CRVAL1', 'CRVAL2', 'CRVAL3', 'CRVAL4', 'CRVAL5']
  red_append, delkeywords, ['CSYER1', 'CSYER2', 'CSYER3', 'CSYER4', 'CSYER5']
  red_append, delkeywords, ['CTYPE1', 'CTYPE2', 'CTYPE3', 'CTYPE4', 'CTYPE5']
  red_append, delkeywords, ['CUNIT1', 'CUNIT2', 'CUNIT3', 'CUNIT3', 'CUNIT5']
  red_append, delkeywords, ['PS1_0', 'PS1_1', 'PS1_2']
  red_append, delkeywords, ['PS2_0', 'PS2_1', 'PS2_2']
  red_append, delkeywords, ['PS3_0', 'PS3_1', 'PS3_2']
  red_append, delkeywords, ['PS5_0', 'PS5_1', 'PS5_2']
  red_append, delkeywords, ['PV1_3', 'PV2_3', 'PV3_3']
  red_append, delkeywords, ['PV5_3']
  for ikey = 0, n_elements(delkeywords)-1 do red_fitsdelkeyword, hdr, delkeywords[ikey]


  
  ;; Now add the new keywords
  coordno = 1                   ; Tabulated coordinate number, for PVi_3 keywords
  
  ;; First spatial dimension, corner coordinates always tabulated 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CTYPE1', 'HPLN-TAB', 'Helioprojective longitude'
  red_fitsaddkeyword, anchor = anchor, hdr, 'CUNIT1', 'arcsec', 'Unit along axis 1'
  red_fitsaddkeyword, anchor = anchor, hdr, 'CNAME1', 'Spatial X'
  red_fitsaddkeyword, anchor = anchor, hdr, 'PS1_0', 'WCS-TAB', 'EXTNAME; EXTVER=EXTLEVEL=1 is default'
  red_fitsaddkeyword, anchor = anchor, hdr, 'PS1_1', ttype, 'TTYPE for column w/coordinates'
  red_fitsaddkeyword, anchor = anchor, hdr, 'PS1_2', 'HPLN-INDEX', 'TTYPE for INDEX'
  red_fitsaddkeyword, anchor = anchor, hdr, 'PV1_3', coordno++, 'Coord. 1 tabulated coordinate number' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CRPIX1', 0, 'Unity transform' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CRVAL1', 0, 'Unity transform' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CDELT1', 1, 'Unity transform' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CSYER1', csyer_spatial_value, csyer_spatial_comment

  ;; Second spatial dimension, corner coordinates always tabulated 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CTYPE2', 'HPLT-TAB', 'Helioprojective latitude' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CUNIT2', 'arcsec', 'Unit along axis 2' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CNAME2', 'Spatial Y' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'PS2_0', 'WCS-TAB', 'EXTNAME; EXTVER=EXTLEVEL=1 is default' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'PS2_1', ttype, 'TTYPE for column w/coordinates' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'PS2_2', 'HPLT-INDEX', 'TTYPE for INDEX' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'PV2_3', coordno++, 'Coord. 2 tabulated coordinate number' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CRPIX2', 0, 'Unity transform' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CRVAL2', 0, 'Unity transform' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CDELT2', 1, 'Unity transform' 
  red_fitsaddkeyword, anchor = anchor, hdr, 'CSYER2', csyer_spatial_value, csyer_spatial_comment

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

end
