; docformat = 'rst'

;+
; Read WCS information from an SST fitscube file.
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
;    filename : in, type=string
; 
;      The name of the file from which the wcs information is to be
;      extracted. 
; 
; 
; :Keywords:
; 
;    coordinates : in, optional, type=structarr
;   
;      An array of structs with the wcs coordinates, { HPLN, HPLT,
;      WAVE, TIME }, as 2 by 2 arrays representing the spatial corners
;      of the FOV.
; 
;    distortions : in, optional
;
;      
; 
; :History:
; 
;   2017-11-03 : MGL. First version.
; 
; 
; 
; 
;-
pro red::fitscube_getwcs, filename $
                          , coordinates = coordinates $
                          , distortions = distortions


  hdr = headfits(filename)
  fxbopen, lun, filename, 'WCS-TAB', bdr

  ;; What coordinates are tabulated?
  ctype = fxpar(hdr, 'CTYPE*')
  TabulateHPLN = strmatch(ctype[0], '*-TAB')
  TabulateHPLT = strmatch(ctype[1], '*-TAB')
  TabulateWAVE = strmatch(ctype[2], '*-TAB')
  TabulateTIME = strmatch(ctype[4], '*-TAB')

  ;; Info from extension header
  tdim1 = fxpar(bdr, 'TDIM1')
  ttype1 = fxpar(bdr, 'TTYPE1')
;  ttype2 = fxpar(bdr, 'TTYPE2')
;  ttype3 = fxpar(bdr, 'TTYPE3')
  fxbread, lun, coord_array, ttype1
;  fxbread, lun, hpln_indx, ttype2
;  fxbread, lun, hplt_indx, ttype3
  fxbclose, lun
  ;; We don't actually need the index arrays. We use the convention
  ;; that hpln_indx and hplt_indx are [1,dim] and the others are
  ;; indgen(dim)+1, where dim is the data array length in the
  ;; respective dimension.
  
  coord_names = strsplit(ttype1, '+', /extract)
  coord_dims = (long(strsplit(strmid(tdim1, 1), ',', /extract)))[1:*]

  if TabulateHPLN then i_hpln = (where(strmatch(coord_names,'HPLN')))[0]
  if TabulateHPLT then i_hplt = (where(strmatch(coord_names,'HPLT')))[0]
  if TabulateWAVE then i_wave = (where(strmatch(coord_names,'WAVE')))[0]
  if TabulateTIME then i_time = (where(strmatch(coord_names,'TIME')))[0]

  if TabulateHPLN then N_hpln = coord_dims[i_hpln] else N_hpln = 1
  if TabulateHPLT then N_hplt = coord_dims[i_hplt] else N_hplt = 1
  if TabulateWAVE then N_wave = coord_dims[i_wave] else N_wave = 1
  if TabulateTIME then N_time = coord_dims[i_time] else N_time = 1

  ;; Reconstruct coord_dims in case of non-tabulated (degenerate)
  ;; dimensions. 
  coord_dims = [N_hpln, N_hplt, N_wave, N_time]
  
  coordinates = replicate({ wave:dblarr(N_hpln, N_hplt) $
                          , hplt:dblarr(N_hpln, N_hplt) $
                          , hpln:dblarr(N_hpln, N_hplt) $
                          , time:dblarr(N_hpln, N_hplt) $
                          }, N_wave, N_time)

  case 1 of
    TabulateWAVE and TabulateTIME : begin
      ;; Both WAVE and TIME
      for j_wave = 0, N_wave-1 do for j_time = 0, N_time-1 do begin
        coordinates[j_wave,j_time].hpln = reform(coord_array[i_hpln,*,*,j_wave,j_time])
        coordinates[j_wave,j_time].hplt = reform(coord_array[i_hplt,*,*,j_wave,j_time])
        coordinates[j_wave,j_time].wave = reform(coord_array[i_wave,*,*,j_wave,j_time])
        coordinates[j_wave,j_time].time = reform(coord_array[i_time,*,*,j_wave,j_time])
      endfor
    end
    TabulateWave : begin
      ;; WAVE, no TIME
      j_time = 0
      for j_wave = 0, N_wave-1 do begin
        coordinates[j_wave,j_time].hpln = reform(coord_array[i_hpln,*,*,j_wave])
        coordinates[j_wave,j_time].hplt = reform(coord_array[i_hplt,*,*,j_wave])
        coordinates[j_wave,j_time].wave = reform(coord_array[i_wave,*,*,j_wave])
     endfor
    end
    TabulateTime : begin
      ;; TIME, no WAVE
      j_wave = 0
      for j_time = 0, N_time-1 do begin
        coordinates[j_wave,j_time].hpln = reform(coord_array[i_hpln,*,*,j_time])
        coordinates[j_wave,j_time].hplt = reform(coord_array[i_hplt,*,*,j_time])
        coordinates[j_wave,j_time].time = reform(coord_array[i_time,*,*,j_time])
      endfor
    end
    else : begin
      ;; Neither WAVE nor TIME
      j_wave = 0
      j_time = 0
      coordinates[j_wave,j_time].hpln = reform(coord_array[i_hpln, *, *])
      coordinates[j_wave,j_time].hplt = reform(coord_array[i_hplt, *, *])
    end
  endcase

  ;; Fill in the non-tabulated coordinates from the main header
  ;; keywords. 
  if ~TabulateWAVE then coordinates.wave += fxpar(hdr, 'CRVAL3') 
  if ~TabulateTIME then coordinates.time += fxpar(hdr, 'CRVAL5') 
  ;; To do this properly, we should check that CUNIT{3,5} are nm and
  ;; s, respectively.

  if arg_present(distortions) then begin

    ;; Assume distortions are always given in lookup form. (We can
    ;; actually check this in the CPDIS* main header keyword.)
    
    dist = mrdfits( filename, 'WCSDVARR', chdr, status = status)

    ;; So far we have only implemented distortions in the wavelength
    ;; coordinate, so we'll assume this is all there could be. But
    ;; returning a struct makes it possible to change this later.

    if status eq 0 then begin
      ;; Leave the distortions keyword undefined if there is no
      ;; WCSDVARR extension.
      dist = reform(wave_distortions $
                    , fxpar(hdr, 'NAXIS1'), fxpar(hdr, 'NAXIS2') $
                    , 1, 1, fxpar(hdr, 'NAXIS5'), /overwrite)
      
      distortions = { WAVE:dist }
    endif

    
  endif
  
  
end
