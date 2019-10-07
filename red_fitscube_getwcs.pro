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
;    distortions : out, optional
;
;      
; 
; :History:
; 
;   2017-11-03 : MGL. First version.
; 
;   2018-12-05 : MGL. Change it from a class method to a regular
;                procedure. 
; 
;   2019-09-25 : MGL. Implement reading multiple wavelength
;                distortions (cavity maps).
; 
;-
pro red_fitscube_getwcs, filename $
                         , coordinates = coordinates $
                         , distortions = distortions

  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

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
  if TabulateTIME then begin
    i_time = (where(strmatch(coord_names,'TIME')))[0]
    if i_time eq 3 and n_elements(coord_dims) eq 3 then coord_dims = [coord_dims, 1]
  endif


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

    ;; Read hierarch keyword dw3
    dw3 = red_fitsgetkeyword(hdr, 'DW3')
    Ndw3 = n_elements(dw3)
    if Ndw3 eq 0 then begin
      print, inam + ' : There is no DW3 keyword.'
      print, inam + ' : No distortions returned.'
      return
    endif
    dw3_keywords = strarr(Ndw3)
    for idw3 = 0, Ndw3-1 do dw3_keywords[idw3] = (dw3[idw3])[0]
    ;; APPLY should be the last DW3 keyword for each distortion, and
    ;; we usually put NAME first.
    name_pos = where(strmatch(dw3_keywords, 'NAME'), Nname) 
    apply_pos = where(strmatch(dw3_keywords, 'APPLY'), Napply)

    ;; Get the scales and offsets for the tuning coordinate
    offsets = fltarr(Napply)
    scales = fltarr(Napply) 
    for i = 0, Napply-1 do begin
      dw3_sub_keywords = dw3_keywords[name_pos[i]:apply_pos[i]]
      dw3_sub = dw3[name_pos[i]:apply_pos[i]]
      pos = where(dw3_sub_keywords eq 'EXTVER')
      extver = (dw3_sub[pos[0]])(1)
      scale_pos = where(dw3_sub_keywords eq 'SCALE3', Nscale)
      offset_pos = where(dw3_sub_keywords eq 'OFFSET3', Noffset)
      if Nscale eq 0 or Noffset eq 0 then begin
        ;; Old WCSDVARR format without SCALE3 and OFFSET3
        indx = indgen( fxpar(hdr, 'NAXIS3') )
        if n_elements(indx) eq 1 then begin
          ;; A single-tuningpoint cube, probably a Stokes cube.
          scales[i]  = 1.
          offsets[i] = 1.
        endif else begin
          eps = 1e-3
          scales[i] = (1. - 2.*eps) / ((indx[-1]+0.5) - (indx[0]-0.5))
          offsets[i] = ( (indx[-1]+0.5)*(.5+eps) - (1.5-eps)*(indx[0]-0.5) ) / (1. - 2.*eps)
        endelse
      endif else begin
        scales[extver-1] = (dw3_sub[scale_pos[0]])(1)
        offsets[extver-1] = (dw3_sub[offset_pos[0]])(1)
      endelse
    endfor                      ; i
    
    ;; Are there distortions extensions? They should be in an
    ;; WCSDVARR extension.
    fits_open, filename, fcb
    free_lun, fcb.unit
    dindx = where(fcb.extname eq 'WCSDVARR', Ndist)
    if Ndist eq 0 then begin
      print, inam + ' : There is no WCSDVARR extension.'
      print, inam + ' : No distortions returned.'
      return
    endif

    if Ndist ne Napply then stop ; Should match!
    
    for idist = 0, Ndist-1 do begin

      wcsdvarr = mrdfits( filename, dindx[idist], chdr, status = status, /silent)

      if status ne 0 then  begin
        print, inam + ' : There was an error reading WCSDVARR extension #'+strtrim(idist, 2)
        print, inam + ' : No distortions returned.'
        undefine, distortions
        return
      endif
      
      if n_elements(distortions) eq 0 then begin
        ;; So far we have only implemented distortions in the
        ;; wavelength coordinate, so we'll assume this is all there
        ;; could be. But returning a struct makes it possible to
        ;; change this later.
        distortions = replicate({ wave:wcsdvarr $
                                  , tun_index:'' $
                                }, Ndist)
      endif

      ;; index=1 for tunings that are affected by this distortion
      dist_index = round((indgen( fxpar(hdr, 'NAXIS3') ) + offsets[idist])*scales[idist] )
      ;; The indices of the tunings that are affected
      tun_index = red_collapserange(where(dist_index eq 1), ld='', rd='')

      ;; Populate the struct.
      distortions[idist].wave      = wcsdvarr
      distortions[idist].tun_index = tun_index

    endfor                      ; idist 
    
  endif 


end
