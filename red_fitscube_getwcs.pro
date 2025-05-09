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
;    distortions : out, optional, type=structarr
;
;      WCS coordinate distortions. Typically cavity maps in the form
;      of distortions to the wavelength coordinate. The number of
;      elements in the array is the number of prefilters requiring
;      separate cavity maps. Each struct has two members: WAVE is the
;      wavelength distortions for each spatial coordinate and scan
;      index. TUN_INDEX is a list (a comma- and dash-separated string,
;      like '0-20') of tuning indices for which the distortions apply.
;      
;    iscan : in, optional, type=integer  
;      
;      Return only wcs info for this scan index. 
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
;   2021-05-03 : MGL. Read the record-valued version of DW3, not the
;                HIERARCH version.
; 
;   2023-03-22 : MGL. New keyword iscan.
; 
;-
pro red_fitscube_getwcs, filename $
                         , coordinates = coordinates $
                         , distortions = distortions $
                         , iscan = iscan

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)                    

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

  if n_elements(iscan) gt 0 && iscan ge N_time then begin
    print, inam + ' : iscan too large: ', iscan, N_time
    stop 
  endif
    
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

  if n_elements(iscan) gt 0 then begin
    ;; The amount of data in coordinates is small so we read them all
    ;; and then select only iscan.
    coordinates = coordinates[*, iscan]
  endif
  
  if arg_present(distortions) then begin

    undefine, distortions
    
    ;; Read hierarch keyword dw3
    dw3 = red_fitsgetkeyword(hdr, 'DW3', field_specifiers = dw3_keywords, count = Ndw3)
    ;;Ndw3 = n_elements(dw3)
    if Ndw3 eq 0 then begin
      print, inam + ' : There is no DW3 keyword.'
      print, inam + ' : No distortions returned.'
      return
    endif

    ;; APPLY should be the last DW3 keyword for each distortion
    last_pos = where(strmatch(dw3_keywords, 'APPLY'), Napply)
    first_pos = [0]
    if n_elements(last_pos) gt 1 then red_append, first_pos, last_pos[0:-2]+1

    ;; Get the scales and offsets for the tuning coordinate
    offsets = fltarr(Napply)
    scales = fltarr(Napply) 
    for i = 0, Napply-1 do begin
      dw3_sub_keywords = dw3_keywords[first_pos[i]:last_pos[i]]
      dw3_sub = dw3[first_pos[i]:last_pos[i]]
      pos = where(dw3_sub_keywords eq 'EXTVER')
      extver = dw3_sub[pos[0]]
      scale_pos = where(dw3_sub_keywords eq 'SCALE3' or dw3_sub_keywords eq 'SCALE.3', Nscale)
      offset_pos = where(dw3_sub_keywords eq 'OFFSET3' or dw3_sub_keywords eq 'OFFSET.3', Noffset)
      if Nscale eq 0 or Noffset eq 0 then begin
        ;; Old WCSDVARR format without SCALE.3 and OFFSET.3
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
        scales[extver-1] = dw3_sub[scale_pos[0]]
        offsets[extver-1] = dw3_sub[offset_pos[0]]
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

      ;;wcsdvarr = mrdfits( filename, dindx[idist], chdr, status =
      ;;status, /silent)
      ;; If iscan is given, read data for just for that scan.
      wcsdvhdr = headfits( filename, ext = 'WCSDVARR')
      wcsdims = fxpar(wcsdvhdr, 'NAXIS*')
      if round(product(wcsdims[2:4]) eq 1) then begin
        if n_elements(iscan) eq 0 || iscan eq 0 then begin
          wcsdvarr = readfits( filename, exten_no = dindx[idist] )
        endif else stop
      endif else begin
        wcsdvarr = readfits( filename, exten_no = dindx[idist], nslice = iscan)
      endelse
              
      ;;if status ne 0 then  begin
      ;;  print, inam + ' : There was an error reading WCSDVARR extension #'+strtrim(idist, 2)
      ;;  print, inam + ' : No distortions returned.'
      ;;  undefine, distortions
      ;;  return
      ;;endif
      
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


cd, '/scratch/mats/2016.09.19/CHROMIS-jan19'
fname = 'cubes_nb/nb_3950_2016-09-19T10:42:01_scans=0-4_corrected_im.fits'
tmpname = 'cubes_nb/tmp.fits'


val1 = red_fitsgetkeyword(fname, 'DW3', field_specifiers = fs1)


red_fitscube_getwcs, fname $
                     , coordinates = coordinates $
                     , distortions = distortions, iscan = 2

stop

file_copy, fname, tmpname, /over

hdr = headfits(tmpname)
indx = where(strmid(hdr, 0, 8) ne 'HIERARCH')
hdr2 = hdr(indx)
red_fitscube_newheader, tmpname, hdr2

val2 = red_fitsgetkeyword(tmpname, 'DW3', field_specifiers = fs2)

red_fitscube_getwcs, tmpname $
                     , coordinates = coordinates2 $
                     , distortions = distortions2

end
