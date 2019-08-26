; docformat = 'rst'

;+
; Correct a fitscube for intensity variation between observation time
; and prefilter calibration time. 
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
; :Params:
; 
;     filename : in, type=string
; 
;       The name of the file containing the fitscube.
; 
; 
; :Keywords:
; 
;     flip : in, optional, type=boolean
;   
;       Produce a flipped version if this keyword is set. 
; 
;     nostatistics : in, optional, type=boolean
;  
;       Do not calculate statistics metadata to put in header keywords
;       DATA*. If statistics keywords already exist, then remove them.
;
; :History:
; 
;   2019-08-26 : MGL. First version.
; 
;-
pro red::fitscube_intensitycorr, filename  $
                                 , flip = flip $
                                 , nostatistics = nostatistics

  inam = red_subprogram(/low, calling = inam1)

  if ~file_test(filename) then stop
  
  ;; Make prpara
  red_make_prpara, prpara, filename
  

  
  red_fitscube_open, filename, fileassoc, fitscube_info, /update

  hdr = fitscube_info.header
  
  ;; Information about processing steps in the formation of the input
  ;; file. 
  prprocs = fxpar(hdr, 'PRPROC*')
  prparas = fxpar(hdr, 'PRPARA*')

  ;; Cube dimensions
  Nx      = fitscube_info.dimensions[0]
  Ny      = fitscube_info.dimensions[1]
  Ntuning = fitscube_info.dimensions[2]
  Nstokes = fitscube_info.dimensions[3]
  Nscans  = fitscube_info.dimensions[4]

  Nframes = long(Ntuning) * long(Nstokes) * long(Nscans)
  
  ;; Check that it is not already intensity corrected.
  pos = where(strmatch(prprocs, inam), Nmatch)
  if Nmatch gt 0 then begin
    print
    print, inam + ' : This file is already intensity corrected:'
    print, filename
    red_fitscube_close, fileassoc, fitscube_info
    return
  endif
  
;  ;; Get scan numbers
;  scannum = red_fitsgetkeyword(filename, 'SCANNUM', variable_values = variable_values)
;  if n_elements(variable_values) gt 0 then begin
;    scannumbers = reform(variable_values.values)
;  endif else begin
;    scannumbers = [scannum]
;  endelse
;  
;  ;; Get medians of the I component of the first scan, to be used for
;  ;; selecting the wavelength points.
;  medi = fltarr(Ntuning)
;  for ituning = 0, Ntuning-1 do begin
;    red_fitscube_getframe, filename, frame, istokes = 0, iscan = 0, ituning = ituning
;    medi[ituning] = median(frame)
;  endfor                        ; ituning
;
;  if n_elements(tuning_selection) gt 0 then begin
;    ppc = tuning_selection
;    ;; Translate negative indices to positive ones
;    negindx = where(ppc lt 0, Nwhere)
;    if Nwhere gt 0 then begin
;      ppc[negindx] = Ntuning + ppc[negindx]
;    endif
;  endif else begin
;    ;; Choose spectral points to use. We want as little signal as
;    ;; possible so continuum points are good. For wide lines we
;    ;; might not have them so pick end points if similar intensity,
;    ;; or just one endpoint if one has significantly higher
;    ;; intensity than the other.
;    print, 'Select spectral points to calculate cross-talk from. Select with left mouse, end with right mouse.'
;    ppc = red_select_spoints(wav, medi)
;    tuning_selection = ppc
;  endelse
;
;  
;  if n_elements(ppc) gt 1 then begin
;    im = 0.
;    for i = 0, n_elements(ppc)-1 do begin
;      red_fitscube_getframe, filename, frame, istokes = 3, iscan = 0, ituning = ppc[i]
;      im += abs(frame)
;    endfor
;  endif else begin
;    red_fitscube_getframe, filename, im, istokes = 3, iscan = 0, ituning = ppc[0]
;  endelse
;
;  if n_elements(mag_mask) eq 0 then begin
;    print, 'Deselect areas with magnetic structures and/or artifacts. End with File-->Quit.'
;    ;;mag_mask = red_select_area(red_histo_opt(im,2.e-3), /noedge, /xroi)
;    mag_mask = red_select_area(red_histo_opt(im,2.e-3), /xroi)
;  end


  ;; Get WB prefilter
  pos = where(strmatch(prprocs, '*make_*_cube'), Nmatch)
  if Nmatch eq 0 then stop
  cube_paras = prparas[pos]
  wbpref = (stregex((json_parse(cube_paras, /tostruct)).dir $
                    , '/([0-9][0-9][0-9][0-9])/', /extract,/subex))[1]
  
  
;  ;; Get name of WB cube from the cube-making parameters.
;  pos = where(strmatch(prprocs, '*make_*_cube'), Nmatch)
;  if Nmatch eq 0 then stop
;
;  make_nb_cube_paras = prparas[pos]
;  wcfile = (json_parse(make_nb_cube_paras, /tostruct)).wcfile
;  whdr = headfits(wcfile)
;
;  ;; Get parameters from the wideband cube
;  if ~file_test(wcfile) then begin
;    print, inam + ' : Cannot find WB cube: '+wcfile
;    return
;  endif
;  
;  fxbopen, bunit, wcfile, 'MWCINFO', bbhdr
;  fxbreadm, bunit, row = 1 $
;              , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01'] $
;            ,   ANG, wcCROP, wcFF, wcGRID, wcND, wcSHIFT, wcTMEAN, wcX01Y01
;  ;; Note that the strarr wfiles cannot be read by fxbreadm! Put it in
;  ;; wbgfiles (WideBand Global).
;  fxbread, bunit, wbgfiles, 'WFILES', 1
;  fxbclose, bunit
;
;  
;  prefilter = fxpar(whdr, 'FILTER1') 
  wbfitfile = 'prefilter_fits/wb/wb_fit_'+wbpref+'.fits'
  if ~file_test(wbfitfile) then begin
    s = ''
    print, inam + ' : No WB fit file : '+wbfitfile
    print
    print, "It looks like you haven't run the a->fit_wb_diskcenter step."
    print, "If you have the calibration data, you can choose to delete the"
    print, "cube now and come back after running that calibration or to"
    print, "continue without intensity correction."
    print
    read, 'Delete [Y/n]', s
    if strmid(s, 0, 1) ne 'n' or strmid(s, 0, 1) ne 'N' then begin
      print
      print, inam+' : Deleting '+filename
      print
      file_delete, filename
      retall
    endif else begin
      print
      print, inam+' : No correction!'
      print
      return
    endelse
  endif
  
  wbpp = readfits(wbfitfile, pphdr)
  wbfitexpr = fxpar(pphdr, 'FITEXPR') ; Read the mpfitexpr fit function
  wb_time_beg = red_time2double(fxpar(pphdr, 'TIME-BEG'))
  wb_time_end = red_time2double(fxpar(pphdr, 'TIME-END'))
  case wbfitexpr of
    'P[0] + X*P[1]'            : begin ; Allow 10 min extrapolation
      wb_time_beg -= 10 * 60 
      wb_time_end += 10 * 60        
    end
    'P[0] + X*P[1] + X*X*P[2]' : begin ; Allow 30 min extrapolation
      wb_time_beg -= 30 * 60 
      wb_time_end += 30 * 60        
    end
    else :                      ; No extrapolation
  endcase
  


  ;; Load prefilter fit results to get the time of the prefilter fit
  ;; calibration data.
  prffiles = file_search('prefilter_fits/*_prefilter.idlsave', count = Nprf)
  time_avg_sum = 0D
  time_avg_n = 0L
  for iprf = 0L, Nprf-1 do begin

    nbpref = (strsplit(file_basename(prffiles[iprf]), '_', /extract))[1]
    if ~(self -> match_prefilters(nbpref, wbpref)) then continue
    
    restore, prffiles[iprf]     ; Restores variable prf which is a struct

    time_avg_sum += n_elements(prf.wav) * prf.time_avg
    time_avg_n   += n_elements(prf.wav)

    xposure = prf.xposure
    
  endfor                        ; iprf
  t_calib = time_avg_sum / time_avg_n

  
  
  ;; We need the WCS time coordinates 
  red_fitscube_getwcs, filename, coordinates = coordinates
  t = reform(coordinates.time[0, 0], Ntuning, Nstokes, Nscans)

  ;; Check that we are not extrapolating (too far).
  if ~file_test(wbfitfile) then begin
    s = ''
    print, inam + ' : No WB fit file : '+wbfitfile
    print
    print, inam + "It looks like your a->fit_wb_diskcenter step didn't find WB data"
    print, "for a long enough time range. Either your prefilter fit calibration data"
    print, "or your data cube have time stamps outside the range:"
    print, 'WB data range (plus margin) : ['+red_timestring(wb_time_beg)+','+red_timestring(wb_time_end)+'].'
    print, 'Cube time range : ['+red_timestring(min(t))+','+red_timestring(max(t))+'].'
    print
    print, "If you have more calibration data, you can choose to delete the"
    print, "cube now and come back after running that calibration again or to"
    print, "continue without intensity correction."
    print
    read, 'Delete [Y/n]', s
    if strmid(s, 0, 1) ne 'n' or strmid(s, 0, 1) ne 'N' then begin
      print
      print, inam+' : Deleting '+filename
      print
      file_delete, filename
      retall
    endif else begin
      print
      print, inam+' : No correction!'
      print
      return
    endelse
  endif


  
  ;; Interpolate to get the WB intensities
  wbints = red_evalexpr(wbfitexpr, t, wbpp) 
  wbintcalib = red_evalexpr(wbfitexpr, t_calib, wbpp) 


  ;; Change intensity to compensate for time difference
  ;; between prefilterfit and data collection. Use the ratio
  ;; of (WB intensity)/(exposure time) for interpolated times
  ;; of calibration data and cube data.
  wbratio = wbintcalib / wbints

  ;; Need also ratio of exposure times to compensate the NB
  ;; data for different exposure time in calibration data and
  ;; in cube data! MOMFBD processing divides output with the
  ;; number of frames, so the single-frame exposure time is
  ;; the appropriate one.
  xpratio = xposure / fxpar(hdr,'TEXPOSUR')
  
  iframe = 0L
  for iscan = 0, Nscans-1 do begin
    for istokes=0, Nstokes-1 do begin
      for ituning = 0, Ntuning-1 do begin

        red_progressbar, iframe, Nframes, 'Correcting intensity of NB data'
        
        ;; Read a frame from the fitscube file
        red_fitscube_getframe, fileassoc, frame $
                               , iscan = iscan, istokes = istokes, ituning = ituning
        
        ;; Write the corrected frame back to the fitscube file
        red_fitscube_addframe, fileassoc, frame $
                               * wbratio[ituning, istokes, iscan] $
                               * xpratio $
                               , iscan = iscan, istokes = istokes, ituning = ituning

        iframe++
        
      endfor                    ; ituning
    endfor                      ; istokes
  endfor                        ; iscan


  ;; Possibly do it also for WB data of scan cubes!


  
  
;  for iscan = 0, Nscans-1 do begin
;
;    if makemask then begin
;      ;; Construct a mask for the padding
;      pad_mask = make_array(Nxx, Nyy, /float, value = 1.) 
;      pad_mask = red_rotation(pad_mask, ang[iscan], wcshift[0,iscan], wcshift[1,iscan], background = 0, full = wcFF)
;      pindx = where(pad_mask le 0.99) ; Pixels that are padding
;    
;      ;; Include the padding mask just in case it rotates into the
;      ;; selected mask.
;      this_mask = mag_mask * pad_mask
;      mindx = where(this_mask)
;    endif else begin
;      mindx = where(mag_mask)
;;      mindx = lindgen(Nx, Ny)
;    endelse 
;      
;    ;;crt = red_get_ctalk(d, idx=ppc, mask=pixmask)
;    crt = dblarr(Nstokes)
;    numerator   = dblarr(Nstokes)
;    denominator = 0d
;    for i = 0, n_elements(ppc)-1 do begin
;      red_fitscube_getframe, filename, im0, istokes = 0, iscan = iscan, ituning = ppc[i] ; Stokes I
;
;      ;;denominator += median(im0[where(this_mask)] *
;      ;;im0[where(this_mask)], /double)
;      
;      ;; Find the centroid by fitting a Gaussian
;      a = red_histo_gaussfit(im0[mindx] * im0[mindx], FWlevel = 0.25)
;      denominator += a[1]
;      
;      for istokes=1, Nstokes-1 do begin
;        red_fitscube_getframe, filename, im, istokes = istokes, iscan = iscan, ituning = ppc[i]
;        ;;numerator[istokes] += median(im0[where(this_mask)] * im[where(this_mask)], /double)
;        a = red_histo_gaussfit(im0[mindx] * im[mindx], FWlevel = 0.25)
;        numerator[istokes] += a[1]
;      endfor                    ; istokes
;    endfor
;    crt = numerator/denominator 
;    
;    print, 'Scan '+strtrim(scannumbers[iscan], 2)+' : crosstalk from I -> Q,U,V =' $
;           , crt[1], ', ', crt[2], ', ', crt[3], format='(A,F8.5,A,F8.5,A,F8.5)'
;    
;    for ituning = 0, Ntuning-1 do begin
;      red_fitscube_getframe, filename, im0, istokes = 0, iscan = 0, ituning = ituning ; Stokes I
;      if makemask then im0[pindx] = median(im0[pindx]) ; Set the padding to median
;      red_fitscube_addframe, filename, im0, istokes = 0, iscan = 0, ituning = ituning ; Write with updated padding
;      for istokes=1, Nstokes-1 do begin
;        ;;d[*,*,tt,ww] -= crt[tt]*d[*,*,0,ww]
;        red_fitscube_getframe, filename, im, istokes = istokes, iscan = iscan, ituning = ituning
;        im -= float(crt[istokes] * im0)
;        if makemask then im[pindx] = median(im[pindx]) 
;        red_fitscube_addframe, filename, im, istokes = istokes, iscan = iscan, ituning = ituning
;      endfor                    ; istokes
;    endfor                      ; ituning
;  endfor                        ; iscan

  ;; Add info about this step
  self -> headerinfo_addstep, hdr $
                              , prstep = 'Intensity correction' $
                              , prpara = prpara $
                              , prproc = inam

  ;; Close the file and write the updated header
  red_fitscube_close, fileassoc, fitscube_info, newheader = hdr

;  if keyword_set(nostatistics) then begin
;
;    ;; Remove any existing statistics keywords from the file. 
;
;    ;; Search for DATA* keywords
;    dindx = where(strmid(hdr, 0, 4) eq 'DATA', Ndata)
;    ;; Loop through the keywords backwards so we don't move
;    ;; them before they are deleted.
;    for idata = Ndata-1, 0, -1 do begin
;      keyword = strtrim(strmid(hdr[dindx[idata]], 0, 8), 2)
;      red_fitsdelkeyword, hdr, keyword
;    endfor                      ; idata
;    
;  endif else begin
;
;    red_fitscube_statistics, filename, /write $
;                             , angles = angles $
;                             , cube_comments = cube_comments $
;                             , full = wcFF $
;                             , grid = wcGRID $
;                             , origNx = Nxx $
;                             , origNy = Nyy $
;                             , shifts = wcSHIFTS 
;
;  endelse
;  
;  if keyword_set(flip) then begin
;    self -> fitscube_flip, filename, flipfile = flipfile
;  endif

  
  ;; For scan cubes, do it also for the WB image.
  fits_info, filename, /SILENT , N_ext = n_ext, EXTNAME=extnames
  if n_ext gt 0 && round(total(strmatch(strtrim(extnames,2),'WBIMAGE'))) eq 1 then begin
    ;; After this, WB data will not turn SI units (like the NB data).
    ;; But it should unify the WB intensity scaling during the day,
    ;; providing a way to check that the scaling is correct by
    ;; comparing WB data from different times.
    wbimage = mrdfits(filename, 'WBIMAGE', ehdr, STATUS=status, /silent)
    wbimage *= mean(wbratio) * xpratio
    modfits, filename, float(wbimage), ehdr, errmsg = errmsg, extname = 'WBIMAGE'
  endif
  
end
