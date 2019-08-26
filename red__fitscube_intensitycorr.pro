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
; :History:
; 
;   2019-08-26 : MGL. First version.
; 
;-
pro red::fitscube_intensitycorr, filename  

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

  ;; Get WB prefilter
  pos = where(strmatch(prprocs, '*make_*_cube'), Nmatch)
  if Nmatch eq 0 then stop
  cube_paras = prparas[pos]
  wbpref = (stregex((json_parse(cube_paras, /tostruct)).dir $
                    , '/([0-9][0-9][0-9][0-9])/', /extract,/subex))[1]

  wbfitfile = 'wb_intensities/wb_fit_'+wbpref+'.fits'
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
  t = reform(coordinates.time[0, 0], Ntuning, Nscans)

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
        red_fitscube_addframe, fileassoc $
                               , frame * wbratio[ituning, iscan] * xpratio $
                               , iscan = iscan, istokes = istokes, ituning = ituning

        iframe++
        
      endfor                    ; ituning
    endfor                      ; istokes
  endfor                        ; iscan


  ;; Possibly do it also for WB data of scan cubes!

  ;; Add info about this step
  self -> headerinfo_addstep, hdr $
                              , prstep = 'INTENSITY-CALIBRATION' $
                              , prpara = prpara $
                              , prproc = inam

  ;; Close the file and write the updated header
  red_fitscube_close, fileassoc, fitscube_info, newheader = hdr

  ;; Remove any old statistics info. 
  red_fitscube_statistics, filename, /remove_only
  
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
