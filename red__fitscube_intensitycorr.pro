; docformat = 'rst'

;+
; Correct a fitscube for intensity variation between observation time
; and prefilter calibration time. 
; 
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
; :Keywords:
; 
;     intensitycorrmethod : in, optional, type="string or boolean", default='fit'
;
;       One of 'old' (for correction based on comparing the current
;       scan with the prefilter calibration data), or 'fit' (for
;       correction based on comparing intensity from
;       fit_wb_diskcenter, evaluated at the current time and at the
;       time of the prefilter calibration). Boolean value TRUE is
;       equivalent to 'fit'. FALSE and 'none' (really any string value
;       that is not 'fit' or 'old') results in no correction.
;
; :History:
; 
;   2019-08-26 : MGL. First version.
; 
;   2019-10-10 : MGL. New keyword nodiskcenter.
; 
;   2020-01-22 : MGL. Remove keyword nodiskcenter, new keyword
;                corrmethod. 
; 
;   2020-10-28 : MGL. Remove statistics calculations.
;
;-
pro red::fitscube_intensitycorr, filename $
                                 , intensitycorrmethod = corrmethod
  

  inam = red_subprogram(/low, calling = inam1)

  if ~file_test(filename) then stop

  case 1 of

    ;; Default
    n_elements(corrmethod) eq 0 : corrmethod = 'fit' 
    
    ;; Undefined strings
    size(corrmethod, /tname) eq 'STRING' : begin 
      if corrmethod ne 'fit' and corrmethod ne 'old' then corrmethod = 'none'
    end

    ;; Other types interpreted as boolean
    else : if corrmethod then corrmethod = 'fit' else corrmethod = 'none'
    
  endcase

  if corrmethod eq 'none' then begin
    ;; No correction, just return without doing anything.
    print, inam + ': No correction'
    return
  endif 

  ;; Modes to store in the processing step metadata
  if corrmethod eq 'old' then begin
    PRMODE = 'LOCAL'
  endif else begin
    PRMODE = 'DISK-CENTER'
  endelse

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
  
  ;; Check that it is not already intensity corrected. (For cubes with
  ;; RESPAPPL, we might consider undoing the old corrections and
  ;; applying new ones.)
  pos = where(strmatch(prprocs, inam), Nmatch)
  if Nmatch gt 0 then begin
    print
    print, inam + ' : This file is already intensity corrected:'
    print, filename
    red_fitscube_close, fileassoc, fitscube_info
    return
  endif
  
  ;; Get info from the cube making step
  pos_makenb   = where(strmatch(prprocs, '*make_nb_cube'  ), Nmakenb  )
  pos_makescan = where(strmatch(prprocs, '*make_scan_cube'), Nmakescan)
  case 1 of
    Nmakenb gt 0 : begin        ; This is a NB cube
      cube_paras = prparas[pos_makenb[0]]
      cube_paras_struct = json_parse(cube_paras, /tostruct)
      wbhdr = headfits(cube_paras_struct.wcfile)
      wbpref = strtrim(fxpar(wbhdr, 'FILTER1'), 2)
      fxbopen, bunit, cube_paras_struct.wcfile, 'MWCINFO', bbhdr
      fxbreadm, bunit, row = 1 $
                , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01', 'DIRECTION'] $
                ,   ANG, wcCROP, wcFF, wcGRID, wcND, wcSHIFT, wcTMEAN, wcX01Y01,   direction
      fxbclose, bunit
    end
    Nmakescan gt 0 : begin      ; This is a SCAN cube
      cube_paras = prparas[pos_makescan[0]]
      cube_paras_struct = json_parse(cube_paras, /tostruct)
      wbpref = (stregex(cube_paras_struct.dir $
                        , '/([0-9][0-9][0-9][0-9])/', /extract,/subex))[1]
      wbimage = mrdfits(filename, 'WBIMAGE', ehdr, STATUS=status, /silent)
      wcTMEAN = median(wbimage)
    end
    else : begin
      print, inam+' : This type of cube is not implemented yet.'
      stop
    end
  endcase

  print, wcTMEAN
;  stop
  

  pfiles = file_search(self.out_dir $
                       + '/prefilter_fits/*_prefilter.idlsave' $
                       , count = Npfiles)
  prefilters = strarr(Npfiles)
  for ipref = 0, Npfiles-1 do begin
    splt = strsplit(file_basename(pfiles[ipref]), '_', /extract)
    prefilters[ipref] = splt[1]
  endfor
  pindx = where(self -> match_prefilters(prefilters, strtrim(fxpar(hdr,'FILTER1'), 2)), Npref)
  if Npref eq 0 then stop
  prefilters = prefilters[pindx]
  ppfiles = pfiles[pindx]
  
;  ;; Prefilter fit data - WB intensity and time. 
;  case strtrim(fxpar(hdr,'INSTRUME'), 2) of
;
;    ;; For CRISP, R and T are simultaneous with the same WB data.
;    'CRISP'   : pfiles = file_search(self.out_dir $
;                                     + '/prefilter_fits/Crisp-?_' $
;                                     + strtrim(fxpar(hdr,'FILTER1'), 2) $
;                                     + '_prefilter.idlsave' $
;                                     , count = Npfiles)
;
;    ;; For CHROMIS, scans are fast so one prefilter should be enough
;    ;; also for Ca II.
;    'CHROMIS' : pfiles = file_search(self.out_dir $
;                                     + '/prefilter_fits/chromis_' $
;                                     + strtrim(fxpar(hdr,'FILTER1'), 2) $
;                                     + '_prefilter.idlsave' $
;                                     , count = Npfiles)
;    else : stop
;  endcase 
;  restore, pfiles[0]            ; Restores variable prf which is a struct
;  time_avg_sum = n_elements(prf.wav) * prf.time_avg
;  time_avg_n   = n_elements(prf.wav)
;  t_calib = time_avg_sum / time_avg_n
;  prefilter_wb = prf.wbint
;  xposure = prf.xposure
  
;  ;; Load prefilter fit results to get the time of the prefilter fit
;  ;; calibration data.
;  prffiles = file_search('prefilter_fits/*_prefilter.idlsave', count = Nprf)
  time_avg_sum = 0D
  time_avg_n = 0L
  for ipref = 0L, n_elements(ppfiles)-1 do begin

;    nbpref = (strsplit(file_basename(prffiles[iprf]), '_', /extract))[1]
;    if ~(self -> match_prefilters(nbpref, wbpref)) then continue
    
    restore, ppfiles[ipref]     ; Restores variable prf which is a struct

    time_avg_sum += n_elements(prf.wav) * prf.time_avg
    time_avg_n   += n_elements(prf.wav)
    
  endfor                        ; iprf
  t_calib = time_avg_sum / time_avg_n
  xposure = prf.xposure
  prefilter_wb = prf.wbint      ; The mean wb intensity from the fitprefilter step.

;  stop
  

  if PRMODE eq 'DISK-CENTER' then begin

    ;; Do DISK-CENTER based correction
    
    wbfitfile = 'wb_intensities/wb_fit_'+wbpref+'.fits'
    if ~file_test(wbfitfile) then begin
      s = ''
      print, inam + ' : No WB fit file : '+wbfitfile
      print
      print, "It looks like you haven't run the a->fit_wb_diskcenter step."
      print
      print, "If you have the calibration data, you can choose to delete the"
      print, "cube now and come back after running that calibration."
      print, ""
      print, "If you do not have the calibration data, do not delete."
      print, "Just continue without intensity correction and then (if you wish)"
      print, "do the old correction afterwards with"
      print, "IDL> a -> fitscube_intensitycorr, filename, intensitycorrmethod = 'old'"
      print
      print, "When making NB cubes from this day in the future, you can also add the "
      print, "intensitycorrmethod = 'old' keyword to make_nb_cube and then not do"
      print, "the fitscube_intensitycorr step afterwards."
      print
      read, 'Delete [Y/n]', s
      if strlowcase(strmid(s, 0, 1)) ne 'n' then begin
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
      else :                    ; No extrapolation
    endcase
    

    ;; We need the WCS time coordinates 
    red_fitscube_getwcs, filename, coordinates = coordinates
    t = reform(coordinates.time[0, 0], Ntuning, Nscans)

    ;; Check that we are not extrapolating (too far).
    if min(t) lt wb_time_beg or max(t) gt wb_time_end then begin
      s = ''
      print, inam + ' : No WB fit file : '+wbfitfile
      print
      print, inam + "It looks like your a->fit_wb_diskcenter step didn't find WB data"
      print, "for a long enough time range. Either your prefilter fit calibration data"
      print, "or your data cube have time stamps outside the range:"
      print, 'WB data range (plus margin) : [' + red_timestring(wb_time_beg) $
             + ',' + red_timestring(wb_time_end) + '].'
      print, 'Cube time range : [' + red_timestring(min(t)) $
             + ',' + red_timestring(max(t)) + '].'
      print
      print, "If you have more calibration data, you can choose to delete the"
      print, "cube now and come back after running that calibration again or to"
      print, "continue without intensity correction."
      print
      read, 'Delete [Y/n]', s
      if strlowcase(strmid(s, 0, 1)) ne 'n' then begin
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


    
    ;; Use the fit to get the WB intensities
    wbints = red_evalexpr(wbfitexpr, t/3600, wbpp) 
    wbintcalib = red_evalexpr(wbfitexpr, t_calib/3600, wbpp)

    ;; Set prref to string representation of wbfitexpr with wbpp
    ;; coefficients.
    prref = 'Median DC WB intensity fit in counts as fcn of x=t/1h : ' + red_renderexpr(wbfitexpr,wbpp) 
            
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

;    oldcorrection = 1d / wcTMEAN
;    fitcorrection = wbratio * xpratio
;    stop
    
    ;; We want the scan-to-scan variations of the "LOCAL" old kind of
    ;; correction. But scaled to the WB intensity calibration.
    correction = mean(wbratio * xpratio) * mean(wcTMEAN) / wcTMEAN

    if 0 then begin
      alt_wbratio = prefilter_wb/wcTMEAN
      alt_correction = replicate(wbratio, nscans)

      print, correction
      print, alt_correction
      print, wbratio
      print, alt_wbratio
      stop
    endif
    
  endif else begin

    ;; The "old" method. Use ratio of *average* WB intensities.

    ;; Calculate wb ratio or get it from calling program? We need the
    ;; median prefilter fit calibration WB intensity over the median
    ;; intensities of the WB cube. However, while wcTMEAN is a
    ;; median() value, prefilter_wb is a mean(). Shouldn't matter very
    ;; much as prefilter_wb is from granulation data so mean() and
    ;; median() should be close.
    case 1 of
      Nmakenb gt 0 : begin      ; This is a NB cube
        wbratio = mean(prefilter_wb/wcTMEAN)
        correction = prefilter_wb/wcTMEAN
      end
      Nmakescan gt 0 : begin    ; This is a SCAN cube
        wbratio = prefilter_wb/wcTMEAN
        correction = prefilter_wb/wcTMEAN
      end
    endcase
    
    ;; Correction for wcTMEAN is normalized to unit mean, so corrects
    ;; for irregularities in the WB intensities. But we need the
    ;; wbratio to get the correction for elevation over the whole day. 
                                ;   correction = wbratio / wcTMEAN
                                ;   stop
;    correction=prefilter_wb/wcTMEAN
;    correction = replicate(wbratio, nscans)
    ;; Set prref to file name of reference WB image from the
    ;; calibration data set.
    prref = 'Mean DC WB median intensity in counts from prefilter fit data : '+strtrim(prefilter_wb, 2) 
    
  endelse 

  ;; Apply the corrections
  print, inam + ' : Corrections = ', correction

;  ;; Does correction have the correct dimensions?
;  help, correction
;  stop
  
  iframe = 0L
  for iscan = 0, Nscans-1 do begin
    for istokes=0, Nstokes-1 do begin
      for ituning = 0, Ntuning-1 do begin

        red_progressbar, iframe, Nframes, /predict $
                         , 'Correcting intensity of NB data'
        
        ;; Read a frame from the fitscube file
        red_fitscube_getframe, fileassoc, frame $
                               , iscan = iscan, istokes = istokes, ituning = ituning
        
        ;; Write the corrected frame back to the fitscube file
        red_fitscube_addframe, fileassoc, frame * correction[iscan] $
                               , iscan = iscan, istokes = istokes, ituning = ituning

        iframe++
        
      endfor                    ; ituning
    endfor                      ; istokes
  endfor                        ; iscan

  ;; Add info about this step
  self -> headerinfo_addstep, hdr $
                              , prstep = 'CALIBRATION-INTENSITY-TEMPORAL' $
                              , prpara = prpara $
                              , prref = prref $
                              , prproc = inam $
                              , prmode = prmode
 
  ;; Close the file and write the updated header
  red_fitscube_close, fileassoc, fitscube_info, newheader = hdr

  ;; The applied correction should be saved as a variable keyword
  ;; RESPAPPL (APPLied RESPonse function).
  red_fitscube_addrespappl, filename, correction, /scans, /update
  
  ;; For scan cubes, do it also for the WB image.
  fits_info, filename, /SILENT , N_ext = n_ext, EXTNAME=extnames
  if n_ext gt 0 && round(total(strmatch(strtrim(extnames,2),'WBIMAGE'))) eq 1 then begin
    ;; After this, WB intensities will not be in SI units (like the NB
    ;; data). But it should (for DISK-CENTER) unify the WB intensity
    ;; scaling during the day, providing a way to check that the
    ;; scaling is correct by comparing WB data from different times.
    wbimage = mrdfits(filename, 'WBIMAGE', ehdr, STATUS=status, /silent)
    wbimage *= mean(correction)
    modfits, filename, float(wbimage), ehdr, errmsg = errmsg, extname = 'WBIMAGE'
  endif

end


a = crispred(/dev)

a -> fitscube_intensitycorr, 'cubes_nb/nb_6302_2016-09-19T09:30:20_scans=12-16_stokes_corrected_im.fits'


end

; Todo:
;
;  1. Integrate wcTMEAN intensity correction based on the current
;  dataset (see make_nb_cube). Mark in the step info "PRMODEn" keyword
;  which kind of WB intensity correction was done: DISK-CENTER or
;  LOCAL, where the former is the new method and LOCAL is the old
;  method that removes center-to-limb variations. Make DISK-CENTER
;  correction the default and LOCAL available with an optional
;  IDL keyword. 
; 
;  2. Remove the wcTMEAN correction from make_nb_cube and instead call
;  fitecube_intensitycorr. 
; 
;  3. Testing: For a day with data data sets at similar mu angles,
;  separated by hours, run make_scan_cube for scans separated in time
;  with both DISK-CENTER and LOCAL. Verify that both corrections give
;  intensities that are consistent over time but that the DISK-CENTER
;  cubes have lower intensities - because the limb darkening is not
;  removed. For DC data, both methods should give similar intensities.
; 
;  4. If possible, test also chromis data with different exposure
;  times.
; 
