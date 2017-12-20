; docformat = 'rst'

;+
; Convert an old LP format science data cube to a fitscube. 
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
; :Params:
; 
;    inname : in, type=string
; 
;       The path to the file to be converted.
; 
; 
; :Keywords:
; 
;    headerdata : in, optional, type=strarr
;
;      FITS header with keywords to be added to the output file.
;
;    headerfile : in, optional, type=string
;
;      The name of a file where headerdata can be found.
;
;    outname : in, optional, type=string, default = inname+'.fits'
;
;      Where to write the output. A spectral cube might also be
;      written, the file name for this will be generated based on
;      outname. 
; 
;    polarimetric : in, optional, type=boolean
;
;      Whether this is a polarimetric cube.
;
; :History:
; 
;    2017-12-06 : MGL. First version.
; 
;-
pro red::fitscube_convertlp, inname $
                             , headerdata = headerdata $
                             , headerfile = headerfile $
                             , outdir = outdir $
                             , outname = outname $
                             , overwrite = overwrite $
                             , polarimetric = polarimetric

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  ;; Make prpara
  red_make_prpara, prpara, headerdata         
  red_make_prpara, prpara, headerfile       
  red_make_prpara, prpara, outdir      
  red_make_prpara, prpara, outname      
  red_make_prpara, prpara, overwrite       
  red_make_prpara, prpara, polarimetric 
  
  dir = file_dirname(inname)+'/'
  iname = file_basename(inname)


  if n_elements(outdir) eq 0 then outdir = dir
  
  if n_elements(outname) eq 0 then begin
    
    iname = file_basename(inname)

    iname = red_strreplace(iname, '.fcube', '')
    iname = red_strreplace(iname, '.icube', '')

    oname = iname + '_fromlp_im.fits'
    oname = outdir + oname
    
  endif

  ;; Already done?
  if file_test(oname) then begin
    if keyword_set(overwrite) then begin
      print, 'Overwriting existing data cube:'
      print, oname
    endif else begin
      print, 'This data cube exists already:'
      print, oname
      return
    endelse
  endif

  if keyword_set(polarimetric) then Nstokes = 4 else Nstokes = 1

  ;; Cameras and detectors
  self->getdetectors
  wbindx = where(strmatch(*self.cameras,'*-W')) ; Matches CRISP and CHROMIS WB cameras
  wbcamera = (*self.cameras)[wbindx[0]]
  wbdetector = (*self.detectors)[wbindx[0]]
  nbindx = where(strmatch(*self.cameras,'*-[NRT]')) ; Matches CRISP and CHROMIS NB cameras
  nbcamera = (*self.cameras)[nbindx[0]]
  nbdetector = (*self.detectors)[nbindx[0]]


  
  file_mkdir, file_dirname(oname)

  red_lp_header, inname, header=header, datatype=datatype, $
                 dims=dims, nx=nx, ny=ny, nt=nt, endian=endian_file

  iname_parts = strsplit(iname, '_', /extract)
  pref = iname_parts[1]
  date_obs = iname_parts[2]
  timestamp = (strsplit(date_obs, 'T', /extract))[1]
  scans = (strsplit(iname_parts[3], '=', /extract))[1]

  
  ;; Scan numbers
  s_array = red_expandrange(scans)
  Nscans = n_elements(s_array)

  nbdir = (*self.data_dirs)[(where(strmatch(*self.data_dirs,'*'+timestamp), Nmatch))[0]] $
          + '/' + nbcamera + '/'
  wbdir = (*self.data_dirs)[(where(strmatch(*self.data_dirs,'*'+timestamp), Nmatch))[0]] $
          + '/' + wbcamera + '/' 
  
  ;; Get tuning info from the spectfile
  spectfile = dir+'spectfile.'+pref+'.idlsave'
  if file_test(spectfile) then begin
    restore, spectfile
    ;; This reads three quantities: SPECT_POS, NORM_FACTOR, NORM_SPECT
    Nwav = n_elements(SPECT_POS)
    is_wb = 0
    rawdir = nbdir
  endif else begin
    Nwav = 1
    is_wb = 1
    rawdir = wbdir
  endelse
 


  ;; Dimensions check!
  if Nt ne Nstokes * Nscans * Nwav then stop
  dims = [Nx, Ny, Nwav, Nstokes, Nscans]

  wbfiles = file_search(wbdir+'*', count = Nwbfiles)
  rawfiles = file_search(rawdir+'*', count = Nrawfiles)
  self->selectfiles, files = rawfiles, states = rawstates $
                     , scan = s_array $
                     , count = Nselected $
                     , selected = selected 
  if Nt ne Nselected then stop

  uindx = uniq(rawstates.tun_wavelength, sort(rawstates.tun_wavelength))
  ulambda = rawstates[uindx].tun_wavelength ; Should now be ordered as spect_pos!
  range = [min([ulambda*1e9, spect_pos/10d]), max([ulambda*1e9, spect_pos/10d])] + [-1, 1]*0.2
  cgplot, ulambda*1e9, spect_pos/10d, psym=-16, color = 'red', /yno $
          , xrange = range,  yrange = range $
          , xtitle = 'From raw data: ulambda / 1 nm', ytitle = 'From spectfile: spect_pos / 1 nm'

  ustates = rawstates[uindx].fullstate

  
  ;; Get times from NB file headers.
  t_array = dblarr(Nscans)               ; WB time
  tbeg_array     = dblarr(Nwav, Nscans)  ; Time beginning for state
  tavg_array     = dblarr(Nwav, Nscans)  ; Time average for state
  tend_array     = dblarr(Nwav, Nscans)  ; Time end for state
  date_beg_array = strarr(Nwav, Nscans)  ; DATE-BEG for state
  date_avg_array = strarr(Nwav, Nscans)  ; DATE-AVG for state
  date_end_array = strarr(Nwav, Nscans)  ; DATE-END for state
  exp_array      = fltarr(Nwav, Nscans)  ; Total exposure time
  sexp_array     = fltarr(Nwav, Nscans)  ; Single exposure time
  nsum_array     = lonarr(Nwav, Nscans)  ; Number of summed exposures
  
  for iscan = 0L, Nscans-1 do begin

    red_progressbar, iscan, Nscans, /predict, 'Reading file headers'

    ;; WB times for log file access
    self->selectfiles, files = wbfiles, states = wbstates $
                       , scan = s_array[iscan] $
                       , count = Nselected $
                       , selected = selected
    these_tavg_array = dblarr(Nselected)
    for iselected = 0L, Nselected-1 do begin
      thishdr = red_readhead(wbfiles[selected[iselected]])
      red_fitspar_getdates, thishdr, date_avg = date_avg 
      these_tavg_array[iselected] = red_time2double((strsplit(date_avg, 'T', /extract))[1])
    endfor                      ; iselected
    t_array[iscan] = mean(these_tavg_array)

    ;; NB times for WCS
    for iwav = 0, Nwav-1 do begin
      self->selectfiles, files = rawfiles, states = rawstates $
                         , scan = s_array[iscan] $
                         , ustat = ustates[iwav] $
                         , count = Nselected $
                         , selected = selected
      these_tbeg_array = dblarr(Nselected)
      these_tavg_array = dblarr(Nselected)
      these_tend_array = dblarr(Nselected)
      these_exp_array = dblarr(Nselected)
      these_sexp_array = dblarr(Nselected)
      these_nsum_array = dblarr(Nselected)
      for iselected = 0L, Nselected-1 do begin
        thishdr = red_readhead(rawfiles[selected[iselected]])
        red_fitspar_getdates, thishdr $
                              , date_beg = date_beg $
                              , date_end = date_end $
                              , date_avg = date_avg 
        these_tbeg_array[iselected] = red_time2double((strsplit(date_beg, 'T', /extract))[1])
        these_tavg_array[iselected] = red_time2double((strsplit(date_avg, 'T', /extract))[1])
        these_tend_array[iselected] = red_time2double((strsplit(date_end, 'T', /extract))[1])
        these_sexp_array[iselected] = fxpar(thishdr, 'XPOSURE')
        these_nsum_array[iselected] = fxpar(thishdr, 'NAXIS3')
        these_exp_array[iselected] = these_sexp_array[iselected]*these_nsum_array[iselected] 
      endfor                    ; iselected
      date_beg_array[iwav, iscan] = date_beg
      date_end_array[iwav, iscan] = date_end
      date_avg_array[iwav, iscan] = date_avg
      tbeg_array[iwav,iscan] = min(these_tbeg_array)
      tavg_array[iwav,iscan] = mean(these_tavg_array)
      tend_array[iwav,iscan] = max(these_tend_array)
      exp_array[iwav,iscan]  = total(these_exp_array)
      sexp_array[iwav,iscan] = mean(these_sexp_array)
      nsum_array[iwav,iscan] = total(these_nsum_array)
    endfor                      ; iwav
    
  endfor                        ; iscan
  
  ;; Make header
  red_mkhdr, hdr, datatype, dims

  ;; Copy some keywords from the last raw data frame.
  anchor = 'DATE'
  red_fitscopykeyword, anchor = anchor, hdr, 'EXTNAME' , thishdr
  red_fitscopykeyword, anchor = anchor, hdr, 'SOLARNET', thishdr
  red_fitscopykeyword, anchor = anchor, hdr, 'OBS_HDU' , thishdr
  red_fitscopykeyword, anchor = anchor, hdr, 'ORIGIN'  , thishdr
  red_fitscopykeyword, anchor = anchor, hdr, 'TELESCOP', thishdr
  red_fitscopykeyword, anchor = anchor, hdr, 'INSTRUME', thishdr
  red_fitscopykeyword, anchor = anchor, hdr, 'CAMERA'  , thishdr
  red_fitscopykeyword, anchor = anchor, hdr, 'DETECTOR', thishdr
  red_fitscopykeyword, anchor = anchor, hdr, 'DATE-OBS', thishdr
  red_fitscopykeyword, anchor = anchor, hdr, 'WAVEUNIT', thishdr
  red_fitscopykeyword, anchor = anchor, hdr, 'OBSERVER', thishdr
  red_fitscopykeyword, anchor = anchor, hdr, 'OBJECT'  , thishdr
  red_fitscopykeyword, anchor = anchor, hdr, 'CADENCE' , thishdr
  
;  red_fitsaddkeyword, anchor = anchor, hdr, 'WAVELNTH', ?, 'Characteristic wavelength'
;  red_fitsaddkeyword, anchor = anchor, hdr, ''
;  red_fitsaddkeyword, anchor = anchor, hdr, ''
;  red_fitsaddkeyword, anchor = anchor, hdr, ''
;  red_fitsaddkeyword, anchor = anchor, hdr, ''
;  red_fitsaddkeyword, anchor = anchor, hdr, ''
  
;  red_fitsaddkeyword, anchor = anchor, hdr, 'FILENAME', oname
;  if is_wb then begin
;    red_fitsaddkeyword, anchor = anchor, hdr, 'CAMERA', wbcamera
;    red_fitsaddkeyword, anchor = anchor, hdr, 'DETECTOR', wbdetector
;  endif else begin
;    red_fitsaddkeyword, anchor = anchor, hdr, 'CAMERA', nbcamera
;    red_fitsaddkeyword, anchor = anchor, hdr, 'DETECTOR', nbdetector
;  endelse
  dateref = self.isodate+'T00:00:00.000000' ; Midnight
  red_fitsaddkeyword, anchor = anchor, hdr, 'DATEREF', dateref, 'Reference time in ISO-8601'

  red_fitsaddkeyword, anchor = anchor, hdr, 'BUNIT', 'dn', 'Units in array: digital number'
  red_fitsaddkeyword, anchor = anchor, hdr, 'BTYPE', 'Intensity', 'Type of data in array'

  ;; Add "global" metadata
  red_metadata_restore, hdr

  ;; Add info about this step
  self -> headerinfo_addstep, hdr $
                              , prstep = 'Convert science data cube from LP format' $
                              , prpara = prpara $
                              , prproc = inam

  self -> fitscube_initialize, oname, hdr, olun, fileassoc, dims 

  wcs = replicate({  wave:dblarr(2,2) $
                   , hplt:dblarr(2,2) $
                   , hpln:dblarr(2,2) $
                   , time:dblarr(2,2) $
                  }, Nwav, Nscans)

  ;; Get metadata from logfiles
  red_logdata, self.isodate, time_r0, r0 = metadata_r0
  red_logdata, self.isodate, time_pig, pig = metadata_pig, rsun = rsun

  ;; Get pointing at center of FOV
  red_wcs_hpl_coords, t_array, metadata_pig, time_pig $
                      , hpln, hplt
  
  ;; The alignment routine (red_aligncube) subtracts the median of the
  ;; cross-correlation measured image shifts but removes no trends. So
  ;; it really tries to make the pointing the same during the whole
  ;; sequence without allowing for drifts. So we should make the
  ;; pointing metadata constant in time, let's use the median:
  hpln = median(hpln)
  hplt = median(hplt)

  ;; But what we want to tabulate is the pointing in the corners of
  ;; the FOV. Assume hpln and hplt are the coordinates of the center
  ;; of the FOV.
  wcs.hpln[0, 0, *, *] = hpln - double(self.image_scale) * (Nx-1)/2.d
  wcs.hpln[1, 0, *, *] = hpln + double(self.image_scale) * (Nx-1)/2.d
  wcs.hpln[0, 1, *, *] = hpln - double(self.image_scale) * (Nx-1)/2.d
  wcs.hpln[1, 1, *, *] = hpln + double(self.image_scale) * (Nx-1)/2.d
  
  wcs.hplt[0, 0, *, *] = hplt - double(self.image_scale) * (Ny-1)/2.d
  wcs.hplt[1, 0, *, *] = hplt - double(self.image_scale) * (Ny-1)/2.d
  wcs.hplt[0, 1, *, *] = hplt + double(self.image_scale) * (Ny-1)/2.d
  wcs.hplt[1, 1, *, *] = hplt + double(self.image_scale) * (Ny-1)/2.d

  for iscan = 0L, Nscans-1 do begin
    for iwav = 0, Nwav-1 do begin
      ;; We rely here on hpln and hplt being the first two tabulated
      ;; coordinates. To make this more general, we should get the
      ;; actual indices from the headers. Maybe later...
      wcs[iwav, iscan].wave = spect_pos[iwav]/10d ; Å --> nm
      wcs[iwav, iscan].time = tavg_array[iwav, iscan]
    endfor                      ; iwav
  endfor                        ; iscan

  ;; Copy data frames
  if ((byte(1L, 0, 1))[0] eq 1) then endian = 'l' else endian='b'
  swap_endian = (datatype gt 1) and (endian ne endian_file)
  openr, ilun, inname, /get_lun, swap_endian=swap_endian
  dat = assoc(ilun, fltarr(Nx,Ny,/nozer), 512) ; Header is 512 bytes
  Nframes = round(product(dims[2:*]))
  for iframe = 0, Nframes-1 do begin
    red_progressbar, iframe, Nframes, /predict, 'Copying frames'
    self -> fitscube_addframe, fileassoc, dat[iframe], iframe = iframe
  endfor                        ; iframe
  free_lun, ilun                ; Close input file

  
  ;; Close output file and make a flipped version.
  self -> fitscube_finish, olun, flipfile = flipfile, wcs = wcs

;  ;; Add cavity maps as WAVE distortions 
;  if ~keyword_set(nocavitymap) then self -> fitscube_addcmap, oname, cavitymaps

  ;; Add some variable keywords
  self -> fitscube_addvarkeyword, oname, 'DATE-BEG', date_beg_array $
                                  , comment = 'Beginning of observation' $
                                  , keyword_value = self.isodate + 'T' + red_timestring(min(tbeg_array)) $
                                  , axis_numbers = [3, 5] 
  self -> fitscube_addvarkeyword, oname, 'DATE-END', date_end_array $
                                  , comment = 'End time of observation' $
                                  , keyword_value = self.isodate + 'T' + red_timestring(max(tend_array)) $
                                  , axis_numbers = [3, 5] 
  self -> fitscube_addvarkeyword, oname, 'DATE-AVG', date_avg_array $
                                  , comment = 'Average time of observation' $
                                  , keyword_value = self.isodate + 'T' + red_timestring(mean(tavg_array)) $
                                  , axis_numbers = [3, 5] 

  ;; Add variable keywords.
  self -> fitscube_addvarkeyword, oname $
                                  , 'SCANNUM', comment = 'Scan number' $
                                  , s_array, keyword_value = s_array[0] $
                                  , axis_numbers = 5

  self -> fitscube_addvarkeyword, oname $
                                  , 'XPOSURE', comment = '[s] Total exposure time' $
                                  , tunit = 's' $
                                  , exp_array, keyword_value = mean(exp_array) $
                                  , axis_numbers = [3, 5] 

  self -> fitscube_addvarkeyword, oname $
                                  , 'TEXPOSUR', comment = '[s] Single-exposure time' $
                                  , tunit = 's' $
                                  , sexp_array, keyword_value = mean(sexp_array) $
                                  , axis_numbers = [3, 5] 

  self -> fitscube_addvarkeyword, oname $
                                  , 'NSUMEXP', comment = 'Number of summed exposures' $
                                  , nsum_array, keyword_value = mean(nsum_array) $
                                  , axis_numbers = [3, 5] 

  ;; Copy some variable-keywords from the ordinary nb cube to the
  ;; flipped version.
  self -> fitscube_addvarkeyword, flipfile, 'SCANNUM',  old_filename = oname, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'ATMOS_R0', old_filename = oname, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'DATE-BEG', old_filename = oname, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'DATE-AVG', old_filename = oname, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'DATE-END', old_filename = oname, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'XPOSURE',  old_filename = oname, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'TEXPOSUR', old_filename = oname, /flipped
  self -> fitscube_addvarkeyword, flipfile, 'NSUMEXP',  old_filename = oname, /flipped

end
