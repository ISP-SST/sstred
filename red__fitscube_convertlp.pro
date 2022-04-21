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
;    cavitymaps  : in, optional, type="fltarr(Nx,Ny,1,1,Nscans)"
; 
;      A 3D cube with cavity maps, each adapted to the corresponding
;      scan in the fitscube file. Unit is nm.
;
;    direction : in, optional, type=integer, default="from config file"
;
;      The relative orientation of reference cameras of different
;      instruments. Note that only if direction is set in the config
;      file will it be assumed to be correct when setting the CSYERR
;      FITS header keyword.
;
;    dorotate : in, optional, type=boolean
; 
;      We will assume multi-scan cubes are de-rotated to compensate
;      for the alt-az field rotation and also rotated to Solar-N up.
;      Set this keyword to assume only the former and therefore make
;      the rotation to the proper orientation.
; 
;    headerdata : in, optional, type=strarr
;
;      FITS header with keywords to be added to the output file.
;
;    headerfile : in, optional, type=string
;
;      The name of a file where headerdata can be found.
;
;    mirrorx : in, optional, type=boolean
;
;      The cube wb data are mirrored in the X direction compared to
;      the raw data.
;
;    mirrory : in, optional, type=boolean
;
;      The cube wb data are mirrored in the Y direction compared to
;      the raw data.
;
;    outname : in, out, optional, type=string, default = inname+'.fits'
;
;      Where to write the output. A spectral cube might also be
;      written, the file name for this will be generated based on
;      outname. If a variable, will be returned with the actual file
;      name.
;
;    overwrite : in, optional, type=boolean
;
;       Don't care if cube is already on disk, overwrite it
;       with a new version.
;
;    point_id : in, optional, type=string, default="dateTtimestamp"
;
;      Value for the POINT_ID header keyword. 
;
;    rotation : in, optional, type=float
;
;      Offset angle to be added to the field rotation angles.
; 
;    scannumbers : in, optional, type=string
; 
;      The scan numbers of the scans in the cube as a dash- and
;      comma-delimited string.
; 
;    wbimagefile : in, optional, type=string
; 
;      Path to a file with a wideband image to be added as a contect
;      image for a single-scan file.
;
; :History:
; 
;    2017-12-06 : MGL. First version.
; 
;    2018-05-30 : MGL. Works also with CRISP data.
; 
;    2021-11-26 : MGL. New default output dir and use
;                 red_fitscube_filename() for the filename. New
;                 keywords point_id, wbimagefile, direction, rotation.
;                 Remove keyword nostatistics (and statistics
;                 calculations). Various minor bug fixes.
; 
;-
pro red::fitscube_convertlp, inname $
                             , cavitymaps = cavitymaps $
                             , direction = direction $
                             , dorotate = dorotate $                                                          
                             , flip = flip $
;                             , headerdata = headerdata $
                             , headerfile = headerfile $
                             , mirrorx = mirrorx $                             
                             , mirrory = mirrory $                             
                             , outdir = outdir $
                             , outname = outname $
                             , overwrite = overwrite $                             
                             , point_id = point_id $
                             , pref = pref $
                             , rotation = rotation $
                             , scannumbers = scannumbers $
                             , timestamp = timestamp $
;                             , wbcubefile = wbcubefile $
                             , wbimagefile = wbimagefile

  ;; Direction and rotation. Cmap.
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)
  
  if n_elements(direction) eq 0 then direction = self.direction
  if n_elements(rotation)  eq 0 then rotation  = self.rotation

  ;; May want to check the ALIGN_CLIP keyword if the cfg files are
  ;; still around. In case there was some mirroring done to the WB
  ;; object (the first one). Then the (default) direction might have
  ;; to be changed.
  
  ;; Alternatively, read a raw WB image and compare its orientation to
  ;; the input cube!

  
  ;; Make prpara
  red_make_prpara, prpara, inname
  red_make_prpara, prpara, dorotate    
  red_make_prpara, prpara, headerfile       
  red_make_prpara, prpara, outdir      
  red_make_prpara, prpara, outname      
  red_make_prpara, prpara, overwrite       
  
  dir = file_dirname(inname)+'/'
  iname = file_basename(inname)
  
  if n_elements(outdir) eq 0 then outdir = 'cubes_converted/'
  file_mkdir, outdir


  ;; Cameras and detectors
  self->getdetectors
  wbindx = where(strmatch(*self.cameras,'*-W')) ; Matches CRISP and CHROMIS WB cameras
  wbcamera = (*self.cameras)[wbindx[0]]
  wbdetector = (*self.detectors)[wbindx[0]]
  nbindx = where(strmatch(*self.cameras,'*-[NRT]')) ; Matches CRISP and CHROMIS NB cameras
  nbcamera = (*self.cameras)[nbindx[0]]
  nbdetector = (*self.detectors)[nbindx[0]]
  
  red_lp_header, inname, header=header, datatype=datatype, $
                 dims=dims, nx=nx, ny=ny, nt=nt, endian=endian_file

  if strmatch(header,'stokes=\[I,Q,U,V]*') then Nstokes = 4 else Nstokes = 1
  
  ;; Split the input file name into parts that may be used to identify
  ;; some parameters
  iname = file_basename(inname)    
  iname_parts = strsplit(iname, '._', /extract)  
  
  ;; Wideband or not?
  case iname_parts[0] of
    'crispex' : is_wb = 0
    'wb'      : is_wb = 1
    else : stop                 ; Ask user?
  endcase

  if is_wb then begin
    filetype = 'wb'
    Ntunings = 1
  endif else filetype = 'nb'

  
  
  if n_elements(pref) eq 0 then begin
    pref = (stregex(iname,'[._]([0-9][0-9][0-9][0-9])[._]',/sub, /extr))[1]
    if strlen(pref) ne 4 then stop
  endif

  if n_elements(timestamp) eq 0 then begin
    timestamp = stregex(iname,'[0-2][0-9]:[0-9][0-9]:[0-9][0-9]',/extr)
    if timestamp eq '' then stop
  endif

  nbdir = (*self.data_dirs)[(where(strmatch(*self.data_dirs,'*'+timestamp), Nmatch))[0]] $
          + '/' + nbcamera + '/'
  wbdir = (*self.data_dirs)[(where(strmatch(*self.data_dirs,'*'+timestamp), Nmatch))[0]] $
          + '/' + wbcamera + '/' 

  if is_wb then rawdir = wbdir else rawdir = nbdir

  
  ;; Get file names.
  wbfiles  = file_search(wbdir+'*',  count = Nwbfiles)
  if is_wb then begin
    rawfiles = wbfiles
  endif else begin
    rawfiles = file_search(rawdir+'*', count = Nrawfiles)
  endelse
  

  if n_elements(scannumbers) ne 0 then begin
    ;; The scan numbers were specified
    if size(scannumbers, /tname) eq 'STRING' then begin
      s_array = red_expandrange(scannumbers)
    endif else begin
      s_array = scannumbers
    endelse
    Nscans = n_elements(s_array)
  endif else begin
    ;; Can we get the scan numbers from the file name?
    scans = (stregex(iname,'scans=([0-9]+-[0-9]+)',/extract,/subexpr))[1]    
    if scans eq '' then scans = (stregex(iname,'scan=([0-9]+)',/extract,/subexpr))[1]    
    if scans ne '' then begin
      s_array = red_expandrange(scans)
      Nscans = n_elements(s_array)
    endif
  endelse


  ;; Establish the dimensions of the file
  case 1 of
    
    n_elements(Nscans) gt 0 and n_elements(Ntunings) gt 0 and n_elements(Nstokes) gt 0 : begin
      ;; Dimensions set already. Do they match the file?
      if Nscans * Ntunings * Nstokes ne Nt then stop
    end

    n_elements(Nscans) gt 0 and n_elements(Ntunings) gt 0 : begin
      Nstokes =  round(Nt / (long(Nscans) * long(Ntunings)))
      if Nstokes ne 4 and Nstokes ne 1 then stop
    end

    n_elements(Nscans) gt 0 and n_elements(Nstokes) gt 0 : begin
      Ntunings = round(Nt / (long(Nscans) * long(Nstokes)))
    end

    n_elements(Ntunings) gt 0 and n_elements(Nstokes) gt 0 : begin
      Nscans = round(Nt / (long(Ntunings) * long(Nstokes)))
      ;; If we had to calculate Nscans this way, we don't know the
      ;; actual scan numbers. Except if Nscans matches the number of
      ;; available scans from the raw data.
      self -> extractstates, wbfiles, wbstates
      stop
      ;; Get unique scannumbers from wbstates. Are there Nscans of them?
      
      scannumbers_available = wbstates[uniq(wbstates.scannumber, sort(wbstates.scannumber))].scannumber
      if n_elements(scannumbers_available) eq Nscans then begin
        scannumbers = scannumbers_available
      endif else begin
        print, inam + ' : Number of scans in the input cube:', Nscans        
        print, inam + ' : Number of scans available in raw data ' $
               + red_collapserange(scannumbers_available) + ':', n_elements(scannumbers)
        print, inam + ' : Please try again and use keyword "scannumbers".'
        return
      endelse
    end

  endcase

  ;; Dimensions set. Do they match the file?
  if long(Ntunings) * long(Nscans) * long(Nstokes) ne Nt then begin
    print, inam + ' : Dimensions do not match.'
    print, 'Nt       : ', Nt    
    print, 'Ntunings : ', Ntunings    
    print, 'Nstokes  : ', Nstokes    
    print, 'Nscans   : ', Nscans
    print, 'Ntunings x Nstokes x Nscans = ', long(Ntunings) * long(Nscans) * long(Nstokes)
    stop
  endif



  ;; Establish the orientation of the file
  if n_elements(mirrorx) eq 0 or n_elements(mirrory) eq 0 then begin
    cfgfiles = file_search(dir+'/../*mfbd*/'+timestamp+'/'+pref+'/cfg/*cfg', count = Ncfg)
    if Ncfg eq 0 then begin
      cfgfiles = file_search(dir+'/../../*mfbd*/'+timestamp+'/'+pref+'/cfg/*cfg', count = Ncfg)
    endif
    if Ncfg gt 0 then begin
      ;; We can check the direction with the align clip of the WB
      ;; object (need the redux routines for this)
      align_clip = long(strsplit(redux_cfggetkeyword(cfgfiles[0],'object0.channel0.align_clip'), ',', /extract))
      mirrorx = align_clip[0] gt align_clip[1]    
      mirrory = align_clip[2] gt align_clip[3]    
    endif
  endif

  if n_elements(point_id) eq 0 then point_id = self.isodate + 'T' + timestamp
  
  if n_elements(outname) eq 0 then begin
    oname = outdir + '/' $
            + red_fitscube_filename(filetype $
                                    , pref $
                                    , timestamp $
                                    , red_collapserange(s_array, /nobrack) $
                                    , point_id $ 
                                    , datatags = ['fromlp', 'im'])
  endif else begin
    oname = outname
  endelse

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
  
  
  dims = [Nx, Ny, Ntunings, Nstokes, Nscans]

  ;; Select for scan numbers and prefilter in WB files, this is the
  ;; filter that is in the file name.
  self->selectfiles, files = wbfiles, states = wbstates $
                     , scan = s_array $
                     , pref = pref $
                     , count = Nselected $
                     , selected = selected
  wbfiles = wbfiles[selected]
  wbstates = wbstates[selected]

  if is_wb then begin
    rawfiles  = wbfiles
    rawstates = wbstates
    uindx = 0
  endif else begin
    ;; Now find the same selection of rawfiles.
    self->selectfiles, files = rawfiles, states = rawstates $
                     , framenumbers = wbstates.framenumber $
                     , count = Nselected $
                     , selected = selected
    rawfiles = rawfiles[selected]
    rawstates = rawstates[selected]
    ;; Find unique tunings based on the fpi_state in the WB states.
    uindx = uniq(wbstates.fpi_state, sort(wbstates.fpi_state))
  endelse
  
  ;; Now get the actual wavelengths from the corresponding rawfile (WB
  ;; or NB, depending on the type of cube) states. Don't trust that
  ;; the file names are ordered the same way. Select for frame numbers.
  self->selectfiles, files = rawfiles, states = rawstates $
                     , scan = s_array $
                     , framenumbers = wbstates[uindx].framenumber $
                     , count = Nselected $
                     , selected = selected

  if Ntunings ne Nselected then begin
    print, inam + ' : Not using all tunings?'
    stop 
  endif

  ;; Sort in wavelength order
  uindx = selected[sort(rawstates[selected].tun_wavelength)]
  ulambda = rawstates[uindx].tun_wavelength 
  ustates = rawstates[uindx].fullstate
  if ~is_wb then fpi_states = rawstates[uindx].fpi_state


  if n_elements(ulambda) ne Ntunings then stop
  

  
  ;; Get times and some other metadata from file headers.
  t_array = dblarr(Nscans)                  ; WB time
  tbeg_array     = dblarr(Ntunings, Nscans) ; Time beginning for state
  tavg_array     = dblarr(Ntunings, Nscans) ; Time average for state
  tend_array     = dblarr(Ntunings, Nscans) ; Time end for state
  date_beg_array = strarr(Ntunings, Nscans) ; DATE-BEG for state
  date_avg_array = strarr(Ntunings, Nscans) ; DATE-AVG for state
  date_end_array = strarr(Ntunings, Nscans) ; DATE-END for state
  exp_array      = fltarr(Ntunings, Nscans) ; Total exposure time
  sexp_array     = fltarr(Ntunings, Nscans) ; Single exposure time
  nsum_array     = lonarr(Ntunings, Nscans) ; Number of summed exposures
  
  for iscan = 0L, Nscans-1 do begin

    red_progressbar, iscan, Nscans, /predict $
                     , 'Processing file headers for scan='+strtrim(s_array[iscan], 2)

    ;; Get WB times for log file access
    self->selectfiles, files = wbfiles, states = wbstates $
                       , scan = s_array[iscan] $
                       , pref = pref $
                       , count = Nselected $
                       , selected = selected
    these_tavg_array = dblarr(Nselected)
    for iselected = 0L, Nselected-1 do begin
      thishdr = red_readhead(wbfiles[selected[iselected]])
      red_fitspar_getdates, thishdr, date_avg = date_avg 
      these_tavg_array[iselected] = red_time2double((strsplit(date_avg, 'T', /extract))[1])
    endfor                      ; iselected
    t_array[iscan] = mean(these_tavg_array)

    ;; Get the characteristic wavelength from the WB. This should be
    ;; the right thing for CRISP but not necessarily for CHROMIS.
    wavelnth = fxpar(thishdr, 'WAVELNTH')
    waveunit = fxpar(thishdr, 'WAVEUNIT')
    
    for iwav = 0, Ntunings-1 do begin
      if is_wb then begin
        ;; Get WB times for WCS.
        self->selectfiles, files = rawfiles, states = rawstates $
                           , scan = s_array[iscan] $
                           , count = Nselected $
                           , selected = selected
      endif else begin
        ;; Get NB times for WCS, should be based on fpi_state only. LC
        ;; states get the same timestamps.
        self->selectfiles, files = rawfiles, states = rawstates $
                           , scan = s_array[iscan] $
                           , fpi_stat = fpi_states[iwav] $
                           , count = Nselected $
                           , selected = selected
      endelse
      
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

  ;; Rotation angles
  ang = red_lp_angles(t_array, self.isodate, /from_log, offset_angle = rotation)
  mang = mean(ang, /nan)
  old_ang = ang - mang          ; These are most likely the angles used for the cube
  ;; So the angles to use now, for the new cube, is the "rotation" parameter
  ;; angle plus "mang". 
  
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
  red_fitscopykeyword, anchor = anchor, hdr, 'OBSERVER', thishdr
  red_fitscopykeyword, anchor = anchor, hdr, 'OBJECT'  , thishdr
  red_fitscopykeyword, anchor = anchor, hdr, 'CADENCE' , thishdr
  
  red_fitsaddkeyword, anchor = anchor, hdr, 'WAVELNTH', wavelnth, 'Characteristic wavelength'
  red_fitsaddkeyword, anchor = anchor, hdr, 'WAVEUNIT', waveunit, 'Unit for WAVELNTH'

  dateref = self.isodate+'T00:00:00.000000' ; Midnight
  red_fitsaddkeyword, anchor = anchor, hdr, 'DATEREF', dateref, 'Reference time in ISO-8601'

  red_fitsaddkeyword, anchor = anchor, hdr, 'BUNIT', 'dn', 'Units in array: digital number'
  red_fitsaddkeyword, anchor = anchor, hdr, 'BTYPE', 'Intensity', 'Type of data in array'
  
  ;; Add "global" metadata
  red_metadata_restore, hdr, anchor = anchor
  
  ;; Add info about this step
  self -> headerinfo_addstep, hdr $
                              , anchor = 'OBS_HDU' $
                              , prstep = 'DATA-CURATION' $
                              , prpara = prpara $
                              , prproc = inam
  
  ;; Add provided header data
  if n_elements(headerfile) gt 0 then red_metadata_restore, hdr, fname = headerfile, anchor = 'DATE'
  
  wcs = replicate({  wave:dblarr(2,2) $
                     , hplt:dblarr(2,2) $
                     , hpln:dblarr(2,2) $
                     , time:dblarr(2,2) $
                  }, Ntunings, Nscans)
;  red_fitscube_getwcs, filename $
;                       , coordinates = coordinates $
;                       , distortions = distortions
;  stop
  
  maxangle = rotation + mang    ; Unless the time-dependent angle is not already applied. 
  ff = [maxangle, 0., 0., 0., 0., reform(ang)]  

  if keyword_set(dorotate) then begin
    ;; Get the size of the derotated frames

    ;; Read one frame
    red_lpcube_getframe, inname, frame, iframe = 0
    
    ;; Get the frame with size implied by ff
    frame = rotate(temporary(frame), direction)
    frame = red_rotation(frame, full = ff, 0, 0, 0)

    ;; Set the spatial dimensions of the output file.
    dims[0:1] = size(frame, /dim)

  endif
  
  self -> fitscube_initialize, oname, hdr, olun, fileassoc, dims, wcs = wcs
  
  ;; Get metadata from logfiles
  red_logdata, self.isodate, time_r0, r0 = metadata_r0
;  red_logdata, self.isodate, time_pig, pig = metadata_pig, rsun = rsun
;  red_logdata, self.isodate, time_turret, turret = metadata_turret, rsun = rsun
  red_logdata, self.isodate, time_diskpos, diskpos = diskpos, rsun = rsun

  
  
  ;; Get pointing at center of FOV
  ;;red_wcs_hpl_coords, t_array, metadata_turret, time_turret $
  ;;                    , hpln, hplt
  red_wcs_hpl_coords, t_array, diskpos, time_diskpos $
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
    for iwav = 0, Ntunings-1 do begin
      ;; We rely here on hpln and hplt being the first two tabulated
      ;; coordinates. To make this more general, we should get the
      ;; actual indices from the headers. Maybe later...
      wcs[iwav, iscan].wave = ulambda[iwav]*1e9 ;spect_pos[iwav]/10d ; Å --> nm
      wcs[iwav, iscan].time = tavg_array[iwav, iscan]
    endfor                      ; iwav
  endfor                        ; iscan

  
  ;; Copy data frames
  Nframes = round(product(dims[2:*]))
  iframe = 0
  for iscan = 0, Nscans-1 do begin
    for ituning = 0, Ntunings-1 do begin
      for istokes = 0, Nstokes-1 do begin

        red_progressbar, iframe, Nframes, /predict, 'Copying frames'
        
        red_lpcube_getframe, inname, frame, Nscans = Nscans $
                             , iscan = iscan, ituning = ituning, istokes = istokes

        red_missing, frame, /inplace, missing_value = !Values.F_NaN
        
        ;; Undo any mirroring done in the momfbd processing
        if keyword_set(mirrorx) then frame = reverse(frame, 1, /over)        
        if keyword_set(mirrory) then frame = reverse(frame, 2, /over)        

        if keyword_set(dorotate) then begin
          ;; Assuming frame is in the WB original orientation:
          frame = rotate(temporary(frame), direction)
          ;; Assuming ang is the combined temporal derotation and the
          ;; rotation angle to Solar N up:
          frame = red_rotation(frame, full = ff, rotation+mang, 0, 0, background = !Values.F_NaN)
        endif
        
        red_fitscube_addframe, fileassoc, frame $
                               , iscan = iscan, ituning = ituning, istokes = istokes

        iframe++

      endfor                    ; istokes
    endfor                      ; ituning   
  endfor                        ; iscan   

  self -> fitscube_finish, olun, wcs = wcs
  
  if n_elements(cavitymaps) gt 0 then begin
    ;; Add cavity maps as WAVE distortions 
        red_fitscube_addcmap, oname, reform(cavitymaps, Nx, Ny, 1, 1, Nscans)
  endif

  if Nscans eq 1 then axis_numbers = 3 else axis_numbers = [3, 5]
  
  ;; Add some variable keywords
  self -> fitscube_addvarkeyword, oname, 'DATE-BEG', date_beg_array $
                                  , comment = 'Beginning of observation' $
                                  , keyword_value = self.isodate + 'T' + red_timestring(min(tbeg_array)) $
                                  , axis_numbers = axis_numbers
  self -> fitscube_addvarkeyword, oname, 'DATE-END', date_end_array $
                                  , comment = 'End time of observation' $
                                  , keyword_value = self.isodate + 'T' + red_timestring(max(tend_array)) $
                                  , axis_numbers = axis_numbers
  self -> fitscube_addvarkeyword, oname, 'DATE-AVG', date_avg_array $
                                  , comment = 'Average time of observation' $
                                  , keyword_value = self.isodate + 'T' + red_timestring(mean(tavg_array)) $
                                  , axis_numbers = axis_numbers

  ;; Add variable keywords.
  self -> fitscube_addvarkeyword, oname $
                                  , 'SCANNUM', comment = 'Scan number' $
                                  , s_array, keyword_value = s_array[0] $
                                  , axis_numbers = 5

  self -> fitscube_addvarkeyword, oname $
                                  , 'XPOSURE', comment = '[s] Total exposure time' $
                                  , tunit = 's' $
                                  , exp_array, keyword_value = mean(exp_array) $
                                  , axis_numbers = axis_numbers 

  self -> fitscube_addvarkeyword, oname $
                                  , 'TEXPOSUR', comment = '[s] Single-exposure time' $
                                  , tunit = 's' $
                                  , sexp_array, keyword_value = mean(sexp_array) $
                                  , axis_numbers = axis_numbers 

  self -> fitscube_addvarkeyword, oname $
                                  , 'NSUMEXP', comment = 'Number of summed exposures' $
                                  , nsum_array, keyword_value = mean(nsum_array) $
                                  , axis_numbers = axis_numbers 
 
  if keyword_set(flip) then begin
    ;; Make a flipped version
    print, 'Flip it!'
    red_fitscube_flip, oname $
                       , flipfile = flipfile $
                       , overwrite = overwrite
  endif

  if n_elements(wbimagefile) ne 0 then begin
    print, inam+'Add a WB image as an image extension'
    wbim = red_readdata(wbimagefile, h = wbhdr, direction = direction)
    ehdr=wbhdr
    red_fitsdelkeyword, ehdr, 'STATE' ; Tends to be empty, but do check?
    ;; Add some more keywords?
    ;; Check dimensions vs nb frames?
    fxaddpar, ehdr, 'XTENSION', 'IMAGE'
    sxdelpar, ehdr, 'SIMPLE'
    if n_elements(wbim_rot) eq 0 then begin
      check_fits, wbim, ehdr, /update
    endif else begin
      check_fits, wbim_rot, ehdr, /update
    endelse
    fxaddpar, ehdr, 'DATE', red_timestamp(/utc, /iso)
    anchor = 'DATE'
    red_fitsaddkeyword, anchor = anchor, ehdr, 'EXTNAME', 'WBIMAGE', 'Wideband image'
    red_fitsaddkeyword, anchor = anchor, ehdr, 'PCOUNT', 0
    red_fitsaddkeyword, anchor = anchor, ehdr, 'GCOUNT', 1
    if  n_elements(wbim_rot) eq 0 then begin
      writefits, oname, wbim, ehdr, /append
    endif else begin
      writefits, oname, wbim_rot, ehdr, /append
    endelse
  endif

  outname = oname

  print
  print, inam + ' : Input: ' + inname
  print, inam + ' : Output: ' + oname
  if keyword_set(flip) then print, inam + ' :    and  ' + flipfile
  
end

cd, '/scratch/mats/convert_lp/2015-04-27/CRISP'

a =  crispred(/dev, /no)

orig_dir = '/scratch/tlibb/2015-04-27/crispex/10:41:20/'

a -> fitscube_convertlp, scannumbers = '0-99' $ ; raw data scans 0-101, stokes cubes 0-99.
                         , '/scratch/tlibb/2015-04-27/calib_tseries/wb.5876.10:41:20.corrected.icube' $
                         , outname = wbcubename, mirrorx = 1, mirrory = 0, /dorot, /over


stop

a -> fitscube_convertlp, scannumbers = '0-99', mirrorx = 1, mirrory = 0, /dorot, /over $
                         , orig_dir+'crispex.stokes.5876.10:41:20.time_corrected.fcube'
;$
;   , wbcubefile = '/scratch/tlibb/2015-04-27/calib_tseries/wb.5876.10:41:20.corrected.icube'

end

