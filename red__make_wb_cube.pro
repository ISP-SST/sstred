; docformat = 'rst'

;+
; Make a de-rotated and de-stretched time-series FITS data cube with
; momfbd-restored wide-band images.
;
; The header should have information that can be used by companion
; method make_nb_cube to make a de-rotated and de-stretched
; time-series cube with narrow-band images.
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
;     dirs : in, type=strarr
; 
;       The directories where the momfbd output is stored. 
; 
; 
; :Keywords:
;
;    align_interactive : in, optional, type=boolean
;
;      Set this keyword to define the alignment FOV by use of the XROI
;      GUI.
;
;    autocrop : in, optional, type=booean
;
;      Try to determine the largest FOV that avoids any bad momfbd
;      subfields along the edges. If this keyword is set, the input
;      value of the crop keyword is ignored and is set to the
;      auto-detected crop parameters.
;
;    clip : in, optional, type=array 
;
;      Successive clips to use when calculating stretch vectors. See
;      red_destretch_tseries. 
;
;    crop : in, out, optional, type=array, default="[0,0,0,0]"
;
;      The array given here will be used to limit the FOV to
;      [xl+crop[0],xh-crop[1],yl+[crop[3],yh-crop[3]]. If /autocrop,
;      then the auto detected crop is returned in this keyword
;      instead.
;
;    direction : in, optional, type=integer, default="from config file"
;
;      The relative orientation of reference cameras of different
;      instruments. Note that only if direction is set in the config
;      file will it be assumed to be correct when setting the CSYERR
;      FITS header keyword.
;
;    integer : in, optional, type=boolean
;
;      Store as integers instead of floats. Uses the BZERO and BSCALE
;      keywords to preserve the intensity scaling.
;
;    interactive : in, optional, type=boolean
;
;      Set this keyword to define the data cube FOV by use of the XROI
;      GUI. If autocrop is set, then use the so defined FOV as an
;      initialization of the FOV in the GUI. Otherwise use the crop
;      keyword (or its default).
;
;    limb_data : in, optional, type=boolean
;
;      Set for data where the limb is in the FOV. Makes sure to
;      measure the median intensity for normalization on disk. Also
;      disables autocrop.
;
;    nametag : in, optional, type=string
;
;      This string is incorporated into the automatically generated
;      output file name.
;
;    negang : in, optional, type=boolean 
;
;      Set this to apply the field rotation angles with the opposite
;      sign.
;
;    nomissing_nans : in, optional, type=boolean 
;
;      Do not set missing-data padding to NaN. (Set it to the median of
;      each frame instead.)
;
;    nochangesize : in, optional, type=boolean
;
;      Do not increase array size to make room for rotation. Useful
;      for approximately circular FOV.
;
;    nostretch : in, optional, type=boolean
;   
;      Compute no temporal stretch vectors if this is set.
;
;    np : in, optional, type=integer, default=3
;
;      Length of subcubes to use for alignment. See red_aligncube.
;
;    oldname : in, optional, type=boolean
;
;      For data from a single datestamp directory, construct the
;      filename as before multiple directories were implemented.
;
;    point_id : in, optional, type=string, default="From first file"
;
;      Value for the POINT_ID header keyword. By default the value if
;      POINT_ID in the temporally first file or, if that does not
;      exist, the value of DATE-OBS in the temporally first file.
;
;    rotation : in, optional, type=float
;
;      Offset angle to be added to the field rotation angles.
;
;    scannos : in, optional, type=strarr, default="*"
;
;       Choose scan numbers to include in the sequence by entering a
;       comma-and-dash delimited string per input directory, like
;       '2-5,7-20,22-30' or the string '*' to include all. Each
;       element in scannos refers to the corresponding element in the
;       dirs parameter and should have the same number of elements. A
;       scalar scannos (like "0" or "*" ) is repeated for all
;       directories. 
;
;    subtract_meanang : in, optional, type=boolean
;
;      Subtract the mean from the derotation angle. 
;
;    tile : in, optional, type=array
;
;       Successive tiles to use when calculating stretch vectors. See
;       red_destretch_tseries. 
;
;    tstep : in, optional, type=integer, default="3 min equivalent" 
;
;      The number of time steps used as a window when doing unsharp
;      masking in the destretching. See red_destretch_tseries.
;
;    xbd : in, optional, type=integer, default=256
;
;      The X size in pixels of the FOV used for alignment. See
;      red_aligncube. 
;
;    ybd : in, optional, type=integer, default=256
;
;      The Y size in pixels of the FOV used for alignment. See
;      red_aligncube.  
;
; 
; 
; :History:
; 
;    2017-08-16 : MGL. First version, based on code from
;                 chromis::polish_tseries. 
; 
;    2017-09-07 : MGL. Add WCS coordinates and some variable-keywords
;                 to the file. Changed red_fitsaddpar -->
;                 red_fitsaddkeyword.   
;
;    2017-09-28 : MGL. WCS coordinates as a single struct parameter to
;                 fitscube_addwcs.
;
;    2017-10-18 : MGL. Use new keyword dimensions in call to method
;                 fitscube_addwcs. 
;
;    2017-10-27 : MGL. New keyword scannos. 
;
;    2017-11-02 : MGL. New keyword autocrop. Remove keywords square
;                 and origsize. 
;
;    2017-11-08 : MGL. New keyword interactive. 
;
;    2017-11-15 : MGL. Scale data to make use of dynamic range.
; 
;    2018-01-12 : MGL. Use red_bad_subfield_crop.
; 
;    2018-02-08 : MGL. Get logged diskpos (pig or turret) rather than
;                 just pig data.
; 
;    2018-05-03 : MGL. New keyword limb_data. 
; 
;    2018-06-14 : MGL. New keyword align_interactive.  
; 
;    2018-06-20 : MGL. Variable-keywords DATA???.
; 
;    2020-03-11 : MGL. New keywords direction and nametag. 
; 
;    2020-03-12 : MGL. New keyword no_subtract_meanang.
; 
;    2020-03-25 : MGL. Support chromis data.
; 
;    2020-03-31 : MGL. New keyword subtract_meanang, remove keyword
;                 no_subtract_meanang. 
; 
;    2020-04-07 : MGL. New keyword rotation, remove keyword
;                 offset_angle. 
; 
;    2020-06-22 : MGL. Append angles to the "full" keyword when
;                 calling red_rotation. 
; 
;    2020-06-29 : MGL. New keyword nomissing_nans.
; 
;    2020-07-08 : MGL. New keyword integer. 
; 
;    2020-10-28 : MGL. Remove statistics calculations.
; 
;    2020-11-09 : MGL. New keyword nostretch.
; 
;    2021-11-23 : MGL. Combine data from multiple datestamp
;                 directories. New keywords oldname and point_id.
;
;    2022-04-08 : OA. Added time-dependent solar coordinates in WCS.
; 
;    2022-09-04 : MGL. New keyword nochangesize.
; 
;    2022-09-26 : MGL. New keyword rotmargin.
;
;-
pro red::make_wb_cube, dirs $
                       , align_interactive = align_interactive $                       
                       , autocrop = autocrop $
                       , clip = clip $
                       , crop = crop $
                       , direction = direction $
                       , integer = integer $
                       , interactive = interactive $
                       , limb_data = limb_data $
                       , nametag = nametag $
                       , nearest = nearest $
                       , negang = negang $
                       , nomissing_nans = nomissing_nans $
                       , nochangesize = nochangesize $
                       , nostretch = nostretch $
                       , np = np $
                       , nthreads = nthreads $
                       , ofile = ofile  $
                       , oldname = oldname $
                       , point_id = point_id $
                       , rotation = rotation $
                       , rotmargin = rotmargin $
                       , scannos = scannos $
                       , subtract_meanang = subtract_meanang $
                       , tile = tile $
                       , tstep = tstep $
                       , xbd = xbd $
                       , ybd = ybd

  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  ;; dirs and scannos should be strarrs of the same length (or scannos
  ;; not given)
  nDirs = n_elements(dirs)
  if nDirs eq 0 then begin
    print, inam + ' : Please specify the directory with momfbd output.'
    retall
  endif

  case n_elements(scannos) of
    nDirs :                                    ; Continue...
    0 : scannos = replicate('*', nDirs)        ; Default is all scans
    1 : scannos = replicate(scannos[0], nDirs) ; Repeat scalar value
    else : stop                                ; We don't want scannos of other lengths
  endcase
  
  ;; Sort the dirs
  indx = sort(dirs)
  dirs = dirs[indx]
  scannos = scannos[indx]
  
  if(~keyword_set(nearest)) then lin = 1 else lin = 0
  
  
  ;; Name of the instrument
  instrument = ((typename(self)).tolower())

  
  if n_elements(direction) eq 0 then direction = self.direction
  if n_elements(rotation)  eq 0 then rotation  = self.rotation
  if n_elements(rotmargin)  eq 0 then rotmargin  = 40
  
  ;; Make prpara
  red_make_prpara, prpara, align_interactive
  red_make_prpara, prpara, clip
  red_make_prpara, prpara, crop
  red_make_prpara, prpara, dirs    
  red_make_prpara, prpara, direction    
  red_make_prpara, prpara, integer
  red_make_prpara, prpara, rotmargin 
  red_make_prpara, prpara, negang  
  red_make_prpara, prpara, nomissing_nans
  red_make_prpara, prpara, nostretch
  red_make_prpara, prpara, np
  red_make_prpara, prpara, point_id 
  red_make_prpara, prpara, rotation
  red_make_prpara, prpara, scannos
  red_make_prpara, prpara, subtract_meanang  
  red_make_prpara, prpara, tile
  red_make_prpara, prpara, tstep
  red_make_prpara, prpara, xbd
  red_make_prpara, prpara, ybd

  if n_elements(clip) eq 0 then clip = [12,  7,  4,  2, 1]
  if n_elements(tile) eq 0 then tile = [12, 24, 36, 48, 64]

  if keyword_set(limb_data) then autocrop = 0

  ;; Camera/detector identification
  self -> getdetectors
  wbindx = where(strmatch(*self.cameras, instrument+'-W', /fold_case))
  wbcamera = (*self.cameras)[wbindx[0]]
  wbdetector = (*self.detectors)[wbindx[0]]

  ;; Get metadata from logfiles
  red_logdata, self.isodate, time_r0, r0 = metadata_r0, ao_lock = ao_lock
  red_logdata, self.isodate, time_pointing, diskpos = metadata_pointing, rsun = rsun
  red_logdata, self.isodate, time_turret, azel = azel

  ;; Search for restored WB images
  case self.filetype of
    'ANA': extension = '.f0'
    'MOMFBD': extension = '.momfbd'
    'FITS': extension = '.fits'
  endcase

  for idir = 0, nDirs-1 do begin
    if scannos[idir] eq '*' then srch = '*' $
    else srch = '*_' + string(red_expandrange(scannos[idir]), format='(I05)') + '_*'
    fls = file_search(dirs[idir] + srch + extension, count = Nfls)
    if Nfls gt 0 then red_append, files, fls
  endfor
  Nfiles = n_elements(files)
  
;  if Nfiles eq 0 then begin
;    print, inam + ' : No files matching regexp: ' + dir + wbdetector + '*' + extension
;    s = ''
;    read, 'Do you want to make a raw wb cube? [yN] ', s
;    if strupcase(strmid(s, 0, 1)) eq 'Y' then begin
;      make_raw = 1
;    endif else retall
;  endif

  if keyword_set(make_raw) then begin
    stop
    ;;files = red_raw_search('data/'+ dir + '/*-W/', instrument = instrument, scannos = scannos, count = Nraw)
    files = self -> raw_search(file_dirname((*self.data_dirs)[0]) + '/' + dir + '/'+instrument.capwords()+'-W/' $
                               , scannos = scannos, count = Nraw)
    if Nraw eq 0 then stop
    self -> extractstates, files, states

    prefilters = states.prefilter
    prefilters = prefilters[uniq(prefilters, sort(prefilters))]
    if n_elements(prefilters) ne 1 then stop ;; Fix this later
    prefilter = prefilters[0]

    allscannos = states.scannumber
    allscannos = allscannos[uniq(allscannos, sort(allscannos))]
    Nscans = n_elements(allscannos)

    datestamp = fxpar(red_readhead(files[0]), 'DATE-OBS')

    ;; Now select the first file from each scan (later to be changed to the best?)
    wfiles = strarr(Nscans)
    wstates = states[0:Nscans-1]
    for iscan = 0, Nscans-1 do begin
      indx = where(strmatch(files, '*_'+string(allscannos[iscan], format = '(i05)')+'_*'))
      wfiles[iscan] = files[indx[0]]
      wstates[iscan] = states[indx[0]]
    endfor                      ; iscan
    self -> get_calib, wstates[0], darkdata = dark, gaindata = gain
  endif else begin
    ;; The file names we want should have no tuning info.
    indx = where(~strmatch(files,'*_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_[+-]*'), Nscans)
    if Nscans eq 0 then stop
    wfiles = files[indx]
    self -> extractstates, wfiles, wstates
    
    ;; We have no special state (or absence of state) to identify
    ;; the global WB images but we do know that their exposure times
    ;; are much larger than the ones corresponding to the individual
    ;; NB states.
;  self -> extractstates, files, states
;  windx = where(states.EXPOSURE gt mean(states.EXPOSURE)*1.5)
;  wstates = states[windx]
;  wfiles = files[windx]
;  Nscans = n_elements(windx)

    prefilter = wstates[0].prefilter
    datestamp = fxpar(red_readhead(wfiles[0]), 'STARTOBS')
  endelse

;  ;; Get a subset of the available scans, either through the scannos
;  ;; keyword or by a selection dialogue.
;  for idir = 0, nDirs-1 do begin
;    if ~(scannos[idir] eq '*') then begin
;      ;; Selected a subset through the scannos keyword
;      uscans = red_expandrange(scannos[idir])
;      match2, uscans, wstates.scannumber, scanindx
;      if max(scanindx eq -1) eq 1 then begin
;        print, inam + ' : You asked for scans ' + scannos + '. However, scans ' $
;               + red_collapserange(uscans[where(scanindx eq -1)], ld = '', rd = '') $
;               + ' are not available.'
;        print, 'Please change the scannos keyword to a subset of ' $
;               + red_collapserange(wstates.scannumber, ld = '', rd = '') $
;               + ' and try again.'
;        retall
;      endif
;      Nscans = n_elements(scanindx)
;      wstates = wstates[scanindx]
;      wfiles  = wfiles[scanindx]
;  endif
;  endfor                        ; idir
  
  uscans = wstates.scannumber
  time = strarr(Nscans)
  date = strarr(Nscans)
  tmean = fltarr(Nscans)

  x01y01 = red_bad_subfield_crop(wfiles, crop $
                                 , autocrop = autocrop  $
                                 , direction = direction $
                                 , interactive = interactive)

  x0 = x01y01[0] & x1 = x01y01[1] & y0 = x01y01[2] & y1 = x01y01[3]
  origNx = x1 - x0 + 1
  origNy = y1 - y0 + 1

  ;; Observations metadata varaibles
  tbeg_array     = dblarr(1, Nscans) ; Time beginning for state
  tavg_array     = dblarr(1, Nscans) ; Time average for state
  tend_array     = dblarr(1, Nscans) ; Time end for state
  date_beg_array = strarr(1, Nscans) ; DATE-BEG for state
  date_avg_array = strarr(1, Nscans) ; DATE-AVG for state
  date_end_array = strarr(1, Nscans) ; DATE-END for state
  exp_array      = fltarr(1, Nscans) ; Total exposure time
  sexp_array     = fltarr(1, Nscans) ; Single exposure time
  nsum_array     = lonarr(1, Nscans) ; Number of summed exposures
  date_obs_array = strarr(Nscans)    ; Datasets for each scan

  ;; Read headers to get obs_time and load the images into a cube
  cub = fltarr(origNx, origNy, Nscans)
  for iscan = 0L, Nscans -1 do begin
    
    red_progressbar, iscan, Nscans, 'Read headers and load the images into a cube'

    im = red_readdata(wfiles[iscan], head = hdr)

    red_missing, im, /inplace, missing_type_wanted = 'nan'
;    red_missing, im, nmissing = Nmissing, indx_missing = indx_missing, indx_data = indx_data        
;    bg = !Values.F_NaN
;    if Nmissing gt 0 then im[indx_missing] = bg
    
    if keyword_set(make_raw) then begin
      if size(im,/n_dim) gt 2 then im = im[*, *, 0]
      im -= dark
      im *= gain
      im = red_fillpix(im, nthreads = 4L)
    endif
    im = rotate(temporary(im), direction)
    
    red_fitspar_getdates, hdr $
                          , date_beg = date_beg $
                          , date_end = date_end $
                          , date_avg = date_avg $
                          , count_avg = hasdateavg $
                          , comment_avg = comment_avg
    date_beg_array[0, iscan] = date_beg
    date_end_array[0, iscan] = date_end
    date_avg_array[0, iscan] = date_avg
    tbeg_array[0, iscan] = red_time2double((strsplit(date_beg,'T',/extract))[1])
    tend_array[0, iscan] = red_time2double((strsplit(date_end,'T',/extract))[1])
    tavg_array[0, iscan] = red_time2double((strsplit(date_avg,'T',/extract))[1])

    date_obs_array[iscan] = fxpar(hdr, 'DATE-OBS')
    
    ;; Exposure time
    exp_array[0, iscan]  = fxpar(hdr, 'XPOSURE')
    sexp_array[0, iscan] = fxpar(hdr, 'TEXPOSUR', count = Ntexposure)
    nsum_array[0, iscan] = fxpar(hdr, 'NSUMEXP', count = Nnsumexp)
    if Ntexposure eq 0 && Nnsumexp eq 0 then begin
      sexp_array[0, iscan] = exp_array[0, iscan]
      nsum_array[0, iscan] = 1
    endif
    
    if hasdateavg then begin
      date_avg_split = strsplit(date_avg, 'T', /extract, count = Nsplit)
      ddate = date_avg_split[0]
      if Nsplit gt 1 then ttime = date_avg_split[1] else undefine, ttime
    endif else undefine, ddate, ttime

    if n_elements(ddate) eq 0 then begin
      print, inam+' : No date and time information for scan '+strtrim(uscans[iscan], 2)
      stop
    endif else begin
      date[iscan] = ddate
      time[iscan] = ttime
    endelse

    cub[*, *, iscan] = (temporary(im))[x0:x1, y0:y1]
    
    ;; Measure time-dependent intensity variation (sun moves in the Sky)
    if keyword_set(limb_data) then begin
      ;; For limb data, calculate the median in a strip of pixels with
      ;; a certain width, from the limb inward.
      strip_width = 5.                                                ; Approximate width, 5 arsec
      strip_width = round(strip_width/float(self.image_scale))        ; Converted to pixels
      if iscan gt 0 then wdelete                                      ; Don't use more than one window for the bimodal plots
      bimodal_threshold = cgOtsu_Threshold(cub[*, *, iscan], /PlotIt) ; Threshold at the limb
      disk_mask = cub[*, *, iscan] gt bimodal_threshold               
      strip_mask = dilate(sobel(disk_mask),replicate(1, strip_width, strip_width))  * disk_mask 
      strip_indx = where(strip_mask, Nmask) ; Indices within the strip
      if Nmask eq 0 then stop
      tmean[iscan] = median((cub[*, *, iscan])[strip_indx])
    endif else begin
      tmean[iscan] = median(cub[*,*,iscan])
    endelse
    
  endfor                        ; iscan

  hdr = red_readhead(wfiles[0]) ; Base cube header on first WB file header
  
  ;; Plot the intensity variations
  red_timeplot, red_time2double(time), tmean/mean(tmean) $
                , xtitle = 'Time', ytitle = 'Intensity correction' $
                , psym=-16, /ynozero $
                , xrange = [red_time2double(time[0]), red_time2double(time[-1])] + 20*[-1, 1]
  
;  ;; Normalize intensity
;  tmean = tmean/mean(tmean)
  for iscan = 0L, Nscans - 1 do cub[*,*,iscan] /= tmean[iscan]

  ;; Set aside non-rotated and non-shifted cube (re-use variable cub1)
  cub1 = cub
  ang = red_lp_angles(time, date[0], /from_log, offset_angle = rotation)
  mang = median(ang)
  if keyword_set(subtract_meanang) then ang -= mang
  if keyword_set(negang) then ang = -ang
;  if n_elements(offset_angle) then ang += offset_angle
  
  ;; De-rotate images in the cube, has to be done before we can
  ;; calculate the alignment
  for iscan = 0L, Nscans -1 do begin
    red_progressbar, iscan, Nscans, inam+' : De-rotating images.'

;    red_missing, cub[*,*,iscan] $    
;                 , nmissing = Nmissing, indx_missing = indx_missing, indx_data = indx_data    
    bg = !Values.F_NaN
;    if Nmissing gt 0 then begin    
;      bg = (cub[*,*,iscan])[indx_missing[0]]    
;    endif else begin    
;      bg = median(cub[*,*,iscan])    
;    endelse    

    cub[*,*,iscan] = red_rotation(cub[*,*,iscan], ang[iscan], nthreads=nthreads, background = bg)
  endfor                        ; iscan

  ;; Align cube
  if n_elements(np) eq 0 then np = 3
;  if(~keyword_set(np)) then begin
;    np = 0L
;    prompt = inam +' : Please introduce the factor to recompute the reference image: '
;    read, np, prompt = prompt
;  endif

  if n_elements(xbd) eq 0 then xbd = round(origNx*0.9)
  if n_elements(ybd) eq 0 then ybd = round(origNy*0.9)

  ;; Define (default) alignment FOV
  xc = origNx/2
  yc = origNy/2
  align_size = [xbd, ybd]

  if keyword_set(align_interactive) then begin
    
    print
    print, 'Use the XROI GUI to either modify an initial alignment ROI/FOV or define a new one from scratch.'
    print, 'Select Quit in the File menu. The last ROI is used.'
    print

    ;; Define default roi
    X_in = xc + [-1,  1, 1, -1]*xbd/2 
    Y_in = yc + [-1, -1, 1,  1]*ybd/2
    roiobject_in = OBJ_NEW('IDLgrROI', X_in, Y_in)
    roiobject_in -> setproperty, name = 'Default'
    
    ;; Fire up the XROI GUI.
    dispim = bytscl(red_histo_opt(total(cub, 3)))
    xroi, dispim, regions_in = [roiobject_in], regions_out = roiobject, /block $
          , tools = ['Translate-Scale', 'Rectangle'] $
          , title = 'Modify or define alignment ROI'
    roiobject[-1] -> getproperty, roi_xrange = roi_x
    roiobject[-1] -> getproperty, roi_yrange = roi_y

    obj_destroy, roiobject_in
    obj_destroy, roiobject

    xc = round(mean(roi_x))
    yc = round(mean(roi_y))

    align_size = [round(roi_x[1])-round(roi_x[0]), round(roi_y[1])-round(roi_y[0])]

  endif 
  
  ;; Calculate the image shifts
  shift = red_aligncube(cub, np, xbd = align_size[0], ybd = align_size[1] $
                        , xc = xc, yc = yc, nthreads=nthreads) ;, cubic = cubic, /aligncube)
  
  if keyword_set(nochangesize) then begin
    ;; Red_rotation.pro only uses keyword full (which is set to ff
    ;; when calling from make_*_cube) if it is an array with at least
    ;; 5 elements. So setting it to -1 is the same as letting it stay
    ;; undefined, but it can still be passed on to make_nb_cube.
    ff = -1
    ff = [0, -rotmargin, rotmargin, -rotmargin, rotmargin, 0]
    ;; Possibly change this to take the shifts into account but not
    ;; the angles. Something like ff = [0, mdx0, mdx1, mdy0, mdy1].
  endif else begin
    ;; Get maximum angle and maximum shift in each direction
    maxangle = max(abs(ang))
    mdx0 = reform(min(shift[0,*]))
    mdx1 = reform(max(shift[0,*]))
    mdy0 = reform(min(shift[1,*]))
    mdy1 = reform(max(shift[1,*]))
    ff = [maxangle, mdx0, mdx1, mdy0, mdy1, reform(ang)]
  endelse
  
  ;; De-rotate and shift cube
  bg = median(cub1)  
  bg = !Values.F_NaN
  dum = red_rotation(cub1[*,*,0], full=ff $
                     , ang[0], shift[0,0], shift[1,0], background = bg, nthreads=nthreads)
  nd = size(dum,/dim)
  nx = nd[0]
  ny = nd[1]
  cub = fltarr([nd, Nscans])
  cub[*,*,0] = temporary(dum)
  for iscan=0, Nscans-1 do begin
    red_progressbar, iscan, Nscans $
                     , inam+' : Making full-size cube, de-rotating and shifting.'

;    red_missing, cub1[*,*,iscan] $    
;                 , nmissing = Nmissing, indx_missing = indx_missing, indx_data = indx_data    
;    bg = !Values.F_NaN
;    if Nmissing gt 0 then begin    
;      bg = (cub1[*,*,iscan])[indx_missing[0]]    
;    endif else begin    
;      bg = median(cub[*,*,iscan])    
;    endelse    

    cub[*,*,iscan] = red_rotation(cub1[*,*,iscan], full=ff $
                                  , ang[iscan], shift[0,iscan], shift[1,iscan] $
                                  , background = bg, nthreads=nthreads)
  endfor                        ; iscan
  
  dts = red_time2double(time)
  if n_elements(tstep) eq 0 then begin
    tstep = fix(round(180. / median(abs(dts[0:Nscans-2] - dts[1:*]))))
  endif
  tstep = tstep < (Nscans-1)

  print, inam + ' : Using the following parameters for de-stretching the time-series: '
  print, '   tstep [~3 m. (?)]= ', tstep
  print, '   scale [pixels / arcsec] = ', 1.0/float(self.image_scale)
  print, '   tile = ['+strjoin(string(tile, format='(I3)'),',')+']'
  print, '   clip = ['+strjoin(string(clip, format='(I3)'),',')+']'

  ;; Calculate stretch vectors
  grid = red_destretch_tseries(cub, 1.0/float(self.image_scale), tile, clip, tstep $
                               , nthreads = nthreads, nostretch = nostretch)

  if ~keyword_set(nostretch) then begin
    for iscan = 0L, Nscans - 1 do begin
      red_progressbar, iscan, Nscans, inam+' : Applying the stretches.'
      ;;imm1 = red_stretch_linear(cub[*,*,iscan], reform(grid[iscan,*,*,*]), nthreads = nthreads)
      imm = red_rotation(cub1[*,*,iscan], full=ff $
                         , ang[iscan], shift[0,iscan], shift[1,iscan] $
                         , background = bg, stretch_grid = reform(grid[iscan,*,*,*]) $
                         , nthreads=nthreads, nearest=nearest)
      
      ;;if keyword_set(make_raw) then
      ;;mindx = where(cub[*,*,iscan] eq bg, Nwhere)
      ;;if keyword_set(make_raw) &&
      ;;if Nwhere gt 0 then imm[mindx] = bg ; Ugly fix, red_stretch destroys the missing data for raws?
      cub[*,*,iscan] = imm
    endfor                      ; iscan
  endif

  ;; Prepare for making output file names
  if keyword_set(make_raw) then begin
    odir = self.out_dir + '/cubes_raw/'
  endif else begin
    odir = self.out_dir + '/cubes_wb/'
  endelse
  file_mkdir, odir
  
;  if keyword_set(oldname) and nDirs eq 1 then begin
;    midpart = prefilter + '_' + datestamp + '_scans=' $
;              + red_collapserange(uscans, ld = '', rd = '')
;  endif else begin
;    ;; Start constructing the output file name
    for idir = 0, nDirs-1 do begin
      indx = where(strmatch(wstates.filename, dirs[idir]+'*'))
      scn = red_collapserange(wstates[indx].scannumber,r='',l='')
      tst = (stregex(wstates[indx[0]].filename,'/([0-2][0-9]:[0-5][0-9]:[0-6][0-9])/',/extra,/sub))[1]
;      red_append, midparts, tst+'='+scn
      red_append, timestamps, tst      
      red_append, scannos_actual, scn
    endfor                      ; idir
;    midpart = prefilter + '_' + datestamp + '_' + strjoin(midparts, '_')
;  endelse
    
  ;; midpart something like : 6302_2016-09-19T09:28:36_09:28:36=0,1_09:30:20=0-4
  
  ;; Save WB results as a fits file
  if keyword_set(make_raw) then begin
    datatype = 'raw'
  endif else begin
    datatype = 'corrected'
  endelse
  if n_elements(nametag) eq 0 then begin
    ;;  ofil = 'wb_'+midpart+'_'+datatype+'_im.fits'
    datatags = [datatype, 'im']
  endif else begin
    ;;   ofil = 'wb_'+midpart+'_'+nametag+'_'+datatype+'_im.fits'
    datatags = [nametag, datatype, 'im']
  endelse

  ;; POINT_ID
  date_obs = fxpar(hdr, 'DATE-OBS', count = Ndate_obs)  
  if n_elements(point_id) eq 0 then begin
    point_id = fxpar(hdr, 'POINT-ID', count = Npoint_id)
    if Npoint_id eq 0 then point_id = date_obs
  endif

  if keyword_set(ofile) then $
    ofil = ofile $
  else $
    ofil = red_fitscube_filename('wb' $
                               , prefilter $                               
                               , timestamps $                               
                               , scannos_actual $                               
                               , point_id $                   
                               , datatags = datatags $                               
                               , oldname = oldname $                               
                              )
  print, inam + ' : saving WB corrected cube -> ' + odir + ofil
  
  ;; Add the wavelength and Stokes dimensions
  dims = [nx, ny, 1, 1, Nscans]
  cub = reform(cub, dims, /overwrite)

  if keyword_set(integer) then begin

    ;; Make it an integer FITS file
    
    ;; Calculate BZERO and BSCALE
    datamin = min(cub, /nan)
    datamax = max(cub, /nan)

    arraymin = -32000.
    arraymax =  32000.

    ;; BZERO and BSCALE
    bscale = (datamax-datamin)/(arraymax-arraymin)
    bzero = datamin - arraymin*bscale
    
    ;; Set keywords for rescaled integers
;  red_fitsaddkeyword, outhdr, 'BITPIX', 16
    red_fitsaddkeyword, hdr, 'BZERO',  bzero
    red_fitsaddkeyword, hdr, 'BSCALE', bscale

    ;; Make it a 2-byte integer array
;  cub = fix(round(cub*30000./max(cub)))
    cub = fix(round((temporary(cub)-BZERO)/BSCALE))
  endif
  
  ;; Make header. Start with header from last input file, it has
  ;; info about MOMFBD processing.
  check_fits, cub, hdr, /update ; Get dimensions right
;  red_fitsaddkeyword, hdr, 'DATE', red_timestamp(/iso) $     ; DATE with time
;                      , AFTER = 'SIMPLE' $
;                      , 'Creation UTC date of FITS header' ;

  ;; Delete some keywords that do not make sense for WB cubes.
  red_fitsdelkeyword, hdr, 'FNUMSUM'    ; May add this later
  red_fitsdelkeyword, hdr, 'FRAMENUM'   ; Not applicable here
  red_fitsdelkeyword, hdr, 'STATE'      ; State info is in WCS and FILTER1 keywords
;  red_fitsdelkeyword, hdr, '' 

  ;; Some old momfbd output could have an old version of the PRSTEP1
  ;; keyword. Repair that here.
  prstp = fxpar(hdr, 'PRSTEP1', comment = prcom)
  if prstp eq 'MOMFBD image restoration' then begin
    fxaddpar, hdr, 'PRSTEP1', 'MOMFBD', strtrim(prcom, 2)
  endif

  anchor = 'DATE'

  ;; Add some keywords
  red_fitsaddkeyword, anchor = anchor, hdr, 'OBS_HDU', 1

  red_fitsaddkeyword, anchor = anchor, hdr, 'POINT_ID', point_id
  red_fitsaddkeyword, anchor = anchor, hdr, 'INFILES', strjoin(wfiles,','), 'Concatenated files'
  
  ;; Add info to headers
  red_fitsaddkeyword, anchor = anchor, hdr, 'BUNIT', 'dn', 'Units in array: digital number'
  red_fitsaddkeyword, anchor = anchor, hdr, 'BTYPE', 'Intensity', 'Type of data in array'

  ;; Cadence info.
  if Nscans gt 1 then begin
    dt = float((red_differential(dts))[1:*])
    red_fitsaddkeyword, anchor = anchor, hdr, 'CADAVG', mean(dt), 'Average of actual cadence'
    ;; If we can access the planned cadence, it should go into keyword CADENCE.
    if Nscans gt 2 then begin
      red_fitsaddkeyword, anchor = anchor, hdr, 'CADMIN', min(dt),      'Minimum of actual cadence'
      red_fitsaddkeyword, anchor = anchor, hdr, 'CADMAX', max(dt),      'Maximum of actual cadence'
      red_fitsaddkeyword, anchor = anchor, hdr, 'CADVAR', stddev(dt)^2, 'Variance of actual cadence'
    endif
  endif

  ;; Add info about this step
  self -> headerinfo_addstep, hdr $
                              , prstep = 'CONCATENATION,SPATIAL-ALIGNMENT,DESTRETCHING' $
                              , prpara = prpara $
                              , prproc = inam

  ;; Make time tabhdu extension with Nscans rows
  s_array = lonarr(Nscans)
  s_array[0] = wstates.scannumber
  t_array = dblarr(1, Nscans)
  t_array[0] = red_time2double(time) ; In [s] since midnight
  
  wcs = replicate({ wave:dblarr(2,2) $
                  , hplt:dblarr(2,2) $
                  , hpln:dblarr(2,2) $
                  , time:dblarr(2,2) $
                  }, 1, Nscans)
  
  ;; TIME reference value, all times are seconds since midnight.
  dateref = fxpar(hdr, 'DATEREF', count = count)
  if count eq 0 then begin
    ;; One should really check also for the existence of MJDREF
    ;; and JDREF but within the pipeline we can be sure we don't
    ;; use them.
    dateref = self.isodate+'T00:00:00.000000' ; Midnight
    red_fitsaddkeyword, hdr, 'DATEREF', dateref, 'Reference time in ISO-8601', after = 'DATE'
  endif

;  help, round(cub)
  print, 'n_elements:', n_elements(cub)
  ;; Write the file
  self -> fitscube_initialize, odir + ofil, hdr, lun, fileassoc, dims 

  ;; The prefilter is the same for the whole cube.
  wcs[*, *].wave = float(prefilter)/10.

  ;; Get pointing at center of FOV
  red_wcs_hpl_coords, t_array[0, *], metadata_pointing, time_pointing $
                      , hpln, hplt
  
  ;; The alignment routine (red_aligncube) subtracts the median of the
  ;; cross-correlation measured image shifts but removes no trends. So
  ;; it really tries to make the pointing the same during the whole
  ;; sequence without allowing for drifts. So we should make the
  ;; pointing metadata constant in time, let's use the median:
;  hpln = median(hpln)
;  hplt = median(hplt)

  ;; Let's correct coordinates from the log with calculated shifts.
  ;; (It's mostly needed in case of 'jumps'.)
  hpln += double(self.image_scale) * shift[0,*]
  hplt += double(self.image_scale) * shift[1,*]

  ;; Let's smooth coordinates.
  dt = (t_array[0,-1] - t_array[0,0]) / 60. ; minutes
  if dt le 15. or Nscans le 3 then fit_expr = 'P[0] + X*P[1]'
  if dt gt 15. and Nscans gt 3 then fit_expr = 'P[0] + X*P[1] + X*X*P[2]'
  pp = mpfitexpr(fit_expr, t_array[0,*], hpln)
  hpln = red_evalexpr(fit_expr, t_array[0,*], pp)
  pp = mpfitexpr(fit_expr, t_array[0,*], hplt)
  hplt = red_evalexpr(fit_expr, t_array[0,*], pp)
  
  ;; But what we want to tabulate is the pointing in the corners of
  ;; the FOV. Assume hpln and hplt are the coordinates of the center
  ;; of the FOV.
  wcs.hpln[0, 0] = hpln - double(self.image_scale) * (Nx-1)/2.d
  wcs.hpln[1, 0] = hpln + double(self.image_scale) * (Nx-1)/2.d
  wcs.hpln[0, 1] = hpln - double(self.image_scale) * (Nx-1)/2.d
  wcs.hpln[1, 1] = hpln + double(self.image_scale) * (Nx-1)/2.d
  
  wcs.hplt[0, 0] = hplt - double(self.image_scale) * (Ny-1)/2.d
  wcs.hplt[1, 0] = hplt - double(self.image_scale) * (Ny-1)/2.d
  wcs.hplt[0, 1] = hplt + double(self.image_scale) * (Ny-1)/2.d
  wcs.hplt[1, 1] = hplt + double(self.image_scale) * (Ny-1)/2.d
  
  
  for iscan = 0, Nscans-1 do begin
    red_fitscube_addframe, fileassoc, cub[*, *, 0, 0, iscan] $
                           , iscan = iscan
    wcs[0, iscan].time = t_array[iscan]    
  endfor                        ; iscan
;  free_lun, lun
;  print, inam + ' : Wrote file '+odir + ofil
  ;; Close fits file 
  self -> fitscube_finish, lun, wcs = wcs, direction = direction
  
  if ~keyword_set(nomissing_nans) and ~keyword_set(integer) then begin
    ;; Set padding pixels to missing-data, i.e., NaN.
    self -> fitscube_missing, odir + ofil, /noflip, missing_type = 'nan'
  endif else begin
    self -> fitscube_missing, odir + ofil, /noflip, missing_type = 'median'
  endelse

  print, inam + ' : Add calibration data to file '+odir + ofil
  fxbhmake, bhdr, 1, 'MWCINFO', 'Info from make_wb_cube'
;  x01y01 = [X0, X1, Y0, Y1]
  fxbaddcol, col, bhdr, ANG,       'ANG'
  fxbaddcol, col, bhdr, CROP,      'CROP'
  fxbaddcol, col, bhdr, FF,        'FF'
  fxbaddcol, col, bhdr, GRID,      'GRID'
  fxbaddcol, col, bhdr, ND,        'ND'
  fxbaddcol, col, bhdr, SHIFT,     'SHIFT'
  fxbaddcol, col, bhdr, TMEAN,     'TMEAN'
  fxbaddcol, col, bhdr, x01y01,    'X01Y01'
  fxbaddcol, col, bhdr, DIRECTION, 'DIRECTION'

  fxbaddcol, wfiles_col, bhdr, WFILES, 'WFILES'

  fxbcreate, bunit, odir + ofil, bhdr

  fxbwritm, bunit $
            , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01', 'DIRECTION'] $
            ,   ANG,   CROP,   FF,   GRID,   ND,   SHIFT,   TMEAN,   x01y01,   direction
  fxbwrite, bunit, wfiles, wfiles_col, 1
  
  fxbfinish, bunit
  
  print, inam + ' : Added calibration data to file '+odir + ofil

  if 0 then begin
    ;; To read the extension:
    fxbopen, bunit, odir+ofil, 'MWCINFO', bbhdr
    fxbreadm, bunit, row = 1 $
              , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01', 'DIRECTION'] $
              ,  rANG,  rCROP,  rFF,  rGRID,  rND,  rSHIFT,  rTMEAN,  rX01Y01,  rdirection
    ;; Note that the strarr wfiles cannot be read by fxbreadm!
    fxbread, bunit, rWFILES, 'WFILES', 1
    fxbclose, bunit
  endif
  ;; Sample saved quantities with comments about their use in make_crispex:
  ;;
  ;; TSTEP           LONG      =            5          [Not used]
  ;; CLIP            INT       = Array[4]              [Not used]
  ;; TILE            INT       = Array[4]              [Not used]
  ;; SCALE           FLOAT     =       26.3852         [Not used]
  ;; ANG             DOUBLE    = Array[5]              [Yes]
  ;; SHIFT           DOUBLE    = Array[2, 5]           [Yes]
  ;; GRID            FLOAT     = Array[5, 2, 40, 25]   [Yes]
  ;; TIME            STRING    = Array[5]              [Not used]
  ;; DATE            STRING    = Array[5]              [Not used]
  ;; WFILES          STRING    = Array[5]              [Yes]
  ;; TMEAN           FLOAT     = Array[5]              [Yes]
  ;; CROP            INT       = Array[4]              [Not used]
  ;; MANG            DOUBLE    =        6.2653633      [Not used]
  ;; X0              LONG      =            0          [Yes]
  ;; X1              LONG      =         1831          [Yes]
  ;; Y0              LONG      =            0          [Yes]
  ;; Y1              LONG      =         1151          [Yes]
  ;; FF              INT       =        0              [Yes] ( = [maxangle, mdx0, mdx1, mdy0, mdy1] if fullframe)
  ;; ND              UNDEFINED = <Undefined>           [Yes] (2d array if fullframe)
  ;;
  ;; (Some of) these quantities should really go into the metadata
  ;; of the FITS file for two reasons: documentation and use in
  ;; make_crispex. The latter so we can avoid using the save file.
  ;;
  ;; Can make_crispex calculate x0,x1,y0,y1 from naxis1 and naxis2 and
  ;; any of the other parameters?

  ;; Add the WCS coordinates
  if self.direction gt 7 || self.direction ne direction then begin
    ;; Direction not set in the config file, assume unknown
    red_fitscube_addwcs, odir + ofil, wcs, dimensions = dims $
                         , /update $
                         , csyer_spatial_value = 120. $ ; 2 arc minutes
                         , csyer_spatial_comment = '[arcsec] Orientation unknown'
  endif else begin
    ;; Direction set in the config file, so orientation should be
    ;; good. But there may still be small rotation arrors and
    ;; substantial translations.
    red_fitscube_addwcs, odir + ofil, wcs, dimensions = dims $
                         , /update $
                         , csyer_spatial_value = 60. $ ; 1 arc minute
                         , csyer_spatial_comment = '[arcsec] Orientation known'
  endelse 

  ;; Add variable keywords.
  self -> fitscube_addvarkeyword, odir + ofil, 'DATE-BEG', date_beg_array $
                                  , anchor = anchor $
                                  , comment = 'Beginning time of observation' $
                                  , keyword_method = 'first' $
                                  , axis_numbers = [3, 5]
  
  self -> fitscube_addvarkeyword, odir + ofil, 'DATE-END', date_end_array $
                                  , anchor = anchor $
                                  , comment = 'End time of observation' $
                                  , keyword_method = 'last' $
                                  , axis_numbers = [3, 5]
  
  self -> fitscube_addvarkeyword, odir + ofil, 'DATE-AVG', date_avg_array $
                                  , anchor = anchor $
                                  , comment = 'Average time of observation' $
                                  , keyword_value = self.isodate + 'T' + red_timestring(mean(tavg_array)) $
                                  , axis_numbers = [3, 5]
  
  self -> fitscube_addvarkeyword, odir + ofil, 'SCANNUM', s_array $
                                  , comment = 'Scan number' $
                                  , anchor = anchor $
                                  , keyword_method = 'first' $
                                  , axis_numbers = 5

  self -> fitscube_addvarkeyword, odir + ofil, 'DATE-OBS', date_obs_array $
                                  , comment = 'Dataset' $
                                  , anchor = anchor $
                                  , keyword_method = 'first' $
                                  , axis_numbers = 5
  
  self -> fitscube_addvarkeyword, odir + ofil, 'XPOSURE', exp_array $
                                  , comment = 'Summed exposure times' $
                                  , anchor = anchor $
                                  , tunit = 's' $
                                  , keyword_method = 'median' $
                                  , axis_numbers = [3, 5] 

  self -> fitscube_addvarkeyword, odir + ofil, 'TEXPOSUR', sexp_array $
                                  , comment = '[s] Single-exposure time' $
                                  , anchor = anchor $
                                  , tunit = 's' $
                                  , keyword_method = 'median' $
                                  , axis_numbers = [3, 5] 

  self -> fitscube_addvarkeyword, odir + ofil, 'NSUMEXP', nsum_array $
                                  , comment = 'Number of summed exposures' $
                                  , anchor = anchor $
                                  , keyword_method = 'median' $
                                  , axis_numbers = [3, 5]
  
  tindx_r0 = where(time_r0 ge min(t_array) and time_r0 le max(t_array), Nt)
  if Nt gt 0 then begin
    r0dims = size(metadata_r0, /dim) ; We have not always had two r0 values.
    if r0dims[0] eq 1 then extra_coordinate1 = [24] else extra_coordinate1 = [24, 8]
    self -> fitscube_addvarkeyword, odir + ofil, 'ATMOS_R0' $
                                    , metadata_r0[*, tindx_r0] $
                                    , anchor = anchor $
                                    , comment = 'Atmospheric coherence length' $
                                    , tunit = 'm' $
                                    , extra_coordinate1 = extra_coordinate1 $     ; WFS subfield sizes 
                                    , extra_labels      = ['WFSZ'] $              ; Axis labels for metadata_r0
                                    , extra_names       = ['WFS subfield size'] $ ; Axis names for metadata_r0
                                    , extra_units       = ['pix'] $               ; Axis units for metadata_r0
                                    , keyword_method = 'mean' $
                                    , time_coordinate = time_r0[tindx_r0] $
                                    , time_unit       = 's'

    self -> fitscube_addvarkeyword, odir + ofil, 'AO_LOCK' $
                                    , ao_lock[tindx_r0] $
                                    , anchor = anchor $
                                    , comment = 'Fraction of time the AO was locking, 2s average' $
                                    , keyword_method = 'mean' $
                                    , time_coordinate = time_r0[tindx_r0] $
                                    , time_unit       = 's'

  endif

  tindx_turret = where(time_turret ge min(t_array) and time_turret le max(t_array), Nt)
  if Nt gt 0 then begin

    self -> fitscube_addvarkeyword, odir + ofil, 'ELEV_ANG' $
                                    , reform(azel[1, tindx_turret]) $
                                    , anchor = anchor $
                                    , comment = 'Elevation angle' $
                                    , keyword_method = 'mean' $
                                    , time_coordinate = time_turret[tindx_turret] $
                                    , time_unit       = 'deg'
    
  end
    
  if 0 then begin
    fname = odir + ofil
    hhh = headfits(fname)
    scn = red_fitsgetkeyword(fname, 'SCANNUM', comment = comment, variable_values = scn_values)
    xps = red_fitsgetkeyword(fname, 'XPOSURE', comment = comment, variable_values = xps_values)
    print, scn, xps
    help, scn_values, xps_values
    print, scn_values.values, xps_values.values
    
    r0 = red_fitsgetkeyword(fname, 'ATMOS_R0', comment = comment, variable_values = r0_values)
    help, r0_values
  endif
  

end

a = crispred(/dev)
dirs = 'momfbd_nopd/09:28:36/6302/cfg/results/'
dirs = 'momfbd_nopd/'+['09:28:36', '09:30:20']+'/6302/cfg/results/'
scannos = ['0,1']
a -> make_wb_cube, dirs, scannos = scannos

stop

if 1 then begin

  a = crispred(/dev)
  dirs = 'momfbd_nopd/'+['09:28:36', '09:30:20']+'/6302/cfg/results/'
  scannos = ['*', '0-4']
  a -> make_wb_cube, dirs, scannos = scannos

  tmp = red_fitsgetkeyword('cubes_wb2/wb_6302_2016-09-19T09:28:36_scans=0,1,0-4_corrected_im.fits', 'SCANNUM',  variable_values = scannum_values)      
  infiles = red_fitsgetkeyword('cubes_wb2/wb_6302_2016-09-19T09:28:36_scans=0,1,0-4_corrected_im.fits', 'INFILES')


  hprint, strsplit(infiles, ',', /extract) + string(reform(scannum_values.values), format = '(i7)')

endif else begin

  a = chromisred(/dev)
  dirs = 'momfbd/'+['10:01:52', '10:03:04', '10:06:06']+'/3950/cfg/results/'
  scannos = ['*', '*', '0-4']
  a -> make_wb_cube, dirs, scannos = scannos

endelse

end
