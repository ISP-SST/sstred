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
;     dir : in, type=string
; 
;       The directory where the momfbd output is stored.
; 
; 
; :Keywords:
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
;    align_interactive : in, optional, type=boolean
;
;      Set this keyword to define the alignment FOV by use of the XROI
;      GUI.
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
;    negang : in, optional, type=boolean 
;
;      Set this to apply the field rotation angles with the opposite
;      sign. 
;
;    np : in, optional, type=integer, default=3
;
;      Length of subcubes to use for alignment. See red_aligncube.
;
;    offset_angle : in, optional, type=float
;
;      Offset angle to be added to the field rotation angles.
;
;    scannos : in, optional, type=string, default="ask"
;
;       Choose scan numbers to include in the sequence by entering a
;       comma-and-dash delimited string, like '2-5,7-20,22-30' or the
;       string '*' to include all.
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
;-
pro red::make_wb_cube, dir $
                       , align_interactive = align_interactive $
                       , autocrop = autocrop $
                       , clip = clip $
                       , crop = crop $
                       , interactive = interactive $
                       , limb_data = limb_data $
                       , negang = negang $
                       , np = np $
                       , offset_angle = offset_angle $
                       , tile = tile $
                       , tstep = tstep $
                       , xbd = xbd $
                       , ybd = ybd $
                       , scannos = scannos 
;                      , ang = ang $
;                      , ext_date = ext_date $
;                      , ext_time = ext_time $
;                      , shift = shift $



  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  if n_elements(dir) eq 0 then begin
    print, inam + ' : Please specify the directory with momfbd output.'
  endif
 
  ;; Make prpara
  red_make_prpara, prpara, dir
  red_make_prpara, prpara, clip
  red_make_prpara, prpara, crop
  red_make_prpara, prpara, negang
  red_make_prpara, prpara, np
  red_make_prpara, prpara, offset_angle
  red_make_prpara, prpara, tile
  red_make_prpara, prpara, tstep
  red_make_prpara, prpara, ybd
  red_make_prpara, prpara, xbd
  red_make_prpara, prpara, align_interactive

  if n_elements(clip) eq 0 then clip = [12,  6,  3,  1]
  if n_elements(tile) eq 0 then tile = [10, 20, 30, 40]

  if keyword_set(limb_data) then autocrop = 0
  
  ;; Camera/detector identification
  self->getdetectors
  wbindx = where(strmatch(*self.cameras,'Crisp-W'))
  wbcamera = (*self.cameras)[wbindx[0]]
  wbdetector = (*self.detectors)[wbindx[0]]
;  nbindx = where(strmatch(*self.cameras,'Chromis-N')) 
;  nbcamera = (*self.cameras)[nbindx[0]]
;  nbdetector = (*self.detectors)[nbindx[0]]
  ;; Should be generalized to multiple NB cameras if CHROMIS gets
  ;; polarimetry. We don't need to identify any PD cameras for
  ;; restored data.

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
  files = file_search(dir + '*' + extension, count = Nfiles)
  if Nfiles eq 0 then begin
    print, inam + ' : No files matching regexp: ' + dir + '*' + extension
    retall
  endif
  ;; We have no special state (or absence of state) to identify
  ;; the global WB images but we do know that their exposure times
  ;; are much larger than the ones corresponding to the individual
  ;; NB states.
  self -> extractstates, files, states
  windx = where(states.EXPOSURE gt mean(states.EXPOSURE)*1.5)
  wstates = states[windx]
  wfiles = files[windx]
  Nscans = n_elements(windx)

  prefilter = wstates[0].prefilter
  datestamp = fxpar(red_readhead(wfiles[0]), 'STARTOBS')

  ;; Get a subset of the available scans, either through the scannos
  ;; keyword or by a selection dialogue.
  if ~(n_elements(scannos) gt 0 && scannos eq '*') then begin
    if n_elements(scannos) gt 0 then begin
      ;; Selected a subset through the scannos keyword
      uscans = red_expandrange(scannos)
      match2, uscans, wstates.scannumber, scanindx
      if max(scanindx eq -1) eq 1 then begin
        print, inam + ' : You asked for scans ' + scannos + '. However, scans ' $
               + red_collapserange(uscans[where(scanindx eq -1)], ld = '', rd = '') $
               + ' are not available.'
        print, 'Please change the scannos keyword to a subset of ' $
               + red_collapserange(wstates.scannumber, ld = '', rd = '') $
               + ' and try again.'
        retall
      endif
      Nscans = n_elements(scanindx)
    endif else begin
      ;; Selection dialogue
      selectionlist = strtrim(wstates[uniq(wstates.scannumber, sort(wstates.scannumber))].scannumber, 2)
      tmp = red_select_subset(selectionlist $
                              , qstring = inam + ' : Select scans:' $
                              , count = Nscans, indx = scanindx)
    endelse
    wstates = wstates[scanindx]
    wfiles  = wfiles[scanindx]
  endif
  
  uscans = wstates.scannumber
  time = strarr(Nscans)
  date = strarr(Nscans)
  tmean = fltarr(Nscans)
  
  red_bad_subfield_crop, wfiles, crop, autocrop = autocrop,  interactive = interactive
  
  hdr = red_readhead(wfiles[0])
  im_dim = fxpar(hdr, 'NAXIS*')
  x0 = crop[0]
  x1 = im_dim[0]-1 - crop[1]
  y0 = crop[2]
  y1 = im_dim[1]-1 - crop[3]
  Nx = x1 - x0 + 1
  Ny = y1 - y0 + 1

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

  ;; Read headers to get obs_time and load the images into a cube
  cub = fltarr(Nx, Ny, Nscans)
  for iscan = 0L, Nscans -1 do begin
    
    red_progressbar, iscan, Nscans, 'Read headers and load the images into a cube'

    im = red_readdata(wfiles[iscan], head = hdr)

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

    ;; Exposure time
    exp_array[0, iscan]  = fxpar(hdr, 'XPOSURE')
    sexp_array[0, iscan] = fxpar(hdr, 'TEXPOSUR')
    nsum_array[0, iscan] = fxpar(hdr, 'NSUMEXP')
    
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

    cub[*, *, iscan] = red_fillpix((temporary(im))[x0:x1, y0:y1], nthreads = 4L)
    
    ;; Measure time-dependent intensity variation (sun moves in the Sky)
    if keyword_set(limb_data) then begin
      ;; For limb data, calculate the median in a strip of pixels with
      ;; a certain width, from the limb inward.
      strip_width = 5.                                                ; Approximate width, 5 arsec
      strip_width = round(strip_width/float(self.image_scale))        ; Converted to pixels
      if iscan gt 0 then wdelete ; Don't use more than one window for the bimodal plots
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
  
  ;; Plot the intensity variations
  cgplot, uscans, tmean/mean(tmean), xtitle = 'Scan number', ytitle = 'Intensity correction', psym=-1, /ynozero

  ;; Normalize intensity
  tmean = tmean/mean(tmean)
  for iscan = 0L, Nscans - 1 do cub[*,*,iscan] /= tmean[iscan]

  ;; Calculate statistics metadata before alignment and rotation so we
  ;; avoid padding.
  DATAMIN  = dblarr(1, Nscans)  ; The minimum data value
  DATAMAX  = dblarr(1, Nscans)  ; The maximum data value
  DATAMEAN = dblarr(1, Nscans)  ; The average data value
  DATARMS  = dblarr(1, Nscans)  ; The RMS deviation from the mean
  DATAKURT = dblarr(1, Nscans)  ; The kurtosis
  DATASKEW = dblarr(1, Nscans)  ; The skewness
  DATAP01  = dblarr(1, Nscans)  ; The 01 percentile
  DATAP10  = dblarr(1, Nscans)  ; The 10 percentile
  DATAP25  = dblarr(1, Nscans)  ; The 25 percentile
  DATAMEDN = dblarr(1, Nscans)  ; The median data value
  DATAP75  = dblarr(1, Nscans)  ; The 75 percentile
  DATAP90  = dblarr(1, Nscans)  ; The 90 percentile
  DATAP95  = dblarr(1, Nscans)  ; The 95 percentile
  DATAP98  = dblarr(1, Nscans)  ; The 98 percentile
  DATAP99  = dblarr(1, Nscans)  ; The 99 percentile
  percentiles = [.01, .10, .25, .50, .75, .90, .95, .98, .99]
  ;; Note that final cube is rescaled as cub =
  ;; fix(round(cub*30000./max(cub))) so we need to do that here as
  ;; well.
  cub1 = double(fix(round(cub*30000./max(cub))))
  for iscan = 0L, Nscans - 1 do begin

    momnt = moment(cub1[*,*,iscan]) ; mean, variance, skewness, kurtosis
    perc = cgpercentiles(cub1[*,*,iscan], percentiles = percentiles)
    
    DATAMIN[0, iscan]  = min(cub1[*,*,iscan])
    DATAMAX[0, iscan]  = max(cub1[*,*,iscan])

    DATAMEAN[0, iscan] = momnt[0]         
    DATARMS[0, iscan]  = sqrt(momnt[1])   
    DATASKEW[0, iscan] = momnt[2]         
    DATAKURT[0, iscan] = momnt[3]         
    
    DATAP01[0, iscan]  = perc[0]          
    DATAP10[0, iscan]  = perc[1]          
    DATAP25[0, iscan]  = perc[2]          
    DATAMEDN[0, iscan] = perc[3]          
    DATAP75[0, iscan]  = perc[4]          
    DATAP90[0, iscan]  = perc[5]          
    DATAP95[0, iscan]  = perc[6]          
    DATAP98[0, iscan]  = perc[7]          
    DATAP99[0, iscan]  = perc[8]          

  endfor                        ; iscan
  
  CUBEMIN  = min(DATAMIN)
  CUBEMAX  = max(DATAMAX)
  
  ;; Accumulate a histogram for the entire cube, use to calculate
  ;; percentiles.
  Nbins = 2L^16                 ; Use many bins!
  binsize = (CUBEMAX - CUBEMIN) / (Nbins - 1.)
  hist = lonarr(Nbins)
  for iscan = 0L, Nscans - 1 do begin
    hist += histogram(cub1[*, *, iscan], min = cubemin, max = cubemax, Nbins = Nbins, /nan)
  endfor                        ; iscan

  momnt = moment(cub1)          ; mean, variance, skewness, kurtosis
  perc = round(cgpercentiles(cub1, percentiles = percentiles))

  CUBEMEAN =  momnt[0]
  CUBERMS  =  sqrt(momnt[1])
  CUBESKEW =  momnt[2]
  CUBEKURT =  momnt[3]

  CUBEP01  =  perc[0]
  CUBEP10  =  perc[1]
  CUBEP25  =  perc[2]
  CUBEMEDN =  perc[3]
  CUBEP75  =  perc[4]
  CUBEP90  =  perc[5]
  CUBEP95  =  perc[6]
  CUBEP98  =  perc[7]
  CUBEP99  =  perc[8]

  if 0 then begin
    ;; Alt computations. These are tested here and needed for the
    ;; narrowband cube where we don't have it all in memory
    ;; simultaneously.
    Npixels = Nx*Ny
    Nframes = Nscans
    altCUBEMEAN = mean(DATAMEAN)
    altCUBERMS  = sqrt(1d/(Nframes*Npixels-1.d) $
                       * ( total( (Npixels-1d)* DATARMS^2 + Npixels* DATAMEAN^2 ) $
                           - Nframes*Npixels * CUBEMEAN^2 ))
    altCUBESKEW = 1.d/(Nframes*Npixels*CUBERMS^3) $
                  * total( Npixels*DATARMS^3*DATASKEW + Npixels*DATAMEAN^3 + 3*(Npixels-1)*DATAMEAN*DATARMS^2) $
                  - CUBEMEAN^3/CUBERMS^3 - 3*(Npixels*Nframes-1d)*CUBEMEAN/(Npixels*Nframes*CUBERMS) 
    altCUBEKURT = 1.d/(Nframes*Npixels*CUBERMS^4) $
                  * total( Npixels*DATARMS^4*(DATAKURT+3) $
                           + 4*Npixels*DATAMEAN*DATARMS^3*DATASKEW $
                           + 6*(Npixels-1)*DATAMEAN^2*DATARMS^2 + Npixels*DATAmean^4 ) $
                  - 4*CUBEMEAN*CUBESKEW/CUBERMS $
                  - 6*(Npixels*Nframes-1d)*CUBEMEAN^2/(Npixels*Nframes*CUBERMS^2) $
                  - CUBEMEAN^4/CUBERMS^4 - 3 

    ;; Base the percentiles on a histogram
    cghistoplot, cub1
    cgplot, /over, [1, 1]*CUBEP01,  !y.crange, color = 'yellow'
    cgplot, /over, [1, 1]*CUBEP10,  !y.crange, color = 'red'
    cgplot, /over, [1, 1]*CUBEP25,  !y.crange, color = 'blue'
    cgplot, /over, [1, 1]*CUBEMEDN, !y.crange, color = 'orange'
    cgplot, /over, [1, 1]*CUBEP75,  !y.crange, color = 'cyan'
    cgplot, /over, [1, 1]*CUBEP90,  !y.crange, color = 'purple'
    cgplot, /over, [1, 1]*CUBEP95,  !y.crange, color = 'green'
    cgplot, /over, [1, 1]*CUBEP98,  !y.crange, color = 'blue'
    cgplot, /over, [1, 1]*CUBEP99,  !y.crange, color = 'red'

    hist_sum = total(hist, /cum) / total(hist)
    
    ALTCUBEP01  = round(cubemin + ((where(hist_sum gt 0.01))[0]+0.5)*binsize)
    ALTCUBEP10  = round(cubemin + ((where(hist_sum gt 0.10))[0]+0.5)*binsize) 
    ALTCUBEP25  = round(cubemin + ((where(hist_sum gt 0.25))[0]+0.5)*binsize) 
    ALTCUBEMEDN = round(cubemin + ((where(hist_sum gt 0.50))[0]+0.5)*binsize) 
    ALTCUBEP75  = round(cubemin + ((where(hist_sum gt 0.75))[0]+0.5)*binsize) 
    ALTCUBEP90  = round(cubemin + ((where(hist_sum gt 0.90))[0]+0.5)*binsize) 
    ALTCUBEP95  = round(cubemin + ((where(hist_sum gt 0.95))[0]+0.5)*binsize) 
    ALTCUBEP98  = round(cubemin + ((where(hist_sum gt 0.98))[0]+0.5)*binsize) 
    ALTCUBEP99  = round(cubemin + ((where(hist_sum gt 0.99))[0]+0.5)*binsize)

    
    print, 'CUBEMEAN',    CUBEMEAN, altCUBEMEAN
    print, 'CUBERMS',     CUBERMS,  altCUBERMS
    print, 'CUBESKEW',    CUBESKEW, altCUBESKEW
    print, 'CUBEKURT',    CUBEKURT, altCUBEKURT
    print, 'ALTCUBEP01',  CUBEP01,  altCUBEP01 
    print, 'ALTCUBEP10',  CUBEP10,  altCUBEP10 
    print, 'ALTCUBEP25',  CUBEP25,  altCUBEP25 
    print, 'ALTCUBEMEDN', CUBEMEDN, altCUBEMEDN 
    print, 'ALTCUBEP75',  CUBEP75,  altCUBEP75 
    print, 'ALTCUBEP90',  CUBEP90,  altCUBEP90 
    print, 'ALTCUBEP95',  CUBEP95,  altCUBEP95 
    print, 'ALTCUBEP98',  CUBEP98,  altCUBEP98 
    print, 'ALTCUBEP99',  CUBEP99,  altCUBEP99

;  print, 'SUM 4th power by frame'
;  rawdatakurt = DATAKURT+3d
;  print, total(total(cub1^4, 1), 1), reform(Npixels*DATARMS^4*rawdatakurt $
;                                            + 4*Npixels*DATAMEAN*DATASKEW*DATARMS^3 $
;                                            + 6*(Npixels-1)*DATAMEAN^2*DATARMS^2 + Npixels*DATAMEAN^4)
  endif
  
  ;; Set aside non-rotated and non-shifted cube (re-use variable cub1)
  cub1 = cub
  
  ang = red_lp_angles(time, date)
  mang = median(ang)
  ang -= mang
  if keyword_set(negang) then ang = -ang
  if n_elements(offset_angle) then ang += offset_angle
  
  ;; De-rotate images in the cube, has to be done before we can
  ;; calculate the alignment
  for iscan = 0L, Nscans -1 do begin
    red_progressbar, iscan, Nscans, inam+' : De-rotating images.'
    cub[*,*,iscan] = red_rotation(cub[*,*,iscan], ang[iscan])
  endfor                        ; iscan

  ;; Align cube
  if n_elements(np) eq 0 then np = 3
;  if(~keyword_set(np)) then begin
;    np = 0L
;    prompt = inam +' : Please introduce the factor to recompute the reference image: '
;    read, np, prompt = prompt
;  endif

  if n_elements(xbd) eq 0 then xbd = round(Nx*0.9)
  if n_elements(ybd) eq 0 then ybd = round(Ny*0.9)

  ;; Define (default) alignment FOV
  xc = Nx/2
  yc = Ny/2
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
    dispim = bytscl(total(cub, 3))
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
                        , xc = xc, yc = yc) ;, cubic = cubic, /aligncube)

  ;; Get maximum angle and maximum shift in each direction
  maxangle = max(abs(ang))
  mdx0 = reform(min(shift[0,*]))
  mdx1 = reform(max(shift[0,*]))
  mdy0 = reform(min(shift[1,*]))
  mdy1 = reform(max(shift[1,*]))
  ff = [maxangle, mdx0, mdx1, mdy0, mdy1]

  ;; De-rotate and shift cube
  dum = red_rotation(cub1[*,*,0], ang[0], shift[0,0], shift[1,0], full=ff)
  nd = size(dum,/dim)
  nx = nd[0]
  ny = nd[1]
  cub = fltarr([nd, Nscans])
  cub[*,*,0] = temporary(dum)
  for iscan=1, Nscans-1 do begin
    red_progressbar, iscan, Nscans $
                     , inam+' : Making full-size cube, de-rotating and shifting.'
    cub[*,*,iscan] = red_rotation(cub1[*,*,iscan], full=ff $
                                  , ang[iscan], shift[0,iscan], shift[1,iscan])
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
  grid = red_destretch_tseries(cub, 1.0/float(self.image_scale), tile, clip, tstep)

  for iscan = 0L, Nscans - 1 do begin
    red_progressbar, iscan, Nscans, inam+' : Applying the stretches.'
    cub[*,*,iscan] = red_stretch(cub[*,*,iscan], reform(grid[iscan,*,*,*]))
  endfor                        ; iscan

  ;; Prepare for making output file names
  odir = self.out_dir + '/cubes_wb/'
  file_mkdir, odir
  midpart = prefilter + '_' + datestamp + '_scans=' $
            + red_collapserange(uscans, ld = '', rd = '')

  ;; Save WB results as a fits file
  ofil = 'wb_'+midpart+'_corrected_im.fits'
  print, inam + ' : saving WB corrected cube -> ' + odir + ofil

  ;; Add the wavelength and Stokes dimensions
  dims = [nx, ny, 1, 1, Nscans]
  cub = reform(cub, dims, /overwrite)

  ;; Make it a 2-byte integer array
  cub = fix(round(cub*30000./max(cub)))

  ;; Make header. Start with header from last input file, it has
  ;; info about MOMFBD processing.
  check_fits, cub, hdr, /update                              ; Get dimensions right
;  red_fitsaddkeyword, hdr, 'DATE', red_timestamp(/iso) $     ; DATE with time
;                      , AFTER = 'SIMPLE' $
;                      , 'Creation UTC date of FITS header' ;
  anchor = 'DATE'

  ;; Add some keywords
  red_fitsaddkeyword, anchor = anchor, hdr, 'OBS_HDU', 1

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
                              , prstep = 'Prepare WB science data cube' $
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
  
  
  for iscan = 0, Nscans-1 do begin
    self -> fitscube_addframe, fileassoc, cub[*, *, 0, 0, iscan] $
                               , iscan = iscan
    wcs[0, iscan].time = t_array[iscan]
  endfor                        ; iscan
  free_lun, lun
  print, inam + ' : Wrote file '+odir + ofil

  print, inam + ' : Add calibration data to file '+odir + ofil
  fxbhmake, bhdr, 1, 'MWCINFO', 'Info from make_wb_cube'
  x01y01 = [X0, X1, Y0, Y1]
  fxbaddcol, col, bhdr, ANG,    'ANG'
  fxbaddcol, col, bhdr, CROP,   'CROP'
  fxbaddcol, col, bhdr, FF,     'FF'
  fxbaddcol, col, bhdr, GRID,   'GRID'
  fxbaddcol, col, bhdr, ND,     'ND'
  fxbaddcol, col, bhdr, SHIFT,  'SHIFT'
  fxbaddcol, col, bhdr, TMEAN,  'TMEAN'
  fxbaddcol, col, bhdr, x01y01, 'X01Y01'

  fxbaddcol, wfiles_col, bhdr, WFILES, 'WFILES'

  fxbcreate, bunit, odir + ofil, bhdr

  fxbwritm, bunit $
            , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01'] $
            ,   ANG,   CROP,   FF,   GRID,   ND,   SHIFT,   TMEAN,   x01y01
  fxbwrite, bunit, wfiles, wfiles_col, 1
  
  fxbfinish, bunit
  
  print, inam + ' : Added calibration data to file '+odir + ofil

  if 0 then begin
    ;; To read the extension:
    fxbopen, bunit, odir+ofil, 'MWCINFO', bbhdr
    fxbreadm, bunit, row = 1 $
            , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01'] $
              ,   rANG,   rCROP,   rFF,   rGRID,   rND,  rSHIFT,   rTMEAN,   rX01Y01
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
  self -> fitscube_addwcs, odir + ofil, wcs, dimensions = dims
  
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
  
  ;; Statistics

  self -> fitscube_addvarkeyword, odir + ofil, 'DATAMIN', DATAMIN $
                                  , comment = 'The minimum data value' $
                                  , anchor = anchor $
                                  , keyword_value = CUBEMIN $
                                  , axis_numbers = [3, 5]


  self -> fitscube_addvarkeyword, odir + ofil, 'DATAMAX', DATAMAX $
                                  , comment = 'The maximum data value' $
                                  , anchor = anchor $
                                  , keyword_value = CUBEMAX $
                                  , axis_numbers = [3, 5]


  self -> fitscube_addvarkeyword, odir + ofil, 'DATAMEAN', DATAMEAN $
                                  , comment = 'The average data value' $
                                  , anchor = anchor $
                                  , keyword_value = CUBEMEAN $
                                  , axis_numbers = [3, 5]

  self -> fitscube_addvarkeyword, odir + ofil, 'DATARMS', DATARMS $
                                  , comment = 'The RMS deviation from the mean' $
                                  , anchor = anchor $
                                  , keyword_value = CUBERMS $
                                  , axis_numbers = [3, 5]

  self -> fitscube_addvarkeyword, odir + ofil, 'DATAKURT', DATAKURT $
                                  , comment = 'The kurtosis' $
                                  , anchor = anchor $
                                  , keyword_value = CUBEKURT $
                                  , axis_numbers = [3, 5]

  self -> fitscube_addvarkeyword, odir + ofil, 'DATASKEW', DATASKEW $
                                  , comment = 'The skewness' $
                                  , anchor = anchor $
                                  , keyword_value = CUBESKEW $
                                  , axis_numbers = [3, 5]

  self -> fitscube_addvarkeyword, odir + ofil, 'DATAP01', DATAP01 $
                                  , comment = 'The 01 percentile' $
                                  , anchor = anchor $
                                  , keyword_value = CUBEP01 $
                                  , axis_numbers = [3, 5]

  self -> fitscube_addvarkeyword, odir + ofil, 'DATAP10', DATAP10 $
                                  , comment = 'The 10 percentile' $
                                  , anchor = anchor $
                                  , keyword_value = CUBEP10 $
                                  , axis_numbers = [3, 5]

  self -> fitscube_addvarkeyword, odir + ofil, 'DATAP25', DATAP25 $
                                  , comment = 'The 25 percentile' $
                                  , anchor = anchor $
                                  , keyword_value = CUBEP25 $
                                  , axis_numbers = [3, 5]

  self -> fitscube_addvarkeyword, odir + ofil, 'DATAMEDN', DATAMEDN $
                                  , comment = 'The median data value' $
                                  , anchor = anchor $
                                  , keyword_value = CUBEMEDN $
                                  , axis_numbers = [3, 5]

  self -> fitscube_addvarkeyword, odir + ofil, 'DATAP75', DATAP75 $
                                  , comment = 'The 75 percentile' $
                                  , anchor = anchor $
                                  , keyword_value = CUBEP75 $
                                  , axis_numbers = [3, 5]

  self -> fitscube_addvarkeyword, odir + ofil, 'DATAP90', DATAP90 $
                                  , comment = 'The 90 percentile' $
                                  , anchor = anchor $
                                  , keyword_value = CUBEP90 $
                                  , axis_numbers = [3, 5]

  self -> fitscube_addvarkeyword, odir + ofil, 'DATAP95', DATAP95 $
                                  , comment = 'The 95 percentile' $
                                  , anchor = anchor $
                                  , keyword_value = CUBEP95 $
                                  , axis_numbers = [3, 5]

  self -> fitscube_addvarkeyword, odir + ofil, 'DATAP98', DATAP98 $
                                  , comment = 'The 98 percentile' $
                                  , anchor = anchor $
                                  , keyword_value = CUBEP98 $
                                  , axis_numbers = [3, 5]

  self -> fitscube_addvarkeyword, odir + ofil, 'DATAP99', DATAP99 $
                                  , comment = 'The 99 percentile' $
                                  , anchor = anchor $
                                  , keyword_value = CUBEP99 $
                                  , axis_numbers = [3, 5]

  
  tindx_r0 = where(time_r0 ge min(t_array) and time_r0 le max(t_array), Nt)
  if Nt gt 0 then begin

    self -> fitscube_addvarkeyword, odir + ofil, 'ATMOS_R0' $
                                    , metadata_r0[*, tindx_r0] $
                                    , anchor = anchor $
                                    , comment = 'Atmospheric coherence length' $
                                    , tunit = 'm' $
                                    , extra_coordinate1 = [24, 8] $               ; WFS subfield sizes 
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
