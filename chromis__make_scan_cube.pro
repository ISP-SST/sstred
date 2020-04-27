; docformat = 'rst'

;+
; Make FITS data cubes with momfbd-restored narrowband images, one
; scan per file.
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
;     clip : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment.
;
;     cmap_fwhm : in, optional, type=float, default=7
;   
;       FWHM in pixels of kernel used for smoothing the cavity map.
;
;     integer : in, optional, type=boolean
;
;       Store as integers instead of floats.
;
;     intensitycorrmethod : in, optional, type="string or boolean", default=FALSE/none
;
;       Indicate whether to do intensity correction based on WB data
;       and with what method. See documentation for red::fitscube_intensitycorr.
;
;     limb_data : in, optional, type=boolean
;
;       Set for data where the limb is in the FOV. Disables autocrop.
;
;     noaligncont : in, optional, type=boolean
;
;       Do not do the align continuum to wideband step.
;
;     nocavitymap : in, optional, type=boolean
;
;       Do not add cavity maps to the WCS metadata.
;       no effect.)
;
;     nostatistics : in, optional, type=boolean
;  
;       Do not calculate statistics metadata to put in header keywords
;       DATA*. If statistics keywords already exist, then remove them.
;
;     odir : in, optional, type=string, detault='cubes_scan/'
;
;       The output directory.
;
;     overwrite : in, optional, type=boolean
;
;       Don't care if cube is already on disk, overwrite it
;       with a new version.
;
;     tile : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment. 
;
; 
; :History:
; 
;    2018-01-19 : MGL. First version, based on code from
;                 chromis::make_wb_cube and chromis::make_nb_cube.
; 
;    2018-01-30 : MGL. Add the corresponding wideband image in an
;                 image extension.
; 
;    2018-02-08 : MGL. Get logged diskpos (pig or turret) rather than
;                 just pig data.
; 
;    2018-05-08 : MGL. New keyword limb_data. 
; 
;    2019-08-26 : MGL. Do integerization, statistics calculations, and
;                 WB intensity correction by calling subprograms. 
; 
;    2020-01-15 : MGL. New keywords intensitycorr and odir. 
; 
;    2020-01-19 : MGL. Rename keyword intensitycorr to
;                 intensitycorrmethod. 
; 
;    2020-04-03 : MGL. New keywords direction and norotation.
; 
;    2020-04-27 : MGL. New keyword rotation.
; 
;-
pro chromis::make_scan_cube, dir $
                             , autocrop = autocrop $
                             , clip = clip $
                             , cmap_fwhm = cmap_fwhm $
                             , crop = crop $
                             , direction = direction $
                             , integer = integer $
                             , intensitycorrmethod = intensitycorrmethod $
                             , interactive = interactive $
                             , limb_data = limb_data $
                             , noaligncont = noaligncont $
                             , nocavitymap = nocavitymap $
                             , norotation = norotation $
                             , nostatistics = nostatistics $
                             , odir = odir $
                             , overwrite = overwrite $
                             , rotation = rotation $
                             , scannos = scannos $
                             , tile = tile 
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  if n_elements(direction) eq 0 then direction = self.direction
  if n_elements(rotation)  eq 0 then rotation  = self.rotation
                             
  ;; Temporarily disable cavity maps by default, can still be be
  ;; written (experimentally) with explicit nocavitymap=0.
  if n_elements(nocavitymap) eq 0 then nocavitymap = 1
  
  ;; Make prpara
  red_make_prpara, prpara, dir
  red_make_prpara, prpara, autocrop
  red_make_prpara, prpara, clip
  red_make_prpara, prpara, crop     
  red_make_prpara, prpara, direction
  red_make_prpara, prpara, integer  
  red_make_prpara, prpara, intensitycorrmethod  
  red_make_prpara, prpara, interactive
  red_make_prpara, prpara, noaligncont 
  red_make_prpara, prpara, nocavitymap  
  red_make_prpara, prpara, norotation       
  red_make_prpara, prpara, overwrite
  red_make_prpara, prpara, rotation       
  red_make_prpara, prpara, tile

  ;; Default keywords
  if n_elements(cmap_fwhm) eq 0 then fwhm = 7.0
  if n_elements(clip) eq 0 then clip = [8, 4,  2,  1,  1  ]
  if n_elements(tile) eq 0 then tile = [8, 16, 32, 64, 128]

  if keyword_set(limb_data) then autocrop = 0

  ;; Output directory
  if(n_elements(odir) eq 0) then odir = self.out_dir + '/cubes_scan/' 

  
  ;; We do currently not correct for the small scale cavity map in
  ;; CHROMIS data. (We should get this from earlier meta data!)
  remove_smallscale = 0       

  ;; Camera/detector identification
  self->getdetectors
  wbindx     = where(strmatch(*self.cameras,'Chromis-W'))
  wbcamera   = (*self.cameras)[wbindx[0]]
  wbdetector = (*self.detectors)[wbindx[0]]
  nbindx     = where(strmatch(*self.cameras,'Chromis-N')) 
  nbcamera   = (*self.cameras)[nbindx[0]]
  nbdetector = (*self.detectors)[nbindx[0]]
  ;; Should be generalized to multiple NB cameras if CHROMIS gets
  ;; polarimetry. We don't need to identify any PD cameras for
  ;; restored data.

  ;; Get metadata from logfiles
  red_logdata, self.isodate, time_r0, r0 = metadata_r0
  red_logdata, self.isodate, time_pointing, diskpos = metadata_pointing, rsun = rsun


  ;; Search for restored images
  case self.filetype of
    'ANA': extension = '.f0'
    'MOMFBD': extension = '.momfbd'
    'FITS': extension = '.fits'
  endcase
  if n_elements(scannos) gt 0 then begin
    files = file_search(dir + '*_'+string(scannos, format = '(i05)')+'_*'+extension, count = Nfiles)      
    if Nfiles eq 0 then begin
      print
      print, inam+' : No momfbd output for scan(s) '+strjoin(scannos, ',')+' found in '+dir
      return
    endif
  endif else begin
    files = file_search(dir + '*'+extension, count = Nfiles)      
    if Nfiles eq 0 then begin
      print
      print, inam+' : No momfbd output found in '+dir
      return
    endif
  endelse
  
  self -> extractstates, files, states
  ;; files & states -> all files in directory

  
  ;; We have no special state (or absence of state) to identify the
  ;; global WB images but we do know that their exposure times are
  ;; much larger than the ones corresponding to the individual NB
  ;; states.
  windx = where(states.EXPOSURE gt mean(states.EXPOSURE)*1.5)
  wstates = states[windx]
  wfiles = files[windx]
  ;; wfiles & wstates -> all global WB files in directory

;  Nscans = n_elements(windx)

  ;; Some info common to all scans
  prefilter = wstates[0].prefilter
  wbghdr = red_readhead(wfiles[0])
  datestamp = fxpar(wbghdr, 'STARTOBS')
  timestamp = (strsplit(datestamp, 'T', /extract))[1]

  red_fitspar_getdates, wbghdr $
                        , date_beg = date_beg $
                        , date_end = date_end $
                        , date_avg = date_avg $
                        , count_avg = hasdateavg $
                        , comment_avg = comment_avg
  
  if hasdateavg then begin
    date_avg_split = strsplit(date_avg, 'T', /extract, count = Nsplit)
    ddate = date_avg_split[0]
    if Nsplit gt 1 then time = date_avg_split[1] else undefine, time
  endif else undefine, ddate, time

  ;; Derotation angle
  ang = (red_lp_angles(time, ddate[0], /from_log, offset_angle = rotation))[0]

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

  ;; wfiles & wstates -> selected and existing global WB files in
  ;; directory. Also uscans is the scans to process and Nscans is the
  ;; number of such scans.



  ;; Establish the FOV, perhaps based on the selected wb files.
  red_bad_subfield_crop, wfiles, crop, autocrop = autocrop,  interactive = interactive
;  hdr = red_readhead(wfiles[0])
  im_dim = fxpar(wbghdr, 'NAXIS*')
  if max(direction eq [1, 3, 4, 6]) eq 1 then begin
    ;; X and Y switched
    y0 = crop[0]
    y1 = im_dim[0]-1 - crop[1]
    x0 = crop[2]
    x1 = im_dim[1]-1 - crop[3]
  endif else begin
    x0 = crop[0]
    x1 = im_dim[0]-1 - crop[1]
    y0 = crop[2]
    y1 = im_dim[1]-1 - crop[3]
  endelse
  Nx = x1 - x0 + 1
  Ny = y1 - y0 + 1
  

  
  ;; Now let's limit the files & states arrays to only the
  ;; scans to process.
  self -> selectfiles, files = files, states = states $
                       , scan = uscans, sel = pindx
  files = files[pindx]
  states = states[pindx]

  ;; files & states -> existing files of selected scans in directory.

  
  ;; Select the nb and wb "per tuning files" by excluding the global
  ;; WB images
  self -> selectfiles, files = files, states = states $
                       , cam = wbcamera, ustat = '' $
                       , sel = wbgindx $
;                       , count = Nscans $
                       , complement = complement, Ncomplement = Ncomplement
  ;; We have no special state (or absence of state) to identify
  ;; the global WB images but we do know that their exposure times
  ;; are much larger than the ones corresponding to the individual
  ;; NB states.
  wbindx = where(states.exposure gt mean(states.exposure)*1.5 $
                 , complement = complement, Ncomplement = Ncomplement) 
  ;; All the per-tuning files and states
  pertuningfiles = files[complement]
  pertuningstates = states[complement]

  ;; pertuningfiles & pertuningstates -> existing files of selected
  ;; scans in directory, excluding the global WB images.

  
    
    
  ;; Continuum alignment only done for Ca II scans (so far). H beta is
  ;; not as wide so should be OK.
  if prefilter eq '3950' && ~keyword_set(noaligncont) then begin
    
    ;; Get wavelength-variable shifts based on continuum vs wideband
    ;; alignment.
    
    aligndir = self.out_dir + '/align/' + timestamp $
               + '/' + prefilter + '/'
    
    nname = aligndir+'scan_numbers.fz'
    sname = aligndir+'continuum_shifts_smoothed.fz'
    
    if ~file_test(nname) || ~file_test(sname) then begin
      print, inam + ' : At least one file missing for aligncont option:'
      print, nname
      print, sname
      retall
    endif
    
    ;; Read the shifts for the continuum images
    align_scannumbers = f0(nname)
    align_shifts = f0(sname)

    ;; Check that we have alignment for all scan numbers
    match2, uscans, align_scannumbers, suba, subb
    missing_indx = where(suba eq -1, Nmissing)
    if Nmissing gt 0 then begin
      print, inam+' : Alignment missing for these scan numbers:'
      print, uscans[missing_indx]
      print, inam+' : Please rerun a -> align_continuum'
      retall
    endif
    
    ;; Select align shifts for the relevant scan numbers.
    nb_shifts = fltarr(2, Nscans)
    nb_shifts[0, *] = align_shifts[0, suba]
    nb_shifts[1, *] = align_shifts[1, suba]
  
;    ;; Use interpolation to get the shifts for the selected scans.
;    nb_shifts = fltarr(2, Nscans)
;    for iscan=0L, Nscans-1 do begin
;      pos = where(align_scannumbers eq uscans[iscan], cccc)
;      if cccc eq 1 then nb_shifts[*, iscan] = align_shifts[*, pos] else begin
;        nb_shifts[0, *] = interpol([reform(align_shifts[0, *])] $
;                                   , [float(align_scannumbers)], [float(uscans)])
;        nb_shifts[1, *] = interpol([reform(align_shifts[1, *])] $
;                                   , [float(align_scannumbers)], [float(uscans)])
;      endelse
;    endfor
    pos = where(~finite(nb_shifts), cccc)
    if cccc gt 0 then nb_shifts[pos] = 0
  endif


  for iscan = 0L, Nscans-1 do begin
    
    ;; This is the loop in which the cubes are written.
    
    ;; Make output file name
    midpart = prefilter + '_' + datestamp + '_scan=' $ 
              + strtrim(uscans[iscan], 2)
    ofile = 'nb_'+midpart+'_corrected.fits'
    filename = odir+ofile

    ;; Already done?
    if file_test(filename) then begin
      if keyword_set(overwrite) then begin
        print, 'Overwriting existing data cube:'
        print, filename
      endif else begin
        print, 'This data cube exists already:'
        print, filename
        continue
      endelse
    endif

    ;; Make the directory if needed.
    file_mkdir, odir

    ;; Unique tuning states, sorted by wavelength
    utunindx = uniq(pertuningstates.fpi_state, sort(pertuningstates.fpi_state))
    Ntuning = n_elements(utunindx)
    sortindx = sort(pertuningstates[utunindx].tun_wavelength)
    ufpi_states = pertuningstates[utunindx[sortindx]].fpi_state
    utunwavelength = pertuningstates[utunindx[sortindx]].tun_wavelength

    wav = utunwavelength
    my_prefilters = pertuningstates[utunindx[sortindx]].prefilter

    ;; Unique nb prefilters
    unbprefindx = uniq(pertuningstates[utunindx].prefilter, sort(pertuningstates[utunindx].prefilter))
    Nnbprefs = n_elements(unbprefindx)
    unbprefs = pertuningstates[utunindx[unbprefindx]].prefilter
;  unbprefsref = dblarr(Nnbprefs)
;
;  for inbpref = 0L, Nnbprefs-1 do begin
;    ;; This is the reference point of the fine tuning for this prefilter:
;    unbprefsref[inbpref] = double((strsplit(pertuningstates[utunindx[unbprefindx[inbpref]]].tuning $
;                                            , '_', /extract))[0])
;  endfor                        ; inbpref
;  
;  unbprefsref *= 1e-10          ; [m]



    ;; Load prefilters
    for inbpref = 0L, Nnbprefs-1 do begin
      pfile = self.out_dir + '/prefilter_fits/chromis_'+unbprefs[inbpref]+'_prefilter.idlsave'
      if ~file_test(pfile) then begin
        print, inam + ' : prefilter file not found: '+pfile
        return
      endif
      
      restore, pfile            ; Restores variable prf which is a struct
      idxpref = where(my_prefilters eq unbprefs[inbpref], count)
      
      if inbpref eq 0 then begin
        units = prf.units
      endif else begin
        if units ne prf.units then begin
          print, inam + ' : Units in ' + pfile + ' do not match those in earlier read files.'
          print, inam + ' : Please rerun the prefilterfit step for these data.'
          retall
        endif
      endelse

      if count eq 1 then begin
        red_append, prefilter_curve, prf.pref
        red_append, prefilter_wav, prf.wav
;        red_append, prefilter_wb, prf.wbint
      endif else begin
        me = median(prf.wav)
        red_append, prefilter_curve, red_intepf(prf.wav-me, prf.pref, wav[idxpref]*1.d10-me)
        red_append, prefilter_wav, wav[idxpref]*1.d10
;        red_append, prefilter_wb, replicate(prf.wbint, count)
      endelse
      
    endfor                      ; inbpref

    rpref = 1.d0/prefilter_curve

    ;; Set up for collecting time and wavelength data
    tbeg_array     = dblarr(Ntuning)   ; Time beginning for state
    tend_array     = dblarr(Ntuning)   ; Time end for state
    tavg_array     = dblarr(Ntuning)   ; Time average for state
    date_beg_array = strarr(Ntuning)   ; DATE-BEG for state
    date_end_array = strarr(Ntuning)   ; DATE-END for state
    date_avg_array = strarr(Ntuning)   ; DATE-AVG for state
    exp_array      = fltarr(Ntuning)   ; Total exposure time
    sexp_array     = fltarr(Ntuning)   ; Single exposure time
    nsum_array     = lonarr(Ntuning)   ; Number of summed exposures

    wcs = replicate({  wave:dblarr(2,2) $
                       , hplt:dblarr(2,2) $
                       , hpln:dblarr(2,2) $
                       , time:dblarr(2,2) $
                    }, Ntuning)

;    ;; Per-tuning files, wb and nb, only for selected scan
;    self -> selectfiles, files = pertuningfiles, states = pertuningstates $
;                         , scan = uscans[iscan] $
;                         , cam = wbcamera $
;                         , sel = wbindx, count = Nwb
;    wbstates = pertuningstates[wbindx]
;    wbfiles = pertuningfiles[wbindx]
;    self -> selectfiles, files = pertuningfiles, states = pertuningstates $
;                         , scan = uscans[iscan] $
;                         , cam = nbcamera $
;                         , sel = nbindx, count = Nnb
;    nbstates = pertuningstates[nbindx]
;    nbfiles = pertuningfiles[nbindx]

    ;; The NB files in this scan, sorted in tuning wavelength order.
    self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                         , cam = nbcamera, scan = uscans[iscan] $
                         , sel = scan_nbindx, count = Nnb
    scan_nbfiles = pertuningfiles[scan_nbindx]
    scan_nbstates = pertuningstates[scan_nbindx]
    sortindx = sort(scan_nbstates.tun_wavelength)
    scan_nbfiles = scan_nbfiles[sortindx]
    scan_nbstates = scan_nbstates[sortindx]

    ;; The WB files in this scan, sorted as the NB files
    self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                         , cam = wbcamera, scan = uscans[iscan] $
                         , sel = scan_wbindx, count = Nwb
    scan_wbfiles = pertuningfiles[scan_wbindx]
    scan_wbstates = pertuningstates[scan_wbindx]
    match2, scan_nbstates.fpi_state, scan_wbstates.fpi_state, sortindx
    scan_wbfiles = scan_wbfiles[sortindx]
    scan_wbstates = scan_wbstates[sortindx]

    
    ;; Do WB stretch correction?
    if Nwb eq Nnb then wbstretchcorr = 1B else wbstretchcorr = 0B

    nbhdr = red_readhead(scan_nbfiles[0]) ; Use for main header
    
    ;; Make FITS header for the NB cube
    hdr = nbhdr                 ; Start with the NB cube header
    red_fitsaddkeyword, hdr, 'BITPIX', -32

    ;; Add info about this step
    prstep = 'Prepare NB science data cube'
    self -> headerinfo_addstep, hdr $
                                , prstep = prstep $
                                , prpara = prpara $
                                , prproc = inam
    
    ;; Read global WB file to use as reference when destretching
    ;; per-tuning wb files and then the corresponding nb files.
    wbim = (red_readdata(wfiles[iscan], head = wbhdr, direction = direction))[x0:x1, y0:y1]

    Nstokes = 1
    if ~keyword_set(norotation) then begin
      ff = [abs(ang),0,0,0,0]
      wbim_rot = red_rotation(wbim, ang, full = ff)
      dims = [size(wbim_rot, /dim), Ntuning, Nstokes, 1] 
    endif else dims = [size(wbim, /dim), Ntuning, Nstokes, 1] 
    
    ;; Add info to headers
    red_fitsaddkeyword, anchor = anchor, hdr, 'BUNIT', units, 'Units in array'
    red_fitsaddkeyword, anchor = anchor, hdr, 'BTYPE', 'Intensity', 'Type of data in array'

    ;; Initialize fits file, set up for writing the data part.
    self -> fitscube_initialize, filename, hdr, lun, fileassoc, dims 
    
    if prefilter eq '3950' and ~keyword_set(noaligncont) then begin
      ;; Interpolate to get the shifts for all wavelengths for
      ;; this scan.

      icont = where(scan_nbstates.prefilter eq '3999')
      xshifts = interpol([0., nb_shifts[0, iscan]] $
                         , [scan_wbstates[icont].tun_wavelength $
                            , scan_nbstates[icont].tun_wavelength]*1e7 $
                         , scan_nbstates.tun_wavelength*1e7)
      yshifts = interpol([0., nb_shifts[1, iscan]] $
                         , [scan_wbstates[icont].tun_wavelength $
                            , scan_nbstates[icont].tun_wavelength]*1e7 $
                         , scan_nbstates.tun_wavelength*1e7)
    endif

    for ituning = 0L, Ntuning - 1 do begin 

;      state = ufpi_states[ituning]


      red_progressbar, ituning, Ntuning $
                       , /predict $
                       , 'Processing scan=' + strtrim(uscans[iscan], 2) 
        
      ;; Collect info about this frame here.
      
      nbhead = red_readhead(scan_nbfiles[ituning])

      ;; DATE-??? keywords
      red_fitspar_getdates, nbhead $
                            , date_beg = date_beg $
                            , date_end = date_end $
                            , date_avg = date_avg
      date_beg_array[ituning] = date_beg
      date_end_array[ituning] = date_end
      date_avg_array[ituning] = date_avg
      tbeg_array[ituning] = red_time2double((strsplit(date_beg,'T',/extract))[1])
      tend_array[ituning] = red_time2double((strsplit(date_end,'T',/extract))[1])
      tavg_array[ituning] = red_time2double((strsplit(date_avg,'T',/extract))[1])

      ;; Wavelength and time
      wcs[ituning, 0].wave = scan_nbstates[ituning].tun_wavelength*1d9
      wcs[ituning, 0].time = tavg_array[ituning, 0]

      ;; Exposure time
      exp_array[ituning]  = fxpar(nbhead, 'XPOSURE')
      sexp_array[ituning] = fxpar(nbhead, 'TEXPOSUR')
      nsum_array[ituning] = fxpar(nbhead, 'NSUMEXP')
      
      ;; Get destretch to anchor camera (residual seeing)
      if wbstretchcorr then begin
        wwi = (red_readdata(scan_wbfiles[ituning] $
                            , direction = direction))[x0:x1, y0:y1]
        grid1 = red_dsgridnest(wbim, wwi, tile, clip)
      endif

      ;; Read image, apply prefilter curve and temporal scaling
      nbim = (red_readdata(scan_nbfiles[ituning] $
                           , direction = direction))[x0:x1, y0:y1] * rpref[ituning] 

      if prefilter eq '3950' and ~keyword_set(noaligncont) then begin
        ;; Apply alignment to compensate for time-variable chromatic
        ;; aberrations.
        nbim = red_shift_sub(nbim, -xshifts[ituning], -yshifts[ituning])
      endif

      ;; Apply destretch to anchor camera and prefilter correction
      if wbstretchcorr then nbim = red_stretch(temporary(nbim), grid1)

      self -> fitscube_addframe, fileassoc $
                                 , red_rotation(temporary(nbim), ang, full = ff) $
                                 , ituning = ituning
      
    endfor                      ; ituning

    
    ;; Get pointing at center of FOV for the different tunings.
    red_wcs_hpl_coords, tavg_array, metadata_pointing, time_pointing $
                        , hpln, hplt
 
    ;; The narrowband cube is aligned to the global wideband image
    ;; which means all narrowband scan positions are aligned to each
    ;; other. So use the median of the coordinates for the different
    ;; tunings.
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

    ;; Close fits file 
    self -> fitscube_finish, lun, wcs = wcs

    ;; Add cavity maps as WAVE distortions 
    if ~keyword_set(nocavitymap) then begin

      pindx = where(scan_nbstates.prefilter ne '3999') ; No cavity map for the Ca II H continuum
      pindx = pindx[uniq(scan_nbstates[pindx].prefilter, sort(scan_nbstates[pindx].prefilter))]
      cprefs = scan_nbstates[pindx].prefilter
      Ncprefs = n_elements(cprefs)
      
      for icprefs = 0, Ncprefs-1 do begin

        cfile = self.out_dir + 'flats/spectral_flats/' $
                + strjoin([scan_nbstates[pindx[icprefs]].detector $
                           , scan_nbstates[pindx[icprefs]].cam_settings $
                           , cprefs[icprefs] $
                           , 'fit_results.sav'] $
                          , '_')

        if ~file_test(cfile) then begin
          print, inam + ' : Error, calibration file not found -> '+cfile
          print, 'Please run the fitprefilter for '+cprefs[icprefs]+' or continue without'
          print, 'cavity map for '+cprefs[icprefs]
          stop
        endif
        restore, cfile               ; The cavity map is in a struct called "fit". 
        cmap = reform(fit.pars[1,*,*]) ; Unit is [Angstrom]
        cmap /= 10.                    ; Make it [nm]
        fit = 0B                       ; Don't need the fit struct anymore.
        
        if keyword_set(remove_smallscale) then begin
          ;; If the small scale is already corrected, then include only the
          ;; low-resolution component in the metadata. The blurring kernel
          ;; should match how the low resolution component was removed when
          ;; making flats.
          npix = 30             ; Can we get this parameter from earlier headers?
          cpsf = red_get_psf(npix*2-1,npix*2-1,double(npix),double(npix))
          cpsf /= total(cpsf, /double)
          cmap = red_convolve(temporary(cmap), cpsf)
          cmap1 = cmap
        endif else begin
          ;; If the small scale is not already corrected, then we still want
          ;; to blur the cavity map slightly.
          npsf = round(fwhm * 7.)
          if((npsf/2)*2 eq npsf) then npsf += 1L
          psf = red_get_psf(npsf, npsf, fwhm, fwhm)
          psf /= total(psf, /double)
          ;; Leave the orignal cmap alone, we might need it later.
          cmap1 = red_convolve(cmap, psf)
        endelse
        
        ;; Read the output of the pinhole calibrations so we can do the same
        ;; to the cavity maps as was done to the raw data in the momfbd
        ;; step. This output is in a struct "alignments" in the save file
        ;; 'calib/alignments.sav'
        restore,'calib/alignments.sav'
        ;; Should be based on state1 or state2 in the struct? make_cmaps
        ;; says "just pick one close to continuum (last state?)".
        indx = where(scan_nbstates[pindx[icprefs]].prefilter eq alignments.state2.prefilter, Nalign)
        case Nalign of
          0    : stop           ; Should not happen!
          1    : amap = invert(      alignments[indx].map           )
          else : amap = invert( mean(alignments[indx].map, dim = 3) )
        endcase
        cmap1 = rdx_img_project(amap, cmap1) ; Apply the geometrical mapping
        cmap1 = red_rotate(cmap1, direction)
        cmap1 = cmap1[x0:x1,y0:y1] ; Clip to the selected FOV

        stop
        ;; Write cmap(s) to the cube
        red_fitscube_addcmap, filename $
                              , reform(red_rotation(cmap1,ang,full=ff),[dims[0:1],1,1,1]) $
                              ;;, reform(cmap1, Nx, Ny, 1, 1, 1) $
                              , cmap_number = icprefs+1 $
                              , prefilter = cprefs[icprefs] $
                              , indx = where(scan_nbstates.prefilter eq cprefs[icprefs])

      endfor                    ; icprefs
      
    endif

    ;; Add some variable keywords
    self -> fitscube_addvarkeyword, filename, 'DATE-BEG', date_beg_array $
                                    , comment = 'Beginning of observation' $
                                    , keyword_value = self.isodate + 'T' + red_timestring(min(tbeg_array)) $
                                    , axis_numbers = [3] 

    self -> fitscube_addvarkeyword, filename, 'DATE-END', date_end_array $
                                    , comment = 'End time of observation' $
                                    , keyword_value = self.isodate + 'T' + red_timestring(max(tend_array)) $
                                    , axis_numbers = [3] 
    self -> fitscube_addvarkeyword, filename, 'DATE-AVG', date_avg_array $
                                    , comment = 'Average time of observation' $
                                    , keyword_value = self.isodate + 'T' + red_timestring(mean(tavg_array)) $
                                    , axis_numbers = [3] 
    
    tindx_r0 = where(time_r0 ge min(tavg_array) and time_r0 le max(tavg_array), Nt)
    if Nt gt 0 then begin
      self -> fitscube_addvarkeyword, filename, 'ATMOS_R0' $
                                      , metadata_r0[*, tindx_r0] $
                                      , comment = 'Atmospheric coherence length' $
                                      , tunit = 'm' $
                                      , extra_coordinate1 = [24, 8] $                ; WFS subfield sizes 
                                      , extra_labels      = ['WFSZ'] $               ; Axis labels for metadata_r0
                                      , extra_names       = ['WFS subfield size'] $  ; Axis names for metadata_r0
                                      , extra_units       = ['pix'] $                ; Axis units for metadata_r0
                                      , keyword_value = mean(metadata_r0[1, tindx_r0]) $
                                      , time_coordinate = time_r0[tindx_r0] $
                                      , time_unit       = 's'
    endif

    self -> fitscube_addvarkeyword, filename $
                                    , 'XPOSURE', comment = 'Summed exposure times' $
                                    , tunit = 's' $
                                    , exp_array, keyword_value = mean(exp_array) $
                                    , axis_numbers = [3] 

    self -> fitscube_addvarkeyword, filename $
                                    , 'TEXPOSUR', comment = '[s] Single-exposure time' $
                                    , tunit = 's' $
                                    , sexp_array, keyword_value = mean(sexp_array) $
                                    , axis_numbers = [3] 

    self -> fitscube_addvarkeyword, filename $
                                    , 'NSUMEXP', comment = 'Number of summed exposures' $
                                    , nsum_array, keyword_value = mean(nsum_array) $
                                    , axis_numbers = [3]

    
    ;; Include the global WB image as an image extension
    ehdr=wbhdr
    fxaddpar, ehdr, 'XTENSION', 'IMAGE'
    sxdelpar, ehdr, 'SIMPLE'
    if keyword_set(norotation) then begin
      check_fits, wbim, ehdr, /update
    endif else begin
      check_fits, wbim_rot, ehdr, /update
    endelse
    fxaddpar, ehdr, 'DATE', red_timestamp(/utc, /iso)
    anchor = 'DATE'
    red_fitsaddkeyword, anchor = anchor, ehdr, 'EXTNAME', 'WBIMAGE', 'Wideband image'
    red_fitsaddkeyword, anchor = anchor, ehdr, 'PCOUNT', 0
    red_fitsaddkeyword, anchor = anchor, ehdr, 'GCOUNT', 1
    if keyword_set(norotation) then begin
      writefits, filename, wbim, ehdr, /append
    endif else begin
      writefits, filename, wbim_rot, ehdr, /append
    endelse
    
    ;; Correct intensity with respect to solar elevation and
    ;; exposure time.
    self -> fitscube_intensitycorr, filename, corrmethod = intensitycorrmethod
    
    if keyword_set(integer) then begin
      ;; Convert to integers
      self -> fitscube_integer, filename $
                                , /delete $
                                , outname = outname $
                                , overwrite = overwrite
      filename = outname
    endif
    
    if ~keyword_set(nostatistics) then begin
      ;; Calculate statistics if not done already
      if keyword_set(norotation) then begin
        red_fitscube_statistics, filename, /write
      endif else begin
        red_fitscube_statistics, filename, /write, full = ff $
                                 , origNx = Nxx $
                                 , origNy = Nyy $
                                 , angles = [ang]
      endelse
    endif
    
    ;; Done with this scan.
    print, inam + ' : Narrowband scan cube stored in:'
    print, filename

  endfor                        ; iscan

  
end
