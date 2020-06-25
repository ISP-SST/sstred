; docformat = 'rst'

;+
; Make a de-rotated and de-stretched time-series FITS data cube with
; momfbd-restored narrow-band images.
;
; It reads all temporal de-rotation and de-stretching information from
; the wideband cube produced by companion method make_wb_cube. Cavity
; maps are added to the WCS metadata as distortions to the wavelength
; coordinate. 
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
;     wcfile : in, type=string
; 
;       The name of the corrected WB cube file.
; 
; 
; :Keywords:
; 
;     clips : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment.
;
;     cmap_fwhm : in, type=float, default=7
;   
;       FWHM in pixels of kernel used for smoothing the cavity map.
;
;     intensitycorrmethod : in, optional, type="string or boolean", default=FALSE/none
;
;       Indicate whether to do intensity correction based on WB data
;       and with what method. See documentation for red::fitscube_intensitycorr.
;
;     integer : in, optional, type=boolean
;
;       Store as integers instead of floats. Uses the BZERO and BSCALE
;       keywords to preserve the intensity scaling.
;
;     noaligncont : in, optional, type=boolean
;
;       Do not do the align continuum to wideband step.
;
;     nocavitymap : in, optional, type=boolean
;
;       Do not add cavity maps to the WCS metadata.
;
;     nostatistics : in, optional, type=boolean
;  
;       Do not calculate statistics metadata to put in header keywords
;       DATA*. 
;
;     overwrite : in, optional, type=boolean
;
;       Don't care if cube is already on disk, overwrite it
;       with a new version.
;
;     tiles : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment. 
;
;     wbsave : in, optional, type=boolean
;
;       Save a cube with the wideband per-tuning align-images. For
;       debugging of alignment with extra wideband objects.
;
; 
; :History:
; 
;    2017-08-17 : MGL. First version, based on code from
;                 chromis::make_crispex. 
;
;    2017-09-07 : MGL. Changed red_fitsaddpar --> red_fitsaddkeyword. 
; 
;    2017-09-08 : MGL. Copy variable-keywords from the WB cube. 
; 
;    2017-09-28 : MGL. Add more variable-keywords. Make a flipped cube
;                 and copy variable keywords to it. 
; 
;    2017-10-20 : MGL. Add a WCS extension.
; 
;    2017-10-27 : MGL. New keyword noaligncont.
; 
;    2017-10-30 : MGL. Incorporate code from red::make_cmaps to add
;                 cavity maps as wavelength distortions to the WCS
;                 metadata. New keyword nocavitymap. Documentation and
;                 cleanup. 
; 
;    2017-11-16 : MGL. New keyword integer. 
; 
;    2018-02-01 : MGL. New keyword wbsave. 
; 
;    2018-03-27 : MGL. Change sign on cmap. 
; 
;    2020-01-17 : MGL. New keywords intensitycorrmethod and odir. 
; 
;    2020-06-16 : MGL. Remove temporal intensity scaling, deprecate
;                 keyword notimecorr.
; 
;-
pro chromis::make_nb_cube, wcfile $
                           , clips = clips $
                           , cmap_fwhm = cmap_fwhm $
                           , integer = integer $
                           , intensitycorrmethod = intensitycorrmethod $
                           , noaligncont = noaligncont $
                           , nocavitymap = nocavitymap $
                           , noflipping = noflipping $
                           , nostatistics = nostatistics $
                           , notimecor = notimecor $
                           , odir = odir $
                           , overwrite = overwrite $
                           , tiles = tiles $
                           , wbsave = wbsave

  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  ;; Deprecated keyword:
  if n_elements(notimecor) gt 0 then begin
    print, inam + ' : Keyword notimecor is deprecated. Use intensitycorrmethod="none" instead.'
    return
  endif
  
  ;; Make prpara
  red_make_prpara, prpara, clips         
  red_make_prpara, prpara, integer
  red_make_prpara, prpara, intensitycorrmethod
  red_make_prpara, prpara, cmap_fwhm
  red_make_prpara, prpara, noaligncont 
  red_make_prpara, prpara, nocavitymap 
;  red_make_prpara, prpara, notimecor 
  red_make_prpara, prpara, np           
  red_make_prpara, prpara, overwrite
  red_make_prpara, prpara, tiles        
  red_make_prpara, prpara, wcfile

  ;; Default keywords
  if n_elements(cmap_fwhm) eq 0 then fwhm = 7.0
  if n_elements(tiles) eq 0 or n_elements(clips) eq 0 then begin
    tiles = [8, 16, 32, 64, 128]
    clips = [8, 4,  2,  1,  1  ]
  endif


  ;; Camera/detector identification
  self -> getdetectors
  wbindx     = where(strmatch(*self.cameras,'Chromis-W'))
  wbcamera   = (*self.cameras)[wbindx[0]]
  wbdetector = (*self.detectors)[wbindx[0]]
  nbindx     = where(strmatch(*self.cameras,'Chromis-N')) 
  nbcamera   = (*self.cameras)[nbindx[0]]
  nbdetector = (*self.detectors)[nbindx[0]]
  ;; Should be generalized to multiple NB cameras if CHROMIS gets
  ;; polarimetry. We don't need to identify any PD cameras for
  ;; restored data.

  ;; We do currently not correct for the small scale cavity map in
  ;; CHROMIS data. (We should get this from earlier meta data!)
  remove_smallscale = 0       
  
  ;; Read the header from the corrected WB cube. Variables begin with
  ;; WC for Wideband Cube. 
  if ~file_test(wcfile) then begin
    print, 'WB cube missing, please run make_wb_cube.'
    print, wcfile
    retall
  endif
  wchead = red_readhead(wcfile)
  ;; Read parameters from the WB cube
  fxbopen, bunit, wcfile, 'MWCINFO', bbhdr
  fxbreadm, bunit, row = 1 $
            , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01', 'DIRECTION'] $
            ,   ANG, wcCROP, wcFF, wcGRID, wcND, wcSHIFT, wcTMEAN, wcX01Y01,   direction
  ;; Note that the strarr wfiles cannot be read by fxbreadm! Put it in
  ;; wbgfiles (WideBand Global).
  fxbread, bunit, wbgfiles, 'WFILES', 1
  fxbclose, bunit

  ;; Default for wb cubes without direction parameter
  if n_elements(direction) eq 0 then direction = 0

  ;; Read wcs extension of wb file to get pointing info
  fxbopen, wlun, wcfile, 'WCS-TAB', wbdr
  ttype1 = fxpar(wbdr, 'TTYPE1')
  fxbread, wlun, wwcs, ttype1
  fxbclose, wlun

  x0 = wcX01Y01[0]
  x1 = wcX01Y01[1]
  y0 = wcX01Y01[2]
  y1 = wcX01Y01[3]
  origNx = x1 - x0 + 1
  origNy = y1 - y0 + 1

  self -> extractstates, wbgfiles, wbgstates
  prefilter = wbgstates[0].prefilter
  
  wchdr0 = red_readhead(wbgfiles[0])
  datestamp = strtrim(fxpar(wchdr0, 'STARTOBS'), 2)
  timestamp = (strsplit(datestamp, 'T', /extract))[1]
  
  search_dir = file_dirname(wbgfiles[0])+'/'
  extension = (strsplit(wbgfiles[0],'.',/extract))[-1]

  files = file_search(search_dir + '*.'+extension, count = Nfiles)      
  
  ;; Find all nb and wb per tuning files by excluding the global WB images 
  self -> selectfiles, files = files, states = states $
                       , cam = wbcamera, ustat = '' $
                       , sel = wbgindx, count = Nscans $
                       , complement = complement, Ncomplement = Ncomplement
  ;; We have no special state (or absence of state) to identify
  ;; the global WB images but we do know that their exposure times
  ;; are much larger than the ones corresponding to the individual
  ;; NB states.
  wbindx = where(states.exposure gt mean(states.exposure)*1.5 $
                  , Nscans, complement = complement, Ncomplement = Ncomplement)
  
  ;; All the per-tuning files and states
  pertuningfiles = files[complement]
  pertuningstates = states[complement]

  ;; Unique tuning states, sorted by wavelength
  utunindx = uniq(pertuningstates.fpi_state, sort(pertuningstates.fpi_state))
  Nwav = n_elements(utunindx)
  sortindx = sort(pertuningstates[utunindx].tun_wavelength)
  ufpi_states = pertuningstates[utunindx[sortindx]].fpi_state
  utunwavelength = pertuningstates[utunindx[sortindx]].tun_wavelength

  wav = utunwavelength
  my_prefilters = pertuningstates[utunindx[sortindx]].prefilter

  ;; Unique nb prefilters
  unbprefindx = uniq(pertuningstates[utunindx].prefilter, sort(pertuningstates[utunindx].prefilter))
  Nnbprefs = n_elements(unbprefindx)
  unbprefs = pertuningstates[utunindx[unbprefindx]].prefilter
  unbprefsref = dblarr(Nnbprefs)

  for inbpref = 0L, Nnbprefs-1 do begin
    ;; This is the reference point of the fine tuning for this prefilter:
    unbprefsref[inbpref] = double((strsplit(pertuningstates[utunindx[unbprefindx[inbpref]]].tuning $
                                            , '_', /extract))[0])
  endfor                        ; inbpref
  
  unbprefsref *= 1e-10          ; [m]

  ;; Get the scan selection from wfiles (from the sav file)
  self -> extractstates, wbgfiles, wbgstates
  uscans = wbgstates.scannumber
  Nscans = n_elements(uscans)

  ;; Per-tuning files, wb and nb, only for selected scans
  self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                       , scan = uscans $
                       , cam = wbcamera $
                       , sel = wbindx, count = Nwb
  wbstates = pertuningstates[wbindx]
  wbfiles = pertuningfiles[wbindx]
  self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                       , scan = uscans $
                       , cam = nbcamera $
                       , sel = nbindx, count = Nnb
  nbstates = pertuningstates[nbindx]
  nbfiles = pertuningfiles[nbindx]
  
  ;; Prepare for making output file names
  if(n_elements(odir) eq 0) then odir = self.out_dir + '/cubes_nb/' 
  midpart = prefilter + '_' + datestamp + '_scans=' $ 
            + red_collapserange(uscans, ld = '', rd = '')
  ofile = 'nb_'+midpart+'_corrected_im.fits'
  filename = odir+ofile

  ;; Already done?
  if file_test(filename) then begin
    if keyword_set(overwrite) then begin
      print, 'Overwriting existing data cube:'
      print, filename
    endif else begin
      print, 'This data cube exists already:'
      print, filename
      return
    endelse
  endif

  file_mkdir, odir

  ;; Load prefilters
  for inbpref = 0L, Nnbprefs-1 do begin
    pfile = self.out_dir + '/prefilter_fits/chromis_'+unbprefs[inbpref]+'_prefilter.idlsave'
    if ~file_test(pfile) then begin
      print, inam + ' : prefilter file not found: '+pfile
      return
    endif
      
    restore, pfile              ; Restores variable prf which is a struct
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
      red_append, prefilter_wb, prf.wbint
    endif else begin
      me = median(prf.wav)
      red_append, prefilter_curve, red_intepf(prf.wav-me, prf.pref, wav[idxpref]*1.d10-me)
      red_append, prefilter_wav, wav[idxpref]*1.d10
      red_append, prefilter_wb, replicate(prf.wbint, count)
    endelse
    
  endfor                        ; inbpref
    
  rpref = 1.d0/prefilter_curve

  ;; Do WB correction?
  if Nwb eq Nnb then wbcor = 1B else wbcor = 0B

  ;; Load WB image and define the image border
;  tmp = red_readdata(wbgfiles[0])

  ;; Spatial dimensions that match the WB cube
  Nx = wcND[0]
  Ny = wcND[1]
  
  ;; Make FITS header for the NB cube
  hdr = wchead                  ; Start with the WB cube header

;  if keyword_set(integer) then begin
;
;    ;; We need to find out BZERO and BSCALE. For now we have to read
;    ;; through all data in the cube. If the momfbd output files could
;    ;; store min and max values in the future, this would be quicker.
;    iprogress = 0
;    Nprogress = Nscans*Nwav
;
;    for iscan = 0L, Nscans-1 do begin
;
;      self -> selectfiles, files = pertuningfiles, states = pertuningstates $
;                           , cam = nbcamera, scan = uscans[iscan] $
;                           , sel = scan_nbindx, count = count
;      scan_nbfiles = pertuningfiles[scan_nbindx]
;
;      if keyword_set(notimecor) then tscl = 1. else tscl = mean(prefilter_wb) / wcTMEAN[iscan]
;
;      for iwav = 0L, Nwav - 1 do begin
;
;        red_progressbar, iprogress, Nprogress $
;                         , /predict $
;                         , 'Calculating BZERO and BSCALE' 
;  
;        mr = momfbd_read(scan_nbfiles[iwav], /img)
;
;        ;; Min and max after scaling to units in file
;        datamax_thisfile = max(mr.patch.img * rpref[iwav] * tscl, min = datamin_thisfile)
;
;        ;; Global min and max
;        if iscan eq 0 and iwav eq 0 then begin
;          datamin = datamin_thisfile
;          datamax = datamax_thisfile
;        endif else begin
;          datamin <= datamin_thisfile
;          datamax >= datamax_thisfile
;        endelse
;
;        iprogress++
;        
;      endfor                    ; iwav
;    endfor                      ; iscan
;
;    ;; Approximate min and max values after scaling to 2-byte
;    ;; integers. Don't use -32768 and 32767, want to leave some margin
;    ;; for altered physical values due to resampling (most likely to
;    ;; lower abs values, but still...).
;    arraymin = -32000.
;    arraymax =  32000.
;
;    ;; BZERO and BSCALE
;    bscale = (datamax-datamin)/(arraymax-arraymin)
;    bzero = datamin - arraymin*bscale
;    
;    ;; Set keywords for rescaled integers
;    red_fitsaddkeyword, hdr, 'BITPIX', 16
;    red_fitsaddkeyword, hdr, 'BZERO', bzero
;    red_fitsaddkeyword, hdr, 'BSCALE', bscale
;    
;  endif else begin
  ;; If not rescaled integers, we just need to set BITPIX for floats
  red_fitsaddkeyword, hdr, 'BITPIX', -32
;  endelse
  
  ;; Add info about this step
  self -> headerinfo_addstep, hdr $
                              , prstep = 'CONCATENATION' $
                              , prpara = prpara $
                              , prproc = inam

  self -> headerinfo_addstep, hdr $
                              , prstep = 'CALIBRATION-INTENSITY-SPECTRAL' $
                              , prpara = prpara $
                              , prref = ['Hamburg FTS spectral atlas (Neckel 1999)' $
                                         , 'Calibration data from '+red_timestring(prf.time_avg, n = 0)] $
                              , prproc = inam

  anchor = 'DATE'

  ;; Add some keywords
  red_fitsaddkeyword, anchor = anchor, hdr, 'OBS_HDU', 1
  
  ;; Add info to headers
  red_fitsaddkeyword, anchor = anchor, hdr, 'BUNIT', units, 'Units in array'
  red_fitsaddkeyword, anchor = anchor, hdr, 'BTYPE', 'Intensity', 'Type of data in array'

  ;; WB and NB data come from different cameras.
  red_fitsaddkeyword, hdr, 'CAMERA',   nbstates[0].camera
  ;; Get DETGAIN, DETOFFS, DETMODEL, DETFIRM from .fitsheader file,
  ;; i.e., red_readhead(.momfbd file). This would have to be handled
  ;; differently with CRISP because each Stokes image is a mix of two
  ;; detectors. 
  mhd = red_readhead(nbfiles[0]) ; Header of momfbd output file
  red_fitscopykeyword, hdr, 'DETECTOR', mhd
  red_fitscopykeyword, hdr, 'DETGAIN',  mhd
  red_fitscopykeyword, hdr, 'DETOFFS',  mhd
  red_fitscopykeyword, hdr, 'DETMODEL', mhd
  red_fitscopykeyword, hdr, 'DETFIRM',  mhd

  ;; Initialize fits file, set up for writing the data part.
  dims = [Nx, Ny, Nwav, 1, Nscans] 
  self -> fitscube_initialize, filename, hdr, lun, fileassoc, dims 

  if keyword_set(wbsave) then begin
    wbfilename = red_strreplace(filename, 'nb_', 'wbalign_')
    self -> fitscube_initialize, wbfilename, hdr, wblun, wbfileassoc, dims 
  endif
  

  ;; Observations metadata varaibles
  tbeg_array     = dblarr(Nwav, Nscans) ; Time beginning for state
  tend_array     = dblarr(Nwav, Nscans) ; Time end for state
  tavg_array     = dblarr(Nwav, Nscans) ; Time average for state
  date_beg_array = strarr(Nwav, Nscans) ; DATE-BEG for state
  date_end_array = strarr(Nwav, Nscans) ; DATE-END for state
  date_avg_array = strarr(Nwav, Nscans) ; DATE-AVG for state
  exp_array      = fltarr(Nwav, Nscans) ; Total exposure time
  sexp_array     = fltarr(Nwav, Nscans) ; Single exposure time
  nsum_array     = lonarr(Nwav, Nscans) ; Number of summed exposures
 
  wcs = replicate({  wave:dblarr(2,2) $
                   , hplt:dblarr(2,2) $
                   , hpln:dblarr(2,2) $
                   , time:dblarr(2,2) $
                  }, Nwav, Nscans)

;  if ~keyword_set(nostatistics) then begin
;    ;; Variables for statistics metadata.
;    DATAMIN  = dblarr(Nwav, Nscans) ; The minimum data value
;    DATAMAX  = dblarr(Nwav, Nscans) ; The maximum data value
;    DATAMEAN = dblarr(Nwav, Nscans) ; The average data value
;    DATARMS  = dblarr(Nwav, Nscans) ; The RMS deviation from the mean
;    DATAKURT = dblarr(Nwav, Nscans) ; The kurtosis
;    DATASKEW = dblarr(Nwav, Nscans) ; The skewness
;    DATAP01  = dblarr(Nwav, Nscans) ; The 01 percentile
;    DATAP10  = dblarr(Nwav, Nscans) ; The 10 percentile
;    DATAP25  = dblarr(Nwav, Nscans) ; The 25 percentile
;    DATAMEDN = dblarr(Nwav, Nscans) ; The median data value
;    DATAP75  = dblarr(Nwav, Nscans) ; The 75 percentile
;    DATAP90  = dblarr(Nwav, Nscans) ; The 90 percentile
;    DATAP95  = dblarr(Nwav, Nscans) ; The 95 percentile
;    DATAP98  = dblarr(Nwav, Nscans) ; The 98 percentile
;    DATAP99  = dblarr(Nwav, Nscans) ; The 99 percentile
;    percentiles = [.01, .10, .25, .50, .75, .90, .95, .98, .99]
;  end

    
  ;; The narrowband cube is aligned to the wideband cube and all
  ;; narrowband scan positions are aligned to each other. So get hpln
  ;; and hplt from the wideband cube wcs coordinates, this should
  ;; apply to all frames in a scan.
  for iscan = 0L, Nscans-1 do begin
    for iwav = 0, Nwav-1 do begin
      ;; We rely here on hpln and hplt being the first two tabulated
      ;; coordinates. To make this more general, we should get the
      ;; actual indices from the headers. Maybe later...
      wcs[iwav, iscan].hpln = reform(wwcs[0,*,*,iscan])
      wcs[iwav, iscan].hplt = reform(wwcs[1,*,*,iscan])
    endfor                      ; iwav
  endfor                        ; iscan
  
  ;; Continuum alignment only done for Ca II scans (so far). H beta is
  ;; not as wide so should be OK.
  if prefilter eq '3950' and ~keyword_set(noaligncont) then begin
    
    ;; Get wavelength-variable shifts based on continuum vs wideband
    ;; alignment.
    
    aligndir = self.out_dir + '/align/' + timestamp $
               + '/' + prefilter + '/'
    
    nname = aligndir+'scan_numbers.fz'
    sname = aligndir+'continuum_shifts_smoothed.fz'
    
    if ~file_test(nname) or ~file_test(sname) then begin
      print, inam + ' : At least one file missing for aligncont option:'
      print, nname
      print, sname
      retall
    endif
    
    ;; Read the shifts for the continuum images
    fzread, align_scannumbers, nname
    fzread, align_shifts, sname, align_header

    ;; Get the wavelengths used for the intra-scan alignment from the
    ;; file header.
    if n_elements(align_header) gt 0 then align_wavelengths = double(strsplit(align_header,/extract))
    if n_elements(align_wavelengths) ne 2 then begin
      ;; Default to WB and NB continuum as used by the align_continuum
      ;; method.
      align_wavelengths = [3.950e-07, 4.000e-07]
    endif
    
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

  
;  if ~keyword_set(nostatistics) then begin
;    ;; Read all images once in order to get statistics. Min and max need
;    ;; to be calculated here so we can build the histogram in the next
;    ;; loop.
;    iprogress = 0
;    Nprogress = Nscans*Nwav
;    for iscan = 0L, Nscans - 1 do begin
;      for iwav = 0L, Nwav - 1 do begin 
;
;        red_progressbar, iprogress, Nprogress $
;                         , /predict $
;                         , 'Calculate statistics'
;
;
;        ;; The NB files in this scan, sorted in tuning wavelength order.
;        self -> selectfiles, files = pertuningfiles, states = pertuningstates $
;                             , cam = nbcamera, scan = uscans[iscan] $
;                             , sel = scan_nbindx, count = count
;        scan_nbfiles = pertuningfiles[scan_nbindx]
;        scan_nbstates = pertuningstates[scan_nbindx]
;
;        if keyword_set(notimecor) then tscl = 1. else tscl = mean(prefilter_wb) / wcTMEAN[iscan]
;        
;        ;; Read image, apply prefilter curve and temporal scaling
;        nbim = (red_readdata(scan_nbfiles[iwav]))[x0:x1, y0:y1] * rpref[iwav] * tscl
;
;        ;; Calculate statistics metadata before alignment and rotation
;        ;; so we avoid padding.
;        momnt = moment(nbim)    ; mean, variance, skewness, kurtosis
;        perc = cgpercentiles(nbim, percentiles = percentiles)
;        
;        DATAMIN[iwav, iscan]  = min(nbim)
;        DATAMAX[iwav, iscan]  = max(nbim)
;
;        DATAMEAN[iwav, iscan] = momnt[0]         
;        DATARMS[iwav, iscan]  = sqrt(momnt[1])   
;        DATASKEW[iwav, iscan] = momnt[2]         
;        DATAKURT[iwav, iscan] = momnt[3]         
;
;        DATAP01[iwav, iscan]  = perc[0]          
;        DATAP10[iwav, iscan]  = perc[1]          
;        DATAP25[iwav, iscan]  = perc[2]          
;        DATAMEDN[iwav, iscan] = perc[3]          
;        DATAP75[iwav, iscan]  = perc[4]          
;        DATAP90[iwav, iscan]  = perc[5]          
;        DATAP95[iwav, iscan]  = perc[6]          
;        DATAP98[iwav, iscan]  = perc[7]          
;        DATAP99[iwav, iscan]  = perc[8]          
;
;        iprogress++
;        
;      endfor                    ; iwav
;    endfor                      ; iscan
;    CUBEMIN  = min(DATAMIN)
;    CUBEMAX  = max(DATAMAX)
;
;    ;; Accumulate a histogram for the entire cube, use to calculate
;    ;; percentiles. Setup here, add contributions in the double loop
;    ;; below. 
;    Nbins = 2L^16               ; Use many bins!
;    binsize = (CUBEMAX - CUBEMIN) / (Nbins - 1.)
;    hist = lonarr(Nbins)
;  endif
  
  iprogress = 0
  Nprogress = Nscans*Nwav
  for iscan = 0L, Nscans-1 do begin

    ;; The NB files in this scan, sorted in tuning wavelength order.
    self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                         , cam = nbcamera, scan = uscans[iscan] $
                         , sel = scan_nbindx, count = count
    scan_nbfiles = pertuningfiles[scan_nbindx]
    scan_nbstates = pertuningstates[scan_nbindx]
    sortindx = sort(scan_nbstates.tun_wavelength)
    scan_nbfiles = scan_nbfiles[sortindx]
    scan_nbstates = scan_nbstates[sortindx]

    ;; The WB files in this scan, sorted as the NB files
    self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                         , cam = wbcamera, scan = uscans[iscan] $
                         , sel = scan_wbindx, count = count
    scan_wbfiles = pertuningfiles[scan_wbindx]
    scan_wbstates = pertuningstates[scan_wbindx]
    match2, scan_nbstates.fpi_state, scan_wbstates.fpi_state, sortindx
    scan_wbfiles = scan_wbfiles[sortindx]
    scan_wbstates = scan_wbstates[sortindx]
    
    ;; Read global WB file to use as reference when destretching
    ;; per-tuning wb files and then the corresponding nb files.
    wb = (red_readdata(wbgfiles[iscan], direction = direction))[x0:x1, y0:y1]
    
    if prefilter eq '3950' and ~keyword_set(noaligncont) then begin
      ;; Interpolate to get the shifts for all wavelengths for
      ;; this scan.
      icont = where(scan_nbstates.prefilter eq '3999')
      xshifts = interpol([0., nb_shifts[0, iscan]] $
                         , align_wavelengths*1e7 $
                         , scan_nbstates.tun_wavelength*1e7)
      yshifts = interpol([0., nb_shifts[1, iscan]] $
                         , align_wavelengths*1e7 $
                         , scan_nbstates.tun_wavelength*1e7)
    endif

    ;;stop
    ;;if keyword_set(notimecor) then tscl = 1. else tscl = mean(wcTMEAN) / wcTMEAN[iscan]
    ;;tscl *= mean(prefilter_wb)
;    tscl = 1.
    
    for iwav = 0L, Nwav - 1 do begin 

      state = ufpi_states[iwav]

      ;; Collect info about this frame here.
      
      nbhead = red_readhead(scan_nbfiles[iwav])

      ;; DATE-??? keywords
      red_fitspar_getdates, nbhead $
                            , date_beg = date_beg $
                            , date_end = date_end $
                            , date_avg = date_avg
      date_beg_array[iwav, iscan] = date_beg
      date_end_array[iwav, iscan] = date_end
      date_avg_array[iwav, iscan] = date_avg
      tbeg_array[iwav, iscan] = red_time2double((strsplit(date_beg,'T',/extract))[1])
      tend_array[iwav, iscan] = red_time2double((strsplit(date_end,'T',/extract))[1])
      tavg_array[iwav, iscan] = red_time2double((strsplit(date_avg,'T',/extract))[1])

      ;; Wavelength and time
      wcs[iwav, iscan].wave = scan_nbstates[iwav].tun_wavelength*1d9
      wcs[iwav, iscan].time = tavg_array[iwav, iscan]

      ;; Exposure time
      exp_array[iwav, iscan]  = fxpar(nbhead, 'XPOSURE')
      sexp_array[iwav, iscan] = fxpar(nbhead, 'TEXPOSUR')
      nsum_array[iwav, iscan] = fxpar(nbhead, 'NSUMEXP')
      
      red_progressbar, iprogress, Nprogress $
                       , /predict $
                       , 'Processing scan=' $
                       + strtrim(uscans[iscan], 2) + ' state=' + state 

      ;; Get destretch to anchor camera (residual seeing)
      if wbcor then begin
        wwi = (red_readdata(scan_wbfiles[iwav], direction = direction))[x0:x1, y0:y1]
        grid1 = red_dsgridnest(wb, wwi, tiles, clips)
      endif

      ;; Read image, apply prefilter curve and temporal scaling
      nbim = (red_readdata(scan_nbfiles[iwav], direction = direction))[x0:x1, y0:y1] * rpref[iwav] ;* tscl

;      if ~keyword_set(nostatistics) then begin
;        ;; Add to the histogram
;        hist += histogram(nbim, min = cubemin, max = cubemax, Nbins = Nbins, /nan)
;      endif
      
      if prefilter eq '3950' and ~keyword_set(noaligncont) then begin
        ;; Apply alignment to compensate for time-variable chromatic
        ;; aberrations.
        nbim = red_shift_sub(nbim, -xshifts[iwav], -yshifts[iwav])
      endif

      ;; Apply destretch to anchor camera and prefilter correction
      if wbcor then nbim = red_stretch(temporary(nbim), grid1)

      ;; Apply derot, align, dewarp based on the output from
      ;; make_wb_cube
      nbim = red_rotation(temporary(nbim), ang[iscan], $
                          wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF)
      nbim = red_stretch(temporary(nbim), reform(wcGRID[iscan,*,*,*]))

;      if keyword_set(integer) then begin
;        self -> fitscube_addframe, fileassoc, fix(round((temporary(nbim)-bzero)/bscale)) $
;                                   , iscan = iscan, ituning = iwav
;      endif else begin
      self -> fitscube_addframe, fileassoc, temporary(nbim) $
                                 , iscan = iscan, ituning = iwav
;      endelse

      if keyword_set(wbsave) then begin
        ;; Same operations as on narrowband image, except for
        ;; "aligncont".
        wbim = wwi              ;* tscl
        wbim = red_stretch(temporary(wbim), grid1)
        wbim = red_rotation(temporary(wbim), ang[iscan], $
                          wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF)
        wbim = red_stretch(temporary(wbim), reform(wcGRID[iscan,*,*,*]))
        self -> fitscube_addframe, wbfileassoc, temporary(wbim) $
                                   , iscan = iscan, ituning = iwav
      endif
      
      iprogress++               ; update progress counter

    endfor                      ; iwav
    
  endfor                        ; iscan

  ;; Close fits file.
  self -> fitscube_finish, lun, wcs = wcs
  if keyword_set(wbsave) then self -> fitscube_finish, wblun, wcs = wcs


  ;; Create cubes for science data and scan-adapted cavity maps.
  cavitymaps = fltarr(Nx, Ny, 1, 1, Nscans)

  if ~keyword_set(nocavitymap) then begin

    ;; Read the original cavity map
    pindx = where(nbstates.prefilter ne '3999') ; No cavity map for the Ca II H continuum
    pindx = pindx[uniq(nbstates[pindx].prefilter, sort(nbstates[pindx].prefilter))]
    cprefs = nbstates[pindx].prefilter
    Ncprefs = n_elements(cprefs)
    
    for icprefs = 0, Ncprefs-1 do begin
      cfile = self.out_dir + 'flats/spectral_flats/' $
              + strjoin([nbstates[pindx[icprefs]].detector $
                         , nbstates[pindx[icprefs]].cam_settings $
                         , cprefs[icprefs] $
                         , 'fit_results.sav'] $
                        , '_')

      if ~file_test(cfile) then begin
        print, inam + ' : Error, calibration file not found -> '+cfile
        print, 'Please run the fitprefilter for '+cprefs[icprefs]+' or continue without'
        print, 'cavity map for '+cprefs[icprefs]
        stop
      endif
      restore, cfile                 ; The cavity map is in a struct called "fit". 
      cmap = reform(fit.pars[1,*,*]) ; Unit is [Angstrom]
      cmap = rotate(temporary(cmap), direction)
      cmap /= 10.               ; Make it [nm]
      cmap = -cmap                   ; Change sign so lambda_correct = lambda + cmap
      fit = 0B                       ; Don't need the fit struct anymore.
      
      if keyword_set(remove_smallscale) then begin
        ;; If the small scale is already corrected, then include only the
        ;; low-resolution component in the metadata. The blurring kernel
        ;; should match how the low resolution component was removed when
        ;; making flats.
        npix = 30               ; Can we get this parameter from earlier headers?
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
      indx = where(nbstates[0].prefilter eq alignments.state2.prefilter, Nalign)
      case Nalign of
        0    : stop             ; Should not happen!
        1    : amap = invert(      alignments[indx].map           )
        else : amap = invert( mean(alignments[indx].map, dim = 3) )
      endcase
      cmap1 = rdx_img_project(amap, cmap1) ; Apply the geometrical mapping
      cmap1 = cmap1[x0:x1,y0:y1]           ; Clip to the selected FOV

      ;; Now make rotated copies of the cavity map
      for iscan = 0L, Nscans-1 do begin

        if ~keyword_set(nocavitymap) then begin
          
          ;; Apply the same derot, align, dewarp as for the science data
          cmap11 = red_rotation(cmap1, ang[iscan], $
                                wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF)
          cmap11 = red_stretch(temporary(cmap11), reform(wcGRID[iscan,*,*,*]))
          
          cavitymaps[0, 0, 0, 0, iscan] = cmap11

          ;; The following block of code is inactive but we want to keep
          ;; it around in case it is needed later. It does further
          ;; operations on the cavity maps based on what is done to the
          ;; science data, such as stretching based on the extra
          ;; per-tuning wideband objects, as well as blurring based on the
          ;; momfbd psfs. It should probably have a boolean keyword that
          ;; activates it. (This code will have to be updated to the
          ;; current pipeline style before it can be used.)
          if 0 then begin
            wb = (red_readdata(wbf[ss]))[x0:x1,y0:y1]
            for ww = 0L, nw - 1 do begin
              
              iwbf = strjoin([self.camwbtag, (strsplit(file_basename(st.ofiles[ww,ss]) $
                                                       , '.', /extract))[1:*]],'.')
              iwbf = file_dirname(st.ofiles[ww,ss]) + '/'+iwbf
              
              ;; Load images
              iwb = (red_readdata(iwbf))[x0:x1,y0:y1]
              im = momfbd_read(st.ofiles[ww,ss])
              
              ;; get dewarp from WB
              igrid = red_dsgridnest(wb, iwb, itiles, iclip)
              
              ;; Convolve CMAP and apply wavelength dep. de-warp
              cmap2 = red_stretch((red_mozaic(red_conv_cmap(cmap, im), /crop))[x0:x1, y0:y1], igrid)
              
              ;; Derotate and shift
              cmap2 = red_rotation(temporary(cmap2), ang[ss], total(shift[0,ss]), $
                                   total(shift[1,ss]), full=wcFF)
              
              ;; Time de-warp
              cmap2 = red_stretch(temporary(cmap2), reform(grid[ss,*,*,*]))
              
            endfor              ; ww
          endif
        endif

      endfor                    ; iscan

      tindx = where(scan_nbstates.prefilter eq cprefs[icprefs], Nt)
      
      ;; Add cavity maps as WAVE distortions
      if Nt gt 0 then begin
        red_fitscube_addcmap, filename, cavitymaps $
                              , cmap_number = icprefs+1 $
                              , prefilter = cprefs[icprefs] $
                              , indx = tindx
      endif
      
    endfor                      ; icprefs
    
  endif 
  
  
  ;; Add some variable keywords
  self -> fitscube_addvarkeyword, filename, 'DATE-BEG', date_beg_array $
                                  , anchor = anchor $
                                  , comment = 'Beginning time of observation' $
                                  , keyword_method = 'first' $
;                                  , keyword_value = self.isodate + 'T' + red_timestring(min(tbeg_array)) $
                                  , axis_numbers = [3, 5] 
  self -> fitscube_addvarkeyword, filename, 'DATE-END', date_end_array $
                                  , anchor = anchor $
                                  , comment = 'End time of observation' $
                                  , keyword_method = 'last' $
;                                  , keyword_value = self.isodate + 'T' + red_timestring(max(tend_array)) $
                                  , axis_numbers = [3, 5] 
  self -> fitscube_addvarkeyword, filename, 'DATE-AVG', date_avg_array $
                                  , anchor = anchor $
                                  , comment = 'Average time of observation' $
                                  , keyword_value = self.isodate + 'T' + red_timestring(mean(tavg_array)) $
                                  , axis_numbers = [3, 5] 
 
  ;; Copy variable-keywords from wb cube file.
  self -> fitscube_addvarkeyword, filename, 'SCANNUM',  old_filename = wcfile $
                                  , anchor = anchor 
  self -> fitscube_addvarkeyword, filename, 'ATMOS_R0', old_filename = wcfile $
                                  , anchor = anchor 
  self -> fitscube_addvarkeyword, filename, 'AO_LOCK', old_filename = wcfile $
                                  , anchor = anchor 
  self -> fitscube_addvarkeyword, filename, 'ELEV_ANG', old_filename = wcfile $
                                  , anchor = anchor 

  self -> fitscube_addvarkeyword, filename, 'XPOSURE', exp_array $
                                  , comment = 'Summed exposure times' $
                                  , anchor = anchor $
                                  , tunit = 's' $
                                  , keyword_method = 'median' $
;                                  , keyword_value = mean(exp_array) $
                                  , axis_numbers = [3, 5] 

  self -> fitscube_addvarkeyword, filename, 'TEXPOSUR', sexp_array $
                                  , comment = '[s] Single-exposure time' $
                                  , anchor = anchor $
                                  , tunit = 's' $
                                  , keyword_method = 'median' $
;                                  , keyword_value = mean(sexp_array) $
                                  , axis_numbers = [3, 5] 

  self -> fitscube_addvarkeyword, filename, 'NSUMEXP', nsum_array $
                                  , comment = 'Number of summed exposures' $
                                  , anchor = anchor $
                                  , keyword_method = 'median' $
;                                  , keyword_value = mean(nsum_array) $
                                  , axis_numbers = [3, 5]

  
  ;; Correct intensity with respect to solar elevation and exposure
  ;; time.
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
    ;; Calculate statistics 
;    red_fitscube_statistics, filename, /write 
    red_fitscube_statistics, filename, /write $
                             , angles = ang $
                             , full = wcFF $
                             , grid = wcGRID $
                             , origNx = origNx $
                             , origNy = origNy $
;                             , percentiles = percentiles $
                             , shifts = wcSHIFT 
  endif

  if ~keyword_set(noflipping) then $
     red_fitscube_flip, filename, flipfile = flipfile $
                        , overwrite = overwrite

  print, inam + ' : Narrowband cube stored in:'
  print, filename
  if ~keyword_set(noflipping) then print, flipfile
  
  if keyword_set(wbsave) then begin
    if ~keyword_set(noflipping) then red_fitscube_flip, wbfilename, flipfile = wbflipfile
    print, inam + ' : Wideband align cube stored in:'
    print, wbfilename
    if ~keyword_set(noflipping) then print, wbflipfile
  endif
  
end
