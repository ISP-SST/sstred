; docformat = 'rst'

;+
; Make FITS data cubes with momfbd-restored images, one scan per file.
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
;     clips : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment.
;
;     cmap_fwhm : in, type=float, default=7
;   
;       FWHM in pixels of kernel used for smoothing the cavity map.
;
;     integer : in, optional, type=boolean
;
;       Store as integers instead of floats.
;
;     noaligncont : in, optional, type=boolean
;
;       Do not do the align continuum to wideband step.
;
;     nocavitymap : in, optional, type=boolean
;
;       Do not add cavity maps to the WCS metadata.
;
;     notimecor : in, optional, type=boolean
;
;       Skip temporal correction of intensities.
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
; 
; :History:
; 
;    2018-01-19 : MGL. First version, based on code from
;                 chromis::make_wb_cube and chromis::make_nb_cube.
; 
;-
pro chromis::make_scan_cube, dir $
                             , autocrop = autocrop $
                             , clip = clip $
                             , crop = crop $
                             , integer = integer $
                             , interactive = interactive $
                             , noaligncont = noaligncont $
                             , nocavitymap = nocavitymap $
                             , overwrite = overwrite $
                             , tile = tile $
                             , scannos = scannos 
                            
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  ;; Temporarily disable cavity maps by default, can still be be
  ;; written (experimentally) with explicit nocavitymap=0.
  if n_elements(nocavitymap) eq 0 then nocavitymap = 1
  
  ;; Make prpara
  red_make_prpara, prpara, dir
  red_make_prpara, prpara, clip
  red_make_prpara, prpara, crop     
  red_make_prpara, prpara, integer
  red_make_prpara, prpara, noaligncont 
  red_make_prpara, prpara, nocavitymap    
  red_make_prpara, prpara, overwrite
  red_make_prpara, prpara, tile

  ;; Default keywords
  if n_elements(cmap_fwhm) eq 0 then fwhm = 7.0
  if n_elements(clip) eq 0 then clip = [8, 4,  2,  1,  1  ]
  if n_elements(tile) eq 0 then tile = [8, 16, 32, 64, 128]

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
  red_logdata, self.isodate, time_pig, pig = metadata_pig, rsun = rsun


  ;; Search for restored images
  case self.filetype of
    'ANA': extension = '.f0'
    'MOMFBD': extension = '.momfbd'
    'FITS': extension = '.fits'
  endcase
  files = file_search(dir + '*'+extension, count = Nfiles)      
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
  datestamp = fxpar(red_readhead(wfiles[0]), 'STARTOBS')
  timestamp = (strsplit(datestamp, 'T', /extract))[1]

  
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
  hdr = red_readhead(wfiles[0])
  im_dim = fxpar(hdr, 'NAXIS*')
  x0 = crop[0]
  x1 = im_dim[0]-1 - crop[1]
  y0 = crop[2]
  y1 = im_dim[1]-1 - crop[3]
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

  
  if ~keyword_set(nocavitymap) then begin

    cavitymap = fltarr(Nx, Ny, 1)

    ;; Read the original cavity map
    istate = 0                  ; Just pick one
    cfile = self.out_dir + 'flats/spectral_flats/' $
            + strjoin([nbstates[istate].detector $
                       ,nbstates[istate].cam_settings $
                       ,nbstates[istate].prefilter $
                       , 'fit_results.sav'] $
                      , '_')

    if ~file_test(cfile) then begin
      print, inam + ' : Error, calibration file not found -> '+cfile
      stop
    endif
    restore, cfile                 ; The cavity map is in a struct called "fit". 
    cmap = reform(fit.pars[1,*,*]) ; Unit is [Angstrom]
    cmap /= 10.                    ; Make it [nm]
    fit = 0B                       ; Don't need the fit struct anymore.
    
    if keyword_set(remove_smallscale) then begin
      ;; If the small scale is already corrected, then include only the
      ;; low-resolution component in the metadata. The blurring kernel
      ;; should match how the low resolution component was removed when
      ;; making flats.
      npix = 30                 ; Can we get this parameter from earlier headers?
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
      0    : stop               ; Should not happen!
      1    : amap = invert(      alignments[indx].map           )
      else : amap = invert( mean(alignments[indx].map, dim = 3) )
    endcase
    cmap1 = rdx_img_project(amap, cmap1) ; Apply the geometrical mapping
    cmap1 = cmap1[x0:x1,y0:y1]           ; Clip to the selected FOV
   
  endif
  
    
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

    ;; Use interpolation to get the shifts for the selected scans.
    nb_shifts = fltarr(2, Nscans)
    for iscan=0L, Nscans-1 do begin
      pos = where(align_scannumbers eq uscans[iscan], cccc)
      if cccc eq 1 then nb_shifts[*, iscan] = align_shifts[*, pos] else begin
        nb_shifts[0, *] = interpol([reform(align_shifts[0, *])] $
                                   , [float(align_scannumbers)], [float(uscans)])
        nb_shifts[1, *] = interpol([reform(align_shifts[1, *])] $
                                   , [float(align_scannumbers)], [float(uscans)])
      endelse
    endfor
    pos = where(~finite(nb_shifts), cccc)
    if cccc gt 0 then nb_shifts[pos] = 0
  endif


  for iscan = 0L, Nscans-1 do begin
    
    ;; This is the loop in which the cubes are written.
    
    ;; Make output file name
    midpart = prefilter + '_' + datestamp + '_scan=' $ 
              + strtrim(uscans[iscan], 2)
    if keyword_set(integer) then begin
      ofile = 'nb_'+midpart+'_corrected_int.fits'
    endif else begin
      ofile = 'nb_'+midpart+'_corrected.fits'
    endelse
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
    tbeg_array     = dblarr(Nwav)   ; Time beginning for state
    tend_array     = dblarr(Nwav)   ; Time end for state
    tavg_array     = dblarr(Nwav)   ; Time average for state
    date_beg_array = strarr(Nwav)   ; DATE-BEG for state
    date_end_array = strarr(Nwav)   ; DATE-END for state
    date_avg_array = strarr(Nwav)   ; DATE-AVG for state
    exp_array      = fltarr(Nwav)   ; Total exposure time
    sexp_array     = fltarr(Nwav)   ; Single exposure time
    nsum_array     = lonarr(Nwav)   ; Number of summed exposures

    wcs = replicate({  wave:dblarr(2,2) $
                       , hplt:dblarr(2,2) $
                       , hpln:dblarr(2,2) $
                       , time:dblarr(2,2) $
                    }, Nwav)

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

    
    ;; Do WB correction?
    if Nwb eq Nnb then wbcor = 1B else wbcor = 0B

    nbhdr = red_readhead(scan_nbfiles[0]) ; Use for main header
    wbhdr = red_readhead(scan_wbfiles[0]) ; Use for header of wb data extension?
    
    ;; Make FITS header for the NB cube
    hdr = nbhdr                 ; Start with the NB cube header
    if keyword_set(integer) then begin

      ;; We need to find out BZERO and BSCALE. For now we have to read
      ;; through all data in the cube. If the momfbd output files could
      ;; store min and max values in the future, this would be quicker.

      for iwav = 0L, Nwav - 1 do begin

;        nbim = (red_readdata(scan_nbfiles[iwav]))[x0:x1, y0:y1] * rpref[iwav] 
;          red_progressbar, iprogress, Nprogress $
;                           , /predict $
;                           , 'Calculating BZERO and BSCALE' 
        
        mr = momfbd_read(scan_nbfiles[iwav], /img)
        ;; Min and max after scaling to units in file
        datamax_thisfile = max(mr.patch.img * rpref[iwav], min = datamin_thisfile)

        ;; Global min and max
        if iwav eq 0 then begin
          datamin = datamin_thisfile
          datamax = datamax_thisfile
        endif else begin
          datamin <= datamin_thisfile
          datamax >= datamax_thisfile
        endelse

      endfor                    ; iwav
        
      ;; Approximate min and max values after scaling to 2-byte
      ;; integers. Don't use -32768 and 32767, want to leave some margin
      ;; for altered physical values due to resampling (most likely to
      ;; lower abs values, but still...).
      arraymin = -32000.
      arraymax =  32000.

      ;; BZERO and BSCALE
      bscale = (datamax-datamin)/(arraymax-arraymin)
      bzero = datamin - arraymin*bscale
      
      ;; Set keywords for rescaled integers
      red_fitsaddkeyword, hdr, 'BITPIX', 16
      red_fitsaddkeyword, hdr, 'BZERO', bzero
      red_fitsaddkeyword, hdr, 'BSCALE', bscale
      
    endif else begin
      ;; If not rescaled integers, we just need to set BITPIX for floats
      red_fitsaddkeyword, hdr, 'BITPIX', -32
    endelse
    
    ;; Add info about this step
    self -> headerinfo_addstep, hdr $
                                , prstep = 'Prepare NB science data cube' $
                                , prpara = prpara $
                                , prproc = inam
    
    ;; Add info to headers
    red_fitsaddkeyword, anchor = anchor, hdr, 'BUNIT', units, 'Units in array'
    red_fitsaddkeyword, anchor = anchor, hdr, 'BTYPE', 'Intensity', 'Type of data in array'

    ;; Initialize fits file, set up for writing the data part.
    dims = [Nx, Ny, Nwav, 1, 1] 
    self -> fitscube_initialize, filename, hdr, lun, fileassoc, dims 

    
    ;; Read global WB file to use as reference when destretching
    ;; per-tuning wb files and then the corresponding nb files.
    wb = (red_readdata(wfiles[iscan]))[x0:x1, y0:y1]
    
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

    for iwav = 0L, Nwav - 1 do begin 

;      state = ufpi_states[iwav]


      red_progressbar, iwav, Nwav $
                       , /predict $
                       , 'Processing scan=' + strtrim(uscans[iscan], 2) 
        
              ;; Collect info about this frame here.
      
      nbhead = red_readhead(scan_nbfiles[iwav])

      ;; DATE-??? keywords
      red_fitspar_getdates, nbhead $
                            , date_beg = date_beg $
                            , date_end = date_end $
                            , date_avg = date_avg
      date_beg_array[iwav] = date_beg
      date_end_array[iwav] = date_end
      date_avg_array[iwav] = date_avg
      tbeg_array[iwav] = red_time2double((strsplit(date_beg,'T',/extract))[1])
      tend_array[iwav] = red_time2double((strsplit(date_end,'T',/extract))[1])
      tavg_array[iwav] = red_time2double((strsplit(date_avg,'T',/extract))[1])

      ;; Wavelength and time
      wcs[iwav, 0].wave = scan_nbstates[iwav].tun_wavelength*1d9
      wcs[iwav, 0].time = tavg_array[iwav, 0]

      ;; Exposure time
      exp_array[iwav]  = fxpar(nbhead, 'XPOSURE')
      sexp_array[iwav] = fxpar(nbhead, 'TEXPOSUR')
      nsum_array[iwav] = fxpar(nbhead, 'NSUMEXP')
      
      ;; Get destretch to anchor camera (residual seeing)
      if wbcor then begin
        wwi = (red_readdata(scan_wbfiles[iwav]))[x0:x1, y0:y1]
        grid1 = red_dsgridnest(wb, wwi, tile, clip)
      endif

      ;; Read image, apply prefilter curve and temporal scaling
      nbim = (red_readdata(scan_nbfiles[iwav]))[x0:x1, y0:y1] * rpref[iwav] 

      if prefilter eq '3950' and ~keyword_set(noaligncont) then begin
        ;; Apply alignment to compensate for time-variable chromatic
        ;; aberrations.
        nbim = red_shift_sub(nbim, -xshifts[iwav], -yshifts[iwav])
      endif

      ;; Apply destretch to anchor camera and prefilter correction
      if wbcor then nbim = red_stretch(temporary(nbim), grid1)

      if keyword_set(integer) then begin
        self -> fitscube_addframe, fileassoc, fix(round((temporary(nbim)-bzero)/bscale)) $
                                   , ituning = iwav
      endif else begin
        self -> fitscube_addframe, fileassoc, temporary(nbim) $
                                   , ituning = iwav
      endelse

    endfor                      ; iwav

    
    ;; Get pointing at center of FOV for the different tunings.
    red_wcs_hpl_coords, tavg_array, metadata_pig, time_pig $
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
    if ~keyword_set(nocavitymap) then self -> fitscube_addcmap, filename, cmap1

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
    
    print, inam + ' : Narrowband cube stored in:'
    print, filename
    
  endfor                        ; iscan

  
end