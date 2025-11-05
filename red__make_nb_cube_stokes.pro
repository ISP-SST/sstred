; docformat = 'rst'

;+
; Make a de-rotated and de-stretched time-series FITS data cube with
; Stokes images based on momfbd-restored narrow-band images.
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
;    clips : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment.
;
;    cmap_fwhm : in, type=float, default=7
;   
;       FWHM in pixels of kernel used for smoothing the cavity map.
;
;    integer : in, optional, type=boolean
;
;       Store as integers instead of floats. Uses the BZERO and BSCALE
;       keywords to preserve the intensity scaling.
;
;    intensitycorrmethod : in, optional, type="string or boolean", default='fit'
;
;       Indicate whether to do intensity correction based on WB data
;       and with what method. See documentation for red::fitscube_intensitycorr.
;
;    nearest : in, optional, type=boolean
;       
;       Use nearest neighbor interpolation (default = bilinear interpolation)
;
;    nocavitymap : in, optional, type=boolean
;
;       Do not add cavity maps to the WCS metadata.
;
;    nocrosstalk : in, optional, type=boolean
;
;       Do not correct the (polarimetric) data cube Stokes components
;       for crosstalk from I to Q, U, V.
;
;    nomissing_nans : in, optional, type=boolean 
;
;      Do not set missing-data padding to NaN. (Set it to the median of
;      each frame instead.)
;
;    noremove_periodic : in, optional, type=boolean
;
;       Do not attempt to use Fourier filtering to remove polcal
;       periodic artifacts. (See crisp::demodulate method.)
;
;    nostretch : in, optional, type=boolean
;   
;      Compute no intrascan stretch vectors if this is set.
;
;    overwrite : in, optional, type=boolean
;
;       Don't care if cube is already on disk, overwrite it
;       with a new version.
;
;    redemodulate : in, optional, type=boolean
;
;       Delete any old (per scan-and-tuning) stokes cubes so they will
;       have to be demodulated from scratch.
;
;    tiles : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment. 
;
;    wbsave : in, optional, type=boolean
;
;       Save a cube with the wideband per-tuning align-images. For
;       debugging of alignment with extra wideband objects.
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
;    2018-10-08 : MGL. New keyword redemodulate.
; 
;    2018-12-21 : MGL. New keyword smooth.
; 
;    2019-03-28 : MGL. New keyword demodulate_only.
; 
;    2019-04-02 : MGL. New keyword nocrosstalk and use the
;                 fitscube_crosstalk method.
; 
;    2019-05-24 : MGL. Do demodulation by use of the make_stokes_cubes
;                 method. Remove keyword demodulate_only.
; 
;    2020-01-16 : MGL. New keyword intensitycorrmethod.
; 
;    2020-02-18 : MGL. Remove CHROMIS keyword noaligncont.
; 
;    2020-03-11 : MGL. Apply rotate() with direction from WB cube. 
; 
;    2020-06-16 : MGL. Remove temporal intensity scaling, deprecate
;                 keyword notimecorr.
; 
;    2020-07-15 : MGL. Remove keyword smooth.
;
;    2020-10-01 : JdlCR. Modified to use new rotation routine that
;                 also applies the stretch_grid.
;
;    2020-10-28 : MGL. Remove statistics calculations.
; 
;    2020-11-09 : MGL. New keyword nostretch.
; 
;    2021-12-02 : MGL. Accept new multi-directory wb cubes.
;
;    2021-12-10 : JdlCR. Make use of the new libgrid routines, now
;                 ported to rdx and maintainable by us.
; 
;    2022-09-10 : MGL. CRISP --> RED.
;
;    2023-04-04 : OA. Added clipping of cavity maps with information
;                 from configuration files (needed to make 'mixed'
;                 cubes).
;
;    2025-02-20 : MGL. Adapt to new camera alignment model.
;
;    2025-05-23 : MGL. New keyword noaligncont.
;
;    2025-06-03 : MGL. Renamed make_nb_cube --> make_nb_cube_stokes,
;                 to be used only to make cubes with all 4 Stokes
;                 components.
;
;-
pro red::make_nb_cube_stokes, wcfile $
                              , ashifts = ashifts $
                              , clips = clips $
                              , cshift_mean = cshift_mean $
                              , date_avg_array = date_avg_array  $
                              , date_beg_array = date_beg_array  $ 
                              , date_end_array = date_end_array  $
                              , exp_array = exp_array $
                              , fileassoc = fileassoc $
                              , filename = filename $
                              , fitpref_time = fitpref_time $
                              , fnumsum_array  = fnumsum_array $
                              , fov_mask = fov_mask $
                              , nocavitymap = nocavitymap $
                              , nsum_array = nsum_array    $  
                              , nthreads = nthreads $
                              , odir = odir $
                              , pertuningfiles = pertuningfiles $
                              , pertuningstates = pertuningstates $
                              , redemodulate = redemodulate $
                              , remove_smallscale = remove_smallscale $
                              , sexp_array = sexp_array     $ 
                              , tiles = tiles $
                              , wbcor = wbcor $
                              , wbfileassoc = wbfileassoc $
                              , wbfilename = wbfilename $
                              , wcs = wcs
                              
;;                              , clips = clips  $
;;                              , cmap_fwhm = cmap_fwhm $
;;                              , integer = integer $
;;                              , intensitycorrmethod = intensitycorrmethod $ 
;;                              , nearest = nearest $
;;                              , noaligncont = noaligncont $  
;;                              , nocrosstalk = nocrosstalk $
;;                              , noflipping = noflipping $
;;                              , nomissing_nans = nomissing_nans $
;;                              , noremove_periodic = noremove_periodic $
;;                              , nostretch = nostretch $
;;                              , nthreads = nthreads $
;;                              , odir = odir $
;;                              , overwrite = overwrite $
;;                              , redemodulate = redemodulate $
;;                              , tiles = tiles $
;;                              , wbsave = wbsave

                              ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  ;; Camera/detector identification
  self -> getdetectors
  wbindx      = where(strmatch(*self.cameras,'*-W'))
  wbcamera    = (*self.cameras)[wbindx[0]]
  wbdetector  = (*self.detectors)[wbindx[0]]
  nbtindx     = where(strmatch(*self.cameras,'*-T')) 
  nbtcamera   = (*self.cameras)[nbtindx[0]]
  nbtdetector = (*self.detectors)[nbtindx[0]]
  nbrindx     = where(strmatch(*self.cameras,'*-R')) 
  nbrcamera   = (*self.cameras)[nbrindx[0]]
  nbrdetector = (*self.detectors)[nbrindx[0]]

  instrument = (strsplit(wbcamera, '-', /extract))[0]
  
  if instrument eq 'Chromis' && self.isodate gt red_dates(tag = 'CHROMIS Ximea') then begin
    ;; This can be remove when Pit has had time to do telescope
    ;; polarimetry calibration for Chromis
    notelmat = 1
  endif 

  ;; Read parameters from the WB cube
  fxbopen, bunit, wcfile, 'MWCINFO', bbhdr
  fxbreadm, bunit, row = 1 $
            , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01', 'DIRECTION'] $
            ,   ANG, wcCROP, wcFF, wcGRID, wcND, wcSHIFT, wcTMEAN, wcX01Y01,   direction
  ;; Note that the strarr wfiles cannot be read by fxbreadm! Put it in
  ;; wbgfiles (WideBand Global).
  fxbread, bunit, wbgfiles, 'WFILES', 1
  fxbclose, bunit
  if self.filetype eq 'MIXED' then wbgfiles = strtrim(wbgfiles, 2)

  ;; Spatial dimensions that match the WB cube
  Nx = wcND[0]
  Ny = wcND[1]

  ;; Get the alignment mode, projective or polywarp.
  alignment_model = red_align_model(wbgfiles[0])

  ;; Don't do any stretching if wcgrid is all zeros.
  nostretch_temporal = total(abs(wcgrid)) eq 0
  sclstr = 0
  if ~nostretch_temporal then sclstr = 1

  x0 = wcX01Y01[0]
  x1 = wcX01Y01[1]
  y0 = wcX01Y01[2]
  y1 = wcX01Y01[3]
  origNx = x1 - x0 + 1
  origNy = y1 - y0 + 1


  ;; Unique tuning states, sorted by wavelength
  nbpertuningstates = pertuningstates[where(~pertuningstates.is_wb)]
  ufpi_states = red_uniquify(nbpertuningstates.fpi_state, indx = xindx)
  Ntuning = n_elements(xindx)

  unbpertuningstates = nbpertuningstates[xindx]
  ufpi_states = red_uniquify(nbpertuningstates.fpi_state, indx = xindx)
  Ntuning = n_elements(xindx)

  unbpertuningstates = nbpertuningstates[xindx]
  utunwavelength = unbpertuningstates.tun_wavelength
  sortindx = sort(utunwavelength)
  utunwavelength = utunwavelength[sortindx]
  ufpi_states = ufpi_states[sortindx]
  unbpertuningstates = unbpertuningstates[sortindx]
;;  
;;  my_prefilters = pertuningstates[xindx].prefilter
;;  wav = utunwavelength

  ;; Unique nb prefilters
  unbprefs = red_uniquify(unbpertuningstates.prefilter, indx = unbprefindx)
  Nnbprefs = n_elements(unbprefs)

;  ;; Unique nb prefilters
;  unbprefs = red_uniquify(pertuningstates[where(~pertuningstates.is_wb)].prefilter, indx = unbprefindx)
;  Nnbprefs = n_elements(unbprefs)

  ;; Get the scan selection from wfiles (from the sav file)
  self -> extractstates, wbgfiles, wbgstates
  uscans = wbgstates.scannumber
  Nscans = n_elements(uscans)


  ;; Define the Stokes file names needed for this nb cube, make the
  ;; Stokes cubes if needed.

  snames = strarr(Nscans, Ntuning)
  for iscan = 0, Nscans-1 do begin
    
    self -> make_stokes_cubes, file_dirname(wbgfiles[iscan]), uscans[iscan] $
       , clips = clips $
       , cmap_fwhm = cmap_fwhm $
       , /nocavitymap $         ; Cavity maps in Stokes cubes aren't used for anything
       , notelmat = notelmat $ 
       , noremove_periodic = noremove_periodic $
       , redemodulate = redemodulate $
       , snames = these_snames $
       , stokesdir = stokesdir $
       , tiles = tiles $
       , nearest = nearest $
       , nthreads = nthreads $
       , fitpref_time = fitpref_time

;    if n_elements(snames) eq 0 then begin
;      Ntuning = n_elements(these_snames)
;      snames = strarr(Nscans, Ntuning)
;    endif
    
    snames[iscan, *] = these_snames
    
  endfor                        ; iscan

  self -> extractstates, snames, sstates


  
;;  if ~keyword_set(nocavitymap) then begin
;;
;;    ;; There should be one cavity map per prefilter. Read them!
;;
;;    cavitymaps_pref = fltarr(origNx, origNy, Nnbprefs)
;;    
;;    for ipref = 0, Nnbprefs-1 do begin
;;
;;      red_progressbar, ipref, Nnbprefs, 'Read cavity maps for '+unbprefs[ipref]
;;      
;;      ;; Read the original cavity map
;;      cfile = self.out_dir + 'flats/spectral_flats/' $
;;              + strjoin([nbtdetector $
;;                         , unbprefs[ipref] $
;;                         , 'fit_results.sav'] $
;;                        , '_')
;;      
;;      if ~file_test(cfile) then begin
;;        red_message, 'Error, calibration file not found -> '+cfile
;;        stop
;;      endif
;;      restore, cfile            ; The cavity map is in a struct called "fit". 
;;      dims = size(fit.pars, /dim)
;;
;;      if dims[0] gt 1 then begin
;;        cmapt = reform(fit.pars[1,*,*]) ; Unit is [Angstrom]
;;        cmapt /= 10.                    ; Make it [nm]
;;        cmapt = -cmapt                  ; Change sign so lambda_correct = lambda + cmap
;;        fit = 0B                        ; Don't need the fit struct anymore.
;;        
;;        cfile = self.out_dir + 'flats/spectral_flats/' $
;;                + strjoin([nbrdetector $
;;                           , unbprefs[ipref] $
;;                           , 'fit_results.sav'] $
;;                          , '_')
;;        
;;        if ~file_test(cfile) then begin
;;          red_message, 'Error, calibration file not found -> '+cfile
;;          stop
;;        endif
;;        restore, cfile                  ; The cavity map is in a struct called "fit". 
;;        cmapr = reform(fit.pars[1,*,*]) ; Unit is [Angstrom]
;;        cmapr /= 10.                    ; Make it [nm]
;;        cmapr = -cmapr                  ; Change sign so lambda_correct = lambda + cmap
;;        fit = 0B                        ; Don't need the fit struct anymore.
;;        
;;        if keyword_set(remove_smallscale) then begin
;;          ;; If the small scale is already corrected, then include only the
;;          ;; low-resolution component in the metadata. The blurring kernel
;;          ;; should match how the low resolution component was removed when
;;          ;; making flats.
;;          npix = 30             ; Can we get this parameter from earlier headers?
;;          cpsf = red_get_psf(npix*2-1,npix*2-1,double(npix),double(npix))
;;          cpsf /= total(cpsf, /double)
;;          cmapr = red_convolve(temporary(cmapr), cpsf)
;;          cmap1r = cmapr
;;          cmapt = red_convolve(temporary(cmapt), cpsf)
;;          cmap1t = cmapt
;;        endif else begin
;;          ;; If the small scale is not already corrected, then we still want
;;          ;; to blur the cavity map slightly.
;;          npsf = round(fwhm * 7.)
;;          if((npsf/2)*2 eq npsf) then npsf += 1L
;;          psf = red_get_psf(npsf, npsf, fwhm, fwhm)
;;          psf /= total(psf, /double)
;;          ;; Leave the orignal cmaps alone, we might need them later.
;;          cmap1r = red_convolve(cmapr, psf)
;;          cmap1t = red_convolve(cmapt, psf)
;;        endelse
;;        
;;        ;; Apply geometrical transformation from the pinhole calibration to the cavity maps.
;;        cmap1t = red_apply_camera_alignment(cmap1t, alignment_model, instrument+'-T' $
;;                                            , pref =  unbprefs[ipref] $
;;                                            , amap = amapt $
;;                                            , /preserve_size)
;;        cmap1r = red_apply_camera_alignment(cmap1r, alignment_model, instrument+'-R' $
;;                                            , pref =  unbprefs[ipref] $
;;                                            , amap = amapr $
;;                                            , /preserve_size)
;;
;;        ;; At this point, the individual cavity maps should be corrected
;;        ;; for camera misalignments, so they should be aligned with
;;        ;; respect to the cavity errors on the etalons. So we can sum
;;        ;; them.
;;        cmap1 = (cmap1r + cmap1t) / 2.
;;
;;        if self.filetype eq 'MOMFBD' then begin
;;          ;; Crop the cavity map to the FOV of the momfbd-restored images.
;;          cmap1 = red_crop_as_momfbd(cmap1, mr)
;;        endif else begin
;;          cmap1 = cmap1[xx0:xx1,yy0:yy1]
;;        endelse
;;        
;;
;;        
;;        ;; Clip to the selected FOV
;;        cmap1r = cmap1r[x0:x1,y0:y1]
;;        cmap1t = cmap1t[x0:x1,y0:y1]
;;
;;        ;; Get the orientation right.
;;        cmap1 = red_rotate(cmap1, direction)
;;
;;;    Are we doing things in the right order?
;;;        
;;;    cmap1 = (cmap1r + cmap1t) / 2.
;;;
;;;    if self.filetype eq 'MOMFBD' then begin
;;;      ;; Crop the cavity map to the FOV of the momfbd-restored images.
;;;      cmap1 = red_crop_as_momfbd(cmap1, mr)
;;;    endif else begin
;;;      cmap1 = cmap1[xx0:xx1,yy0:yy1]
;;;    endelse
;;;    
;;;    ;; Get the orientation right.
;;;    cmap1 = red_rotate(cmap1, direction)
;;;
;;;    ;; Clip to the selected FOV
;;;    cmap1 = cmap1[x0:x1,y0:y1]
;;
;;        
;;        cavitymaps_pref[*, *, ipref] = cmap1
;;        
;;      endif                     ;else undefine, cmap1
;;    endfor                      ; ipref
;;  endif                         ;else undefine, cmap1
  
  ;; Unique tuning states, sorted by wavelength
  utunindx = uniq(sstates.fpi_state, sort(sstates.fpi_state))
  Ntuning = n_elements(utunindx)
  sortindx = sort(sstates[utunindx].tun_wavelength)
  utuning = sstates[utunindx[sortindx]].fpi_state

  Nstokes = 4

;;  ;; Create cubes for science data and scan-adapted cavity maps.
;;  cavitymaps = fltarr(Nx, Ny, 1, 1, Nscans)

;  if ~keyword_set(nocavitymap) then begin
;
;    ;; Read the original cavity map
;    cfile = self.out_dir + 'flats/spectral_flats/' $
;            + strjoin([nbtstates[0].detector $
;                       , prefilter $
;                       , 'fit_results.sav'] $
;                      , '_')
;
;    if ~file_test(cfile) then begin
;      red_message, 'Error, calibration file not found -> '+cfile
;      stop
;    endif
;    restore, cfile                  ; The cavity map is in a struct called "fit". 
;    cmapt = reform(fit.pars[1,*,*]) ; Unit is [Angstrom]
;;    cmapt = rotate(temporary(cmapt), direction)
;    cmapt /= 10.                ; Make it [nm]
;    cmapt = -cmapt              ; Change sign so lambda_correct = lambda + cmap
;    fit = 0B                    ; Don't need the fit struct anymore.
;    
;    cfile = self.out_dir + 'flats/spectral_flats/' $
;            + strjoin([nbrstates[0].detector $
;                       , prefilter $
;                       , 'fit_results.sav'] $
;                      , '_')
;
;    if ~file_test(cfile) then begin
;      red_message, 'Error, calibration file not found -> '+cfile
;      stop
;    endif
;    restore, cfile                  ; The cavity map is in a struct called "fit". 
;    cmapr = reform(fit.pars[1,*,*]) ; Unit is [Angstrom]
;;    cmapr = rotate(temporary(cmapr), direction)
;    cmapr /= 10.                ; Make it [nm]
;    cmapr = -cmapr              ; Change sign so lambda_correct = lambda + cmap
;    fit = 0B                    ; Don't need the fit struct anymore.
;    
;    if keyword_set(remove_smallscale) then begin
;      ;; If the small scale is already corrected, then include only the
;      ;; low-resolution component in the metadata. The blurring kernel
;      ;; should match how the low resolution component was removed when
;      ;; making flats.
;      npix = 30                 ; Can we get this parameter from earlier headers?
;      cpsf = red_get_psf(npix*2-1,npix*2-1,double(npix),double(npix))
;      cpsf /= total(cpsf, /double)
;      cmapr = red_convolve(temporary(cmapr), cpsf)
;      cmap1r = cmapr
;      cmapt = red_convolve(temporary(cmapt), cpsf)
;      cmap1t = cmapt
;    endif else begin
;      ;; If the small scale is not already corrected, then we still want
;      ;; to blur the cavity map slightly.
;      npsf = round(fwhm * 7.)
;      if((npsf/2)*2 eq npsf) then npsf += 1L
;      psf = red_get_psf(npsf, npsf, fwhm, fwhm)
;      psf /= total(psf, /double)
;      ;; Leave the orignal cmaps alone, we might need them later.
;      cmap1r = red_convolve(cmapr, psf)
;      cmap1t = red_convolve(cmapt, psf)
;    endelse
;    
;    ;; Apply geometrical transformation from the pinhole calibration to the cavity maps.
;    cmap1t = red_apply_camera_alignment(cmap1t, alignment_model $
;                                        , *self.cameras,instrument+'-T', amap = amapt)
;    cmap1r = red_apply_camera_alignment(cmap1r, alignment_model $
;                                        , *self.cameras,instrument+'-R', amap = amapr)
;
;    ;; At this point, the individual cavity maps should be corrected
;    ;; for camera misalignments, so they should be aligned with
;    ;; respect to the cavity errors on the etalons. So we can sum
;    ;; them.
;    cmap1 = (cmap1r + cmap1t) / 2.
;
;    if self.filetype eq 'MOMFBD' then begin
;      ;; Crop the cavity map to the FOV of the momfbd-restored images.
;      cmap1 = red_crop_as_momfbd(cmap1, mr)
;    endif else begin
;      cmap1 = cmap1[xx0:xx1,yy0:yy1]
;    endelse
;    
;    ;; Get the orientation right.
;    cmap1 = red_rotate(cmap1, direction)
;
;    ;; Clip to the selected FOV
;    cmap1 = cmap1[x0:x1,y0:y1]
;
;  endif                         ; nocavitymap

  
  iprogress = 0
  Nprogress = Nscans*Ntuning

  cshift_mean = fltarr(2, Nnbprefs, Nscans)
  nSummed = fltarr(Nnbprefs, Nscans)
  idxpref = intarr(Ntuning)

  for iscan = 0L, Nscans-1 do begin

    ;; Read global WB file to use as reference when destretching
    ;; per-tuning wb files and then the corresponding nb files.
    wb = (red_readdata(wbgfiles[iscan], direction = direction))[x0:x1, y0:y1]
    indx_nan = where(~finite(wb),Nnan, complement=data_indx)
    if Nnan then wb[indx_nan] = median(wb[data_indx])
;    ts = (strsplit(wbgfiles[iscan],'/',/extract))[1]



    ;; Scaling parameter with respect to mean wb intensity (mainly
    ;; varying elevation)
;  if keyword_set(notimecor) then begin
;    tscl = replicate(1., Nscans)
;  endif else begin
;    tscl = mean(nbt_prefilter_wb) * mean(wcTMEAN) / wcTMEAN
;  endelse

    
    for ituning = 0L, Ntuning - 1 do begin 

      red_progressbar, iprogress, Nprogress $
                       , /predict $
                       , 'Processing scan=' $
                       + strtrim(uscans[iscan], 2) + ' tuning=' + utuning[ituning] 

      ;; Read the pre-made Stokes cube
      
      tmp = red_readdata(snames[iscan, ituning], head = stokhdr)
      nbdims = size(tmp, /dim)
      if max(direction eq [1, 3, 4, 6]) eq 1 then begin
        nbdims = nbdims[[1, 0, 3]] ; X and Y switched
      endif else begin
        nbdims = nbdims[[0, 1, 3]] 
      endelse 
      nbim = fltarr(nbdims)
      for istokes = 0, Nstokes-1 do begin
        nbim[0, 0, istokes] = rotate(tmp[*, *, 0, istokes], direction)
      endfor                               ; istokes
      nbim = reform(nbim[x0:x1, y0:y1, *]) ;* tscl[iscan]
      
      ;; The Stokes cube is already wbcorrected so we should not
      ;; repeat that here.

      ;; A CHROMIS scan can involve more than a single NB prefilter
      ipref = where(unbprefs eq sstates[iscan,ituning].prefilter)
      
      ;; Get some metadata from the Stokes cube header.
      red_fitspar_getdates, stokhdr $
                            , date_beg = date_beg $
                            , date_avg = date_avg $
                            , date_end = date_end 

      date_beg_array[ituning, iscan] = date_beg
;;      date_avg_array[ituning, iscan] = date_avg
      date_end_array[ituning, iscan] = date_end

      ;; Exposure time
      exp_array[ituning, iscan]  = red_fitsgetkeyword(stokhdr, 'XPOSURE')
      sexp_array[ituning, iscan] = red_fitsgetkeyword(stokhdr, 'TEXPOSUR')
      nsum_array[ituning, iscan] = red_fitsgetkeyword(stokhdr, 'NSUMEXP')
 ;     stop
      fnumsum_array[ituning, iscan] = red_fitsgetkeyword(stokhdr, 'FNUMSUM') ;fnumsum                             <-------------------------
      
      ;; Apply derot, align, dewarp based on the output from
      ;; make_wb_cube

      for istokes = 0, Nstokes-1 do begin

        if n_elements(fov_mask) gt 0 then nbim[*, *, istokes] = nbim[*, *, istokes] * fov_mask

        ;; Missing pixels, if any, are set to zero at this point
        frame = nbim[*, *, istokes]
        
        red_missing, frame, /inplace $
                     , nmissing = Nmissing, indx_missing = indx_missing, indx_data = indx_data
        ;; Now the missing pixels should be NaN
        
        if Nmissing gt 0 then begin
          bg = (nbim[*, *, istokes])[indx_missing[0]]
        endif else begin
          bg = median(nbim[*, *, istokes])
        endelse
        
        frame = red_rotation(frame, ang[iscan], $
                             wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF , $
                             background = bg, nearest = nearest, $
                             stretch_grid = reform(wcGRID[iscan,*,*,*])*sclstr, nthreads=nthreads)
        ;; Missing pixels should still be NaN
        
        ;;if Nwhere gt 0 then frame[mindx] = bg ; Ugly fix, red_stretch destroys the missing data?
        
        red_fitscube_addframe, fileassoc, frame $
                               , iscan = iscan, ituning = ituning, istokes = istokes
      endfor                    ; istokes
      
      if keyword_set(wbsave) then begin
        ;; Same operations as on narrowband image.
;        wbim = wwi * tscl[iscan]
;        wbim = red_stretch(temporary(wbim), grid1)
        wbim = red_rotation(temporary(wbim), ang[iscan] $
                            , wcSHIFT[0,iscan], wcSHIFT[1,iscan] $
                            , full=wcFF $
                            , stretch_grid = reform(wcGRID[iscan,*,*,*])*sclstr $
                            , nthreads=nthreads, nearest = nearest)

        red_fitscube_addframe, wbfileassoc, temporary(wbim) $
                               , iscan = iscan, ituning = ituning, istokes = istokes
      endif
      
      iprogress++               ; update progress counter
      
    endfor                      ; ituning

;;    if ~keyword_set(nocavitymap) then begin
;;      if ~keyword_set(nomissing_nans) then bg=!Values.F_NaN
;;      ;; Apply the same derot, align, dewarp as for the science data
;;      cmap1 = cavitymaps_pref[*, *, ipref]
;;      cmap11 = red_rotation(cmap1, ang[iscan] $
;;                            , wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF $
;;                            , stretch_grid = reform(wcGRID[iscan,*,*,*])*sclstr $
;;                            , nthreads=nthreads, background=bg)
;;      
;;      cavitymaps[0, 0, 0, 0, iscan] = cmap11
;;      
;;    endif 
;;    
  endfor                        ; iscan
  

;;  ;; Add cavity maps as WAVE distortions. Close the file first.
;;  lun = (size(fileassoc,/struc)).file_lun
;;  free_lun, lun
;;  if ~keyword_set(nocavitymap) then red_fitscube_addcmap, filename, cavitymaps
;;  ;; Should we open the file again after?
  
end

