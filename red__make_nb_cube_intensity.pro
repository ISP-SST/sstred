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
;     intensitycorrmethod : in, optional, type="string or boolean", default='fit'
;
;       Indicate whether to do intensity correction based on WB data
;       and with what method. See documentation for red::fitscube_intensitycorr.
;
;     integer : in, optional, type=boolean
;
;       Store as integers instead of floats. Uses the BZERO and BSCALE
;       keywords to preserve the intensity scaling.
;
;     nearest : in, optional, type=boolean
;       
;       Use nearest neighbor interpolation (default = bilinear interpolation)
;
;     noaligncont : in, optional, type=boolean
;
;       Do not do the align continuum to wideband step.
;
;     nocavitymap : in, optional, type=boolean
;
;       Do not add cavity maps to the WCS metadata.
;
;     nomissing_nans : in, optional, type=boolean 
;
;       Do not set missing-data padding to NaN. (Set it to the median
;       of each frame instead.)
;
;     nostretch : in, optional, type=boolean
;   
;       Compute no intrascan stretch vectors if this is set.
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
;    2020-10-01 : JdlCR. Apply all corrections at once and use
;                 internal interpolation routines
; 
;    2020-11-09 : MGL. New keyword nostretch.
;
;    2021-02-07 : JdlCR. Bugfix concerning the continuum alignment
;                 shifts. They should be applied separately. Now we are
;                 applying them with the WB rubbersheet correction or
;                 (if WB rubersheet correcions is not used) as a shift
;                 of the rotation axis. 
; 
;    2021-12-02 : MGL. Accept new multi-directory wb cubes.
;
;    2021-12-10 : JdlCR. Make use of the new libgrid routines, now
;                 ported to rdx and maintainable by us.
;
;    2021-04-03 : OA. Added clipping of cavity maps with information
;                 from configuration files (needed to make 'mixed' cubes).
;
;    2025-05-14 : MGL. Adapt to new camera alignment model.
;
;    2025-05-23 : MGL. New keyword redemodulate. Hand over to
;                 red::make_nb_cube for CHROMIS in polarimetric mode.
;
;    2025-06-04 : MGL. Renamed chromis::make_nb_cube -->
;                 make_nb_cube_intensity, to be used only to make
;                 non-Stokes cubes.
; 
;-
pro red::make_nb_cube_intensity, wcfile $
                                 , ashifts = ashifts $
                                 , clips = clips $
                                 , cshift_mean = cshift_mean $
                                 , fileassoc = fileassoc $
                                 , filename = filename $
                                 , fov_mask = fov_mask $
                                 , nocavitymap = nocavitymap $
                                 , nthreads = nthreads $
                                 , odir = odir $
                                 , pertuningfiles = pertuningfiles $
                                 , pertuningstates = pertuningstates $
                                 , redemodulate = redemodulate $
                                 , remove_smallscale = remove_smallscale $
                                 , prefilter_curve = prefilter_curve $
                                 , tiles = tiles $
                                 , wbfilename = wbfilename $
                                 , wbfileassoc = wbfileassoc $
                                 , wbcor =  wbcor $
                                 , wcs = wcs $
;
                                 , date_beg_array = date_beg_array  $ 
                                 , date_end_array = date_end_array  $
                                 , date_avg_array = date_avg_array  $
                                 , exp_array      = exp_array       $
                                 , sexp_array     = sexp_array     $ 
                                 , nsum_array     = nsum_array    $  
                                 , fnumsum_array  = fnumsum_array




;;                                 , clips = clips $
;;                                 , cmap_fwhm = cmap_fwhm $
;;                                 , fitpref_time = fitpref_time $
;;                                 , integer = integer $
;;                                 , intensitycorrmethod = intensitycorrmethod $
;;                                 , nearest = nearest $
;;                                 , noaligncont = noaligncont $  
;;                                 , nocavitymap = nocavitymap $
;;                                 , noflipping = noflipping $
;;                                 , nomissing_nans = nomissing_nans $
;;                                 , nostretch = nostretch $
;;                                 , nthreads = nthreads $
;;                                 , odir = odir $
;;                                 , overwrite = overwrite $
;;                                 , redemodulate = redemodulate $
;;                                 , tiles = tiles $
;;                                 , wbsave = wbsave

  ;; We need to support a few different cases: 1) CHROMIS data without
  ;; polarimetry; 2) CRISP and CHROMIS polarimetric data but just
  ;; summing LC states instead of demodulating.
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  
  
  ;; Camera/detector identification
  self -> getdetectors
  wbindx     = where(strmatch(*self.cameras,'*-W'))
  wbcamera   = (*self.cameras)[wbindx[0]]
  wbdetector = (*self.detectors)[wbindx[0]]
  nbindx     = where(strmatch(*self.cameras,'*-[NTR]')) 
  nbcameras   = (*self.cameras)[nbindx]
  nbdetectors = (*self.detectors)[nbindx]
  Nnbcams = n_elements(nbcameras)

  instrument = (strsplit(wbcamera, '-', /extract))[0]

  polarimetric_data = self -> polarimetric_data() ; Polarimetric data --> non-Stokes cube

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
  Ntunings = n_elements(xindx)

  unbpertuningstates = nbpertuningstates[xindx]
  utunwavelength = unbpertuningstates.tun_wavelength
  sortindx = sort(utunwavelength)
  utunwavelength = utunwavelength[sortindx]
  ufpi_states = ufpi_states[sortindx]
  unbpertuningstates = unbpertuningstates[sortindx]

  ;; Unique nb prefilters
  unbprefs = red_uniquify(unbpertuningstates.prefilter, indx = unbprefindx)
  Nnbprefs = n_elements(unbprefs)

  ;; Get the scan selection 
  uscans = red_uniquify(pertuningstates.scannumber)
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
                               , cam = nbcameras $
                               , sel = nbindx, count = Nnb
  nbstates = pertuningstates[nbindx]
  nbfiles = pertuningfiles[nbindx]
  
  rpref = 1.d0/prefilter_curve

  ;; Spatial dimensions that match the WB cube
  Nx = wcND[0]
  Ny = wcND[1]

  iprogress = 0
  Nprogress = Nscans*Ntunings*Nnbcams

  cshift_mean = fltarr(2, Nnbprefs, Nscans)
  nSummed = fltarr(Nnbprefs, Nscans)
  idxpref = intarr(Ntunings)

  for iscan = 0L, Nscans-1 do begin

    ts = (strsplit(wbgfiles[iscan],'/',/extract))[1]

    ;; The NB files in this scan, sorted in tuning wavelength order.
    self -> selectfiles,  files = pertuningfiles, states = pertuningstates $ 
                                  , cam = nbcameras, scan = uscans[iscan], timestamps = ts $
                                  , sel = scan_nbindx, count = count
    scan_nbfiles = pertuningfiles[scan_nbindx]
    scan_nbstates = pertuningstates[scan_nbindx]
    sortindx = sort(scan_nbstates.tun_wavelength)
    scan_nbfiles = scan_nbfiles[sortindx]
    scan_nbstates = scan_nbstates[sortindx]

    ;; Read global WB file to use as reference when destretching
    ;; per-tuning wb files and then the corresponding nb files.
    wb = (red_readdata(wbgfiles[iscan], direction = direction))[x0:x1, y0:y1]
    

    ;;stop
    ;;if keyword_set(notimecor) then tscl = 1. else tscl = mean(wcTMEAN) / wcTMEAN[iscan]
    ;;tscl *= mean(prefilter_wb)
;    tscl = 1.
    
    for ituning = 0L, Ntunings - 1 do begin 

      ;; Loop over tunings, potentially sum multiple images after
      ;; destretching. Also combine things like exposure times, etc.

      state = ufpi_states[ituning]

      idxprefilter = where(scan_nbstates[ituning].prefilter eq unbprefs)
      idxpref[ituning] = idxprefilter

      nbim = 0.
      wbim = 0.

      undefine, these_texposur, these_xposure
      undefine, these_date_beg, these_date_avg, these_date_end
      
      for icam = 0, Nnbcams-1 do begin

        red_progressbar, iprogress, Nprogress $
                         , /predict $
                         , 'Processing scan=' + strtrim(uscans[iscan], 2) $
                         + ' state=' + state $
                         + ' ' + nbcameras[icam]

        these_nbindx = where(scan_nbstates.fpi_state eq ufpi_states[ituning] $
                             and scan_nbstates.camera eq nbcameras[icam]  $
                             , count)

        if polarimetric_data then begin
          ulc = scan_nbstates[these_nbindx].lc
          Nlc = n_elements(ulc)
          if Nlc gt 1 then begin
            ;; Make sure the LC states are ordered.
            sindx = sort(ulc)
            ulc = ulc[sindx]
            these_nbindx = these_nbindx[sindx]
          endif
        endif else begin
          ulc = ''
          Nlc = 1
        endelse

        
        for ilc = 0, Nlc-1 do begin

          ;; This should be a single NB file
          this_nbindx  = these_nbindx[ilc]
          this_nbfile  = scan_nbfiles[this_nbindx]
          this_nbstate = scan_nbstates[this_nbindx]

          ;; Now find the matching WB file. The selecfiles method
          ;; checks for cameras being equal but with PD data, the
          ;; restored WB camera is actually both the W and D cameras.
          ;; So we need to use strmatch instead. (Perhaps change
          ;; selectfiles to use strmatch?)
          if polarimetric_data then begin
            this_wbindx = where(pertuningstates.fpi_state eq ufpi_states[ituning] $
                                and strmatch(pertuningstates.camera, '*'+wbcamera+'*') $
                                and pertuningstates.lc eq ulc[ilc] $
                                and pertuningstates.scannumber eq uscans[iscan] $
                                , count)
          endif else begin
            this_wbindx = where(pertuningstates.fpi_state eq ufpi_states[ituning] $
                                and strmatch(pertuningstates.camera, '*'+wbcamera+'*') $
                                and pertuningstates.scannumber eq uscans[iscan] $
                                , count)
          endelse
          this_wbfile  = pertuningfiles[this_wbindx]
          this_wbstate = pertuningfiles[this_wbindx]

          if 0 then begin
            print
            print, 'Matching nb and wb files:'
            print, this_nbfile
            print, this_wbfile
            print
          endif
          
          ;; Get destretch to anchor camera (residual seeing)
          if wbcor then begin
            this_wbim = (red_readdata(this_wbfile, direction = direction))[x0:x1, y0:y1]
            wb_grid = rdx_cdsgridnest(wb, this_wbim, tiles, clips)
          endif

          ;; Read NB image(s), apply prefilter curve (rpref) and
          ;; temporal scaling
          this_nbim = (red_readdata(this_nbfile, direction = direction, head = nbhdr))[x0:x1, y0:y1] $
                      * rpref[ituning] ;* tscl
;          bg = median(this_nbim)

          if n_elements(fov_mask) gt 0 then this_nbim *= fov_mask                    
          
          ;; Accumulate header keywords here
          red_append, these_xposure,  fxpar(nbhdr, 'XPOSURE')
          red_append, these_texposur, fxpar(nbhdr, 'TEXPOSUR')

          ;; DATE-??? keywords
          red_fitspar_getdates, nbhdr $
                                , date_beg = date_beg $
                                , date_avg = date_avg $
                                , date_end = date_end 
          red_append, these_date_beg, date_beg
          red_append, these_date_avg, date_avg
          red_append, these_date_end, date_end
          
;          if prefilter eq '3950' and ~keyword_set(noaligncont) then begin
;            ;; Apply alignment to compensate for time-variable chromatic
;            ;; aberrations.
;                                ;nbim = red_shift_sub(nbim, -xshifts[ituning], -yshifts[ituning], )
;            ;; add shift to total shift
          idx_shift = wcSHIFT[0,iscan]
          idy_shift = wcSHIFT[1,iscan]
          cshift = ashifts[*, ituning, iscan]
;            cshift = [xshifts[ituning], yshifts[ituning]]
                                ;          endif else begin
;            ;; With ashifts, this is actually the same as above, since
;            ;; ashift[*,*] = 0 for non-3950 data:
;            idx_shift = wcSHIFT[0,iscan] 
;            idy_shift = wcSHIFT[1,iscan]
;            cshift = [0,0]
;          endelse

          cshift_mean[*, idxprefilter, iscan] += cshift
          nSummed[idxprefilter, iscan] += 1
          
          ;; If needed, accumulate destretch to anchor camera and
          ;; prefilter correction
          if wbcor then begin
            grid1 = wb_grid
            grid1[0,*,*] += cshift[0]
            grid1[1,*,*] += cshift[1]
            this_nbim = rdx_cstretch(temporary(this_nbim), grid1, nthreads=nthreads)
            unrotated_shifts = [0,0]
          endif else begin
            unrotated_shifts = cshift
          endelse

          
          ;; Apply field derotation, alignments, dewarping based on
          ;; the output from make_wb_cube and the per-tuning WB files.

          red_missing, this_nbim $
                       , nmissing = Nmissing, indx_missing = indx_missing, indx_data = indx_data
          if Nmissing gt 0 then begin
            bg = this_nbim[indx_missing[0]]
          endif else begin
            bg = median(this_nbim)
          endelse
          
          this_nbim = red_rotation(temporary(this_nbim), ang[iscan] $
                                   , idx_shift, idy_shift, full=wcFF $
                                   , background = bg, nearest = nearest, nthreads = nthreads $
                                   , stretch_grid = reform(wcGRID[iscan,*,*,*])*sclstr $
                                   , unrotated_shifts = unrotated_shifts)

          ;; Accumulate while normalizing with number of summed images.
          nbim += this_nbim / (Nlc*Nnbcams)
          
          if keyword_set(wbsave) then begin
            ;; Same operations as on narrowband image, except for
            ;; "aligncont".
            this_wbim = rdx_cstretch(temporary(this_wbim), wb_grid, nthreads=nthreads)
            bg = median(this_wbim)
            this_wbim = red_rotation(temporary(this_wbim), ang[iscan] $
                                     , wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF $
                                     , background = bg,  stretch_grid=reform(wcGRID[iscan,*,*,*])*sclstr $
                                     , nthreads=nthreads, nearest = nearest)
            ;; No need to use unrotated_shifts as it only concerns the NB
            ;; channel.
            
            wbim += this_wbim / (Nlc*Nnbcams)
            
          endif
          
        endfor                  ; ilc

        iprogress++             ; update progress counter
        
      endfor                    ; icam
      
      ;; The metadata below should take into account when multiple NB
      ;; images are summed. This is not done yet!   <------------------------- !!!!!!!!!!!!
      
      ;; Collect info about this frame here.
      these_nbindx = where(scan_nbstates.fpi_state eq ufpi_states[ituning]  $
                           , count)
      nbhead = red_sumheaders(scan_nbfiles[these_nbindx], nbim)
      fns = red_fitsgetkeyword_multifile(scan_nbfiles[these_nbindx],'FNUMSUM')
      undefine, fnumsum
      for i = 0, n_elements(fns)-1 do red_append, fnumsum, rdx_str2ints(fns[i])
      red_append, fnumsum_total, fnumsum
      Nsumexp = n_elements(fnumsum)
      fnumsum = rdx_ints2str(red_uniquify(fnumsum))
;      fxaddpar, nbhead, 'FNUMSUM', fnumsum, 'Raw frame numbers'
      ;; fxaddpar, nbhead, 'NSUMEXP', Nsumexp, 'Number of summed exposures'
      
      ;;nbhead = red_readhead(scan_nbfiles[ituning])

      ;; DATE-??? keywords
;      red_fitspar_getdates, nbhead $
;                            , date_beg = date_beg $
;                            , date_end = date_end $
;                            , date_avg = date_avg
      date_beg_array[ituning, iscan] = min(these_date_beg)
      date_end_array[ituning, iscan] = max(these_date_end)

      ;; Exposure time
      exp_array[ituning, iscan]     = total(these_xposure)   ; fxpar(nbhead, 'XPOSURE')
      sexp_array[ituning, iscan]    = median(these_texposur) ; fxpar(nbhead, 'TEXPOSUR')
      nsum_array[ituning, iscan]    = Nsumexp                ; fxpar(nbhead, 'NSUMEXP')
      fnumsum_array[ituning, iscan] = fnumsum
      
      red_fitscube_addframe, fileassoc, temporary(nbim) $
                             , iscan = iscan, ituning = ituning
      
      if keyword_set(wbsave) then begin
        red_fitscube_addframe, wbfileassoc, temporary(wbim) $
                               , iscan = iscan, ituning = ituning
      endif
      
    endfor                      ; ituning

    cshift_mean[0,*,iscan] /= nSummed[*,iscan]
    cshift_mean[1,*,iscan] /= nSummed[*,iscan]
    
  endfor                        ; iscan
  
  ;; String version of all frame numbers in the cube, sorted.
  fnumsum_total = rdx_ints2str(red_uniquify(fnumsum_total))

;;  ;; Create cubes for science data and scan-adapted cavity maps.
;;  cavitymaps = fltarr(Nx, Ny, 1, 1, Nscans)

  ;; Close fits file.
  lun = (size(fileassoc,/struc)).file_lun
  free_lun, lun

;;  if ~keyword_set(nocavitymap) then begin
;;
;;    ;; Read the original cavity map
;;    pindx = where(nbstates.prefilter ne '3999') ; No cavity map for the Ca II H continuum
;;    pindx = pindx[uniq(nbstates[pindx].prefilter, sort(nbstates[pindx].prefilter))]
;;    cprefs = nbstates[pindx].prefilter
;;    Ncprefs = n_elements(cprefs)
;;    
;;    for icprefs = 0, Ncprefs-1 do begin
;;      cfile = self.out_dir + 'flats/spectral_flats/' $
;;              + strjoin([nbstates[pindx[icprefs]].detector $
;;                         , nbstates[pindx[icprefs]].cam_settings $
;;                         , cprefs[icprefs] $
;;                         , 'fit_results.sav'] $
;;                        , '_')
;;
;;      if ~file_test(cfile) then begin
;;        print, inam + ' : Error, calibration file not found -> '+cfile
;;        print, 'Please run the fitprefilter for '+cprefs[icprefs]+' or continue without'
;;        print, 'cavity map for '+cprefs[icprefs]
;;        stop
;;      endif
;;      restore, cfile                 ; The cavity map is in a struct called "fit". 
;;      cmap = reform(fit.pars[1,*,*]) ; Unit is [Angstrom]
;;      cmap /= 10.                    ; Make it [nm]
;;      cmap = -cmap                   ; Change sign so lambda_correct = lambda + cmap
;;      fit = 0B                       ; Don't need the fit struct anymore.
;;      
;;      if keyword_set(remove_smallscale) then begin
;;        ;; If the small scale is already corrected, then include only the
;;        ;; low-resolution component in the metadata. The blurring kernel
;;        ;; should match how the low resolution component was removed when
;;        ;; making flats.
;;        npix = 30               ; Can we get this parameter from earlier headers?
;;        cpsf = red_get_psf(npix*2-1,npix*2-1,double(npix),double(npix))
;;        cpsf /= total(cpsf, /double)
;;        cmap = red_convolve(temporary(cmap), cpsf)
;;        cmap1 = cmap
;;      endif else begin
;;        ;; If the small scale is not already corrected, then we still want
;;        ;; to blur the cavity map slightly.
;;        npsf = round(cmap_fwhm * 7.)
;;        if((npsf/2)*2 eq npsf) then npsf += 1L
;;        psf = red_get_psf(npsf, npsf, cmap_fwhm, cmap_fwhm)
;;        psf /= total(psf, /double)
;;        ;; Leave the orignal cmap alone, we might need it later.
;;        cmap1 = red_convolve(cmap, psf)
;;      endelse
;;
;;      ;; Read the output of the pinhole calibrations so we can do the same
;;      ;; to the cavity maps as was done to the raw data in the momfbd
;;      ;; step. This output is in a struct "alignments" in the save file
;;      ;; 'calib/alignments.sav'
;;;      restore,'calib/alignments.sav'
;;;      ;; Should be based on state1 or state2 in the struct? make_cmaps
;;;      ;; says "just pick one close to continuum (last state?)".
;;;      indx = where(nbstates[0].prefilter eq alignments.state2.prefilter, Nalign)
;;;      case Nalign of
;;;        0    : stop             ; Should not happen!
;;;        1    : amap = invert(      alignments[indx].map           )
;;;        else : amap = invert( mean(alignments[indx].map, dim = 3) )
;;;      endcase
;;
;;      cmap1 = red_apply_camera_alignment(cmap1 $
;;                                         , alignment_model, instrument+'-N' $
;;                                         , prefilter = nbstates[0].prefilter $
;;                                         , amap = amap $
;;                                         , /preserve_size)
;;;      cmap1 = rdx_img_project(amap, cmap1, /preserve) ; Apply the geometrical mapping      
;;
;;      if self.filetype eq 'MOMFBD' then begin
;;        ;; Crop the cavity map to the FOV of the momfbd-restored images.
;;        mr = momfbd_read(wbgfiles[0],/nam)
;;        cmap1 = red_crop_as_momfbd(cmap1, mr)
;;      endif else begin ;; Crop with information from the cfg file
;;        spl = strsplit(wbgfiles[0],'/',/extract)
;;        cw = where(strmatch(spl,'*cfg*'))
;;        cfg_dir=strjoin(spl[0:cw],'/')
;;        cfg_file = cfg_dir+'/'+'momfbd_reduc_'+wbgstates[0].prefilter+'_'+$
;;                   string(wbgstates[0].scannumber,format='(I05)')+'.cfg'
;;        cfg = redux_readcfg(cfg_file)
;;        num_points = long(redux_cfggetkeyword(cfg, 'NUM_POINTS'))
;;        margin = num_points/8
;;        sim_xy = redux_cfggetkeyword(cfg, 'SIM_XY', count = cnt)
;;        if cnt gt 0 then begin
;;          sim_xy = rdx_str2ints(sim_xy)
;;          indx = indgen(n_elements(sim_xy)/2)*2
;;          indy = indx+1
;;          sim_x = sim_xy[indx]
;;          sim_y = sim_xy[indy]   
;;        endif else begin
;;          sim_x = rdx_str2ints(redux_cfggetkeyword(cfg, 'SIM_X'))
;;          sim_y = rdx_str2ints(redux_cfggetkeyword(cfg, 'SIM_Y'))
;;        endelse
;;        xx0 = min(sim_x) + margin - num_points/2 
;;        xx1 = max(sim_x) - margin + num_points/2 - 1
;;        yy0 = min(sim_y) + margin - num_points/2 
;;        yy1 = max(sim_y) - margin + num_points/2 - 1
;;        cmap1 = cmap1[xx0:xx1,yy0:yy1]
;;      endelse
;;      ;; Get the orientation right
;;      cmap1 = red_rotate(cmap1, direction)
;;      ;; Clip to the selected FOV
;;      cmap1 = cmap1[x0:x1,y0:y1]
;;
;;      if file_test(cfg_dir+'/fov_mask.fits') then begin
;;        fov_mask = readfits(cfg_dir+'/fov_mask.fits')
;;        if self.filetype eq 'MOMFBD' then $
;;           fov_mask = red_crop_as_momfbd(fov_mask, mr) $
;;        else $
;;           fov_mask = fov_mask[xx0:xx1,yy0:yy1]
;;        fov_mask = red_rotate(fov_mask, direction)
;;      endif
;;      
;;      
;;      ;; Now make rotated copies of the cavity map
;;      for iscan = 0L, Nscans-1 do begin
;;
;;        if ~keyword_set(nocavitymap) then begin
;;
;;          if ~keyword_set(nomissing_nans) then bg=!Values.F_NaN
;;          ;; Apply the same derot, align, dewarp as for the science data
;;          cmap11 = red_rotation(cmap1, ang[iscan], $
;;                                wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF, $
;;                                stretch_grid=reform(wcGRID[iscan,*,*,*])*sclstr, $
;;                                nthreads=nthreads, nearest=nearest, $
;;                                unrotated_shifts = cshift_mean[*,icprefs,iscan], $
;;                                background=bg)          
;;
;;          cavitymaps[0, 0, 0, 0, iscan] = cmap11
;;
;;          ;; The following block of code is inactive but we want to keep
;;          ;; it around in case it is needed later. It does further
;;          ;; operations on the cavity maps based on what is done to the
;;          ;; science data, such as stretching based on the extra
;;          ;; per-tuning wideband objects, as well as blurring based on the
;;          ;; momfbd psfs. It should probably have a boolean keyword that
;;          ;; activates it. (This code will have to be updated to the
;;          ;; current pipeline style before it can be used.)
;;          if 0 then begin
;;            wb = (red_readdata(wbf[ss]))[x0:x1,y0:y1]
;;            for ww = 0L, nw - 1 do begin
;;              
;;              iwbf = strjoin([self.camwbtag, (strsplit(file_basename(st.ofiles[ww,ss]) $
;;                                                       , '.', /extract))[1:*]],'.')
;;              iwbf = file_dirname(st.ofiles[ww,ss]) + '/'+iwbf
;;              
;;              ;; Load images
;;              iwb = (red_readdata(iwbf))[x0:x1,y0:y1]
;;              im = momfbd_read(st.ofiles[ww,ss])
;;              
;;              ;; get dewarp from WB
;;              igrid = rdx_cdsgridnest(wb, iwb, itiles, iclip)
;;              
;;              ;; Convolve CMAP and apply wavelength dep. de-warp
;;              cmap2 = rdx_cstretch((red_mozaic(red_conv_cmap(cmap, im), /crop))[x0:x1, y0:y1], igrid, $
;;                                   nthreads=nthreads)
;;              
;;              ;; Derotate and shift
;;              cmap2 = red_rotation(temporary(cmap2), ang[ss], total(shift[0,ss]), $
;;                                   total(shift[1,ss]), full=wcFF, stretch_grid=reform(grid[ss,*,*,*]), $
;;                                   nearest=nearest, nthreads=nthreads,unrotated_shifts = cshift_mean[*,icprefs,iscan])
;;              
;;              ;; Time de-warp
;;              ;;cmap2 = red_stretch(temporary(cmap2), reform(grid[ss,*,*,*]))
;;              
;;            endfor              ; ww
;;          endif
;;        endif
;;
;;      endfor                    ; iscan
;;
;;      tindx = where(scan_nbstates.prefilter eq cprefs[icprefs], Nt)
;;      
;;      ;; Add cavity maps as WAVE distortions
;;      if Nt gt 0 then begin
;;        red_fitscube_addcmap, filename, cavitymaps $
;;                              , cmap_number = icprefs+1 $
;;                              , prefilter = cprefs[icprefs] $
;;                              , indx = tindx
;;      endif
;;      
;;    endfor                      ; icprefs
;;    
;;  endif
  
  
end

