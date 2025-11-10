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
;     intensitycorrmethod : in, optional, type="string or boolean", default='fit'
;
;       Indicate whether to do intensity correction based on WB data
;       and with what method. See documentation for red::fitscube_intensitycorr.
;
;     nearest : in, optional, type=boolean
;       
;       Use nearest neighbor interpolation (default = bilinear interpolation)
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
                                 , date_avg_array = date_avg_array  $
                                 , date_beg_array = date_beg_array  $ 
                                 , date_end_array = date_end_array  $
                                 , exp_array = exp_array $
                                 , fileassoc = fileassoc $
                                 , filename = filename $
                                 , fnumsum_array = fnumsum_array $
                                 , fov_mask = fov_mask $
                                 , nearest = nearest $
                                 , nsum_array = nsum_array    $  
                                 , nthreads = nthreads $
                                 , pertuningfiles = pertuningfiles $
                                 , pertuningstates = pertuningstates $
                                 , prefilter_curve = prefilter_curve $
                                 , remove_smallscale = remove_smallscale $
                                 , sexp_array  = sexp_array     $ 
                                 , tiles = tiles $
                                 , wbcor = wbcor $
                                 , wbfileassoc = wbfileassoc $
                                 , wbfilename = wbfilename $
                                 , wbsave = wbsave $
                                 , wcs = wcs

  ;; We need to support a few different cases: 1) CHROMIS data without
  ;; polarimetry; 2) CRISP and CHROMIS polarimetric data but just
  ;; summing LC states instead of demodulating.
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  
  
  ;; Camera/detector identification
  self -> cameras $
     , instrument = instrument $
     , Nnbcams = Nnbcams $
     , nb_detectors = nbdetectors $
     , nb_cameras = nbcameras $
     , wb_detector = wbdetector $
     , wb_camera = wbcamera

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
          
          cshift = ashifts[*, ituning, iscan]
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
                                   , wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF $
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
  
  ;; Close fits file.
  lun = (size(fileassoc,/struc)).file_lun
  free_lun, lun
  
end

