; docformat = 'rst'

;+
; Make Stokes cubes by demodulating momfbd output files.
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
;     dir : in, type=string
; 
;       The path to the momfbd-restored data.
;
;     scanno : in, type=integer
;
;       The scan number for which to make Stokes cubes.
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
;     nocavitymap : in, optional, type=boolean
;
;       Do not add cavity maps to the WCS metadata.
;
;     noremove_periodic : in, optional, type=boolean
;
;       Do not attempt to use Fourier filtering to remove polcal
;       periodic artifacts. (See crisp::demodulate method.)
;
;     nopolarimetry : in, optional, type=boolean
;
;       For a polarimetric dataset, don't make a Stokes cube.
;       Instead combine all LC states for both cameras into a single
;       NB image per tuning, producing a cube similar to that for a
;       data set without polarimetry. (For a nonpolarlimetric dataset,
;       no effect.)
;
;     nostatistics : in, optional, type=boolean
;  
;       Do not calculate statistics metadata to put in header keywords
;       DATA*. 
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
;     redemodulate : in, optional, type=boolean
;
;       Delete any old (per scan-and-tuning) stokes cubes so they will
;       have to be demodulated from scratch.
;
;     tiles : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment. 
;
;     smooth : in, optional, type=varies, default=5
;
;       How to smooth the modulation matrices? Set to the string
;       'momfbd' to smooth subfield by subfield using the PSFs
;       estimated by momfbd. Set to a number to smooth by a Gaussian
;       kernel of that width. 
;
;     nthreads : in, optional, type = int
;
;       number of threads to use in red_stretch
;
;     nearest : in, optional, type = boolean
;
;       perform nearest neighbor interpolation (default = bilinear)
;
; :History:
; 
;    2019-03-21 : MGL. First version.
; 
;    2020-07-15 : MGL. Remove keyword smooth.
;
;    2020-10-01 : JdlCR. Use the new stretching routines and allow for
;                 optional nearest neighbor interpolation. 
; 
;    2022-09-10 : MGL. CRISP --> RED.
;
;    2025-02-20 : MGL. Adapt to new camera alignment model.
;                                
;-
pro red::make_stokes_cubes, dir, scanno $
                            , clips = clips $
                            , cmap_fwhm = cmap_fwhm $
                            , nocavitymap = nocavitymap $
                            , overwrite = overwrite $
                            , noremove_periodic = noremove_periodic $
                            , notelmat = notelmat $
                            , redemodulate = redemodulate $
                            , snames = snames $
                            , stokesdir = stokesdir $
                            , tiles = tiles $
                            , nthreads = nthreads $
                            , nearest = nearest $
                            , fitpref_time = fitpref_time

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)                                      

  if n_elements(dir) eq 0 then begin
    red_message, 'Please specify the directory with momfbd output.'
    retall
  endif

  if n_elements(scanno) ne 1 then begin
    ;; Make it accept only a single scan number so we can return the
    ;; Stokes file names for that scan.
    red_message, 'Please specify a single scan number and try again.'
    return
  endif
  
  ;; Either keyword means to overwrite existing Stokes cubes with new
  ;; demodulations. 
  redemodulate = keyword_set(overwrite) or keyword_set(redemodulate)

  undefine, snames
  
  ;; Smooth with Gaussian kernel by default
  smooth_by_kernel = 5          ; Default width
  smooth_by_subfield = 0

  ;; Store Stokes cubes in separate directories for different smooth
  ;; options:
  stokesdir = dir + '/stokes_sbs'+strtrim(smooth_by_subfield,2) $
              + '_sbk'+strtrim(smooth_by_kernel,2)+'/'
  
  file_mkdir, stokesdir

  ;; Default keywords
  if n_elements(cmap_fwhm) eq 0 then fwhm = 7.0
  if n_elements(tiles) eq 0 or n_elements(clips) eq 0 then begin
    tiles = [12, 16, 32, 64, 72]
    clips = [12, 8,  4,  2 , 1]
  endif

  ;; Camera/detector identification
  self->getdetectors
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

  
  case instrument of
    'Crisp' : begin
;; We currently do correct for the small scale cavity map in CRISP
      ;; data. (We should get this from earlier meta data!)
      remove_smallscale = 1
    end
    'Chromis' : begin
      remove_smallscale = 0
    end
  endcase

  ;; Get metadata from logfiles
  red_logdata, self.isodate, time_r0, r0 = metadata_r0, ao_lock = ao_lock
  red_logdata, self.isodate, time_pointing, diskpos = metadata_pointing, rsun = rsun
  red_logdata, self.isodate, time_turret, azel = azel

  
  
  case self.filetype of
    'ANA': extension = '.f0'
    'MOMFBD': extension = '.momfbd'
    'FITS': extension = '.fits'
    'MIXED' : extension = '.{fits,momfbd}'
  endcase

  wbgfile = file_search(dir + '/*_'+string(scanno, format = '(i05)')+'_[0-9][0-9][0-9][0-9]' + extension $
                        , count = Nfiles)
  if Nfiles ne 1 then stop
  self -> extractstates, wbgfile, wbgstate

  ;; Get the alignment mode, projective or polywarp.
  alignment_model = red_align_model(wbgfile)

;  files = file_search(dir + '*_'+string(scanno, format = '(i05)')+'_*_[0-9][0-9][0-9][0-9]' + extension, count = Nfiles) ; wbdetector? ; ; ;
;  self -> extractstates, files, states

;  indx = where(strmatch(files, '*'+wbdetector+'_*_[0-9][0-9][0-9][0-9]'+extension), Nmatch)
;  wbgfile = files[indx[0]]
;  wbgstate = states[indx[0]]  

  
;  ;; Get the FOV based on all scans.
;  wfiles = file_search(dir + wbdetector+'*' + extension, count = Nwfiles) 
;  ;; The global WB file name should have no tuning info.
;  indx = where(~strmatch(wfiles,'*_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_[+-]*'), Nmatch)
;  if Nmatch eq 0 then stop
;  wbgfiles = wfiles[indx]

  
  hdr = red_readhead(wbgfile)
  im_dim = fxpar(hdr, 'NAXIS*')
  x0 = 0
  x1 = im_dim[0]-1
  y0 = 0
  y1 = im_dim[1]-1
  Nx = x1 - x0 + 1
  Ny = y1 - y0 + 1

  prefilter = strtrim(fxpar(hdr,'FILTER1'), 2)

;  wchdr0 = red_readhead(wbgfile)
  datestamp = strtrim(fxpar(hdr, 'STARTOBS'), 2)
  timestamp = (strsplit(datestamp, 'T', /extract))[1]
  


  ;; All other files for this scan
  pertuningfiles = file_search(dir + '/*_'+string(scanno, format = '(i05)') $
                               + '_*[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_[+-]*' $
                               + extension $
                               , count = Nfiles)
  self -> extractstates, pertuningfiles, pertuningstates

  if Nfiles eq 0 then begin
    red_message, 'No files for scan ' + strtrim(scanno, 2)
    return
  endif
  

  
  ;; Unique tuning states, sorted by wavelength
;  ufpi_states = red_uniquify(pertuningstates[where(~pertuningstates.is_wb)].fpi_state, indx = xindx)
;  Nwav = n_elements(xindx)
;  utunwavelength = pertuningstates[xindx].tun_wavelength
;  sortindx = sort(utunwavelength)
;  utunwavelength = utunwavelength[sortindx]
;  ufpi_states = ufpi_states[sortindx]
;  my_prefilters = pertuningstates[xindx].prefilter
;  wav = utunwavelength

  ;; Unique nb prefilters
  unbprefs = red_uniquify(pertuningstates[where(~pertuningstates.is_wb)].prefilter, indx = unbprefindx)
  Nnbprefs = n_elements(unbprefs)
;  unbprefsref = dblarr(Nnbprefs)
;  
;  
;  for inbpref = 0L, Nnbprefs-1 do begin
;    ;; This is the reference point of the fine tuning for this prefilter:
;    unbprefsref[inbpref] = double((strsplit(ufpi_states[0], '_', /extract))[0])
;  endfor                        ; inbpref
;  
;  unbprefsref *= 1e-10          ; [m]
  
  ;; Load prefilters
  if ~keyword_set(fitpref_time) then fitpref_time='_'

  for ipref = 0, Nnbprefs-1 do begin

    ;; Select per-tuning files for the selected scan and the current
    ;; prefilter.

    self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                                 , cam = nbtcamera $
                                 , sel = nbtindx, count = Nnbt $
                                 , prefilter = unbprefs[ipref]
    nbtstates = pertuningstates[nbtindx]
    nbtfiles = pertuningfiles[nbtindx]
    
    self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                                 , cam = nbrcamera $
                                 , sel = nbrindx, count = Nnbr $
                                 , prefilter = unbprefs[ipref]
    nbrstates = pertuningstates[nbrindx]
    nbrfiles = pertuningfiles[nbrindx]

    if instrument eq 'Chromis' then begin
      ;; For CHROMIS WB, we also filter for the fpi_states from the
      ;; NBT file selection. This is because filtering for the
      ;; prefilter does not work for CHROMIS, as selectfiles does not
      ;; check that NB filters are identical to the one given in the
      ;; prefilter keyword, only that it "matches" the WB filter. So
      ;; in this sense, 3634, 3969, and 3999 are all the same.
      self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                                   , cam = wbcamera $
                                   , sel = wbindx, count = Nwb $
                                   , fpi_state = red_uniquify(nbtstates.fpi_state)
    endif else begin
      ;; For CRISP and CRISP2 filtering for the prefilter should be
      ;; ok. And doing it the CHROMIS way didn't work for some reason.
      self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                                   , cam = wbcamera $
                                   , sel = wbindx, count = Nwb $
                                   , prefilter = unbprefs[ipref]
    endelse
    wbstates = pertuningstates[wbindx]
    wbfiles = pertuningfiles[wbindx]

    help, wbstates, nbrstates, nbtstates
    if n_elements(wbstates) ne n_elements(nbrstates) then stop

    ;; T camera

    pfile = self.out_dir + '/prefilter_fits/'+instrument+'-T_'+unbprefs[ipref]+fitpref_time+'prefilter.idlsave'
    if ~file_test(pfile) then begin
      red_message, 'Prefilter file not found: '+pfile
      return
    endif
    restore, pfile              ; Restores variable prf which is a struct
    prft = prf
    
    nbt_units = prft.units
;  nbt_prefilter_curve = prf.pref
;  nbt_prefilter_wav = prf.wav
;  nbt_prefilter_wb = prf.wbint
    
    ;; R camera
    
    pfile = self.out_dir + '/prefilter_fits/'+instrument+'-R_'+unbprefs[ipref]+fitpref_time+'prefilter.idlsave'
    if ~file_test(pfile) then begin
      red_message, 'Prefilter file not found: '+pfile
      return
    endif
    restore, pfile              ; Restores variable prf which is a struct
    prfr = prf
    
    nbr_units = prfr.units  
;  nbr_prefilter_curve = prf.pref
;  nbr_prefilter_wav = prf.wav
;  nbr_prefilter_wb = prf.wbint

    if nbr_units ne nbt_units then begin
      red_message, 'Units for the T and R cameras do not match. ' $
                   + 'Please rerun the prefilterfit step for these data.'
      retall
    endif
    units = nbr_units
    
    
    if ~keyword_set(nocavitymap) then begin
      
      ;; Read the original cavity map
      cfile = self.out_dir + 'flats/spectral_flats/' $
              + strjoin([nbtdetector $
                         , unbprefs[ipref] $
                         , 'fit_results.sav'] $
                        , '_')
      
      if ~file_test(cfile) then begin
        red_message, 'Error, calibration file not found -> '+cfile
        stop
      endif
      restore, cfile            ; The cavity map is in a struct called "fit". 
      dims = size(fit.pars, /dim)

      if dims[0] gt 1 then begin
        cmapt = reform(fit.pars[1,*,*]) ; Unit is [Angstrom]
        cmapt /= 10.                    ; Make it [nm]
        cmapt = -cmapt                  ; Change sign so lambda_correct = lambda + cmap
        fit = 0B                        ; Don't need the fit struct anymore.
        
        cfile = self.out_dir + 'flats/spectral_flats/' $
                + strjoin([nbrdetector $
                           , unbprefs[ipref] $
                           , 'fit_results.sav'] $
                          , '_')
        
        if ~file_test(cfile) then begin
          red_message, 'Error, calibration file not found -> '+cfile
          stop
        endif
        restore, cfile                  ; The cavity map is in a struct called "fit". 
        cmapr = reform(fit.pars[1,*,*]) ; Unit is [Angstrom]
        cmapr /= 10.                    ; Make it [nm]
        cmapr = -cmapr                  ; Change sign so lambda_correct = lambda + cmap
        fit = 0B                        ; Don't need the fit struct anymore.
        
        if keyword_set(remove_smallscale) then begin
          ;; If the small scale is already corrected, then include only the
          ;; low-resolution component in the metadata. The blurring kernel
          ;; should match how the low resolution component was removed when
          ;; making flats.
          npix = 30             ; Can we get this parameter from earlier headers?
          cpsf = red_get_psf(npix*2-1,npix*2-1,double(npix),double(npix))
          cpsf /= total(cpsf, /double)
          cmapr = red_convolve(temporary(cmapr), cpsf)
          cmap1r = cmapr
          cmapt = red_convolve(temporary(cmapt), cpsf)
          cmap1t = cmapt
        endif else begin
          ;; If the small scale is not already corrected, then we still want
          ;; to blur the cavity map slightly.
          npsf = round(fwhm * 7.)
          if((npsf/2)*2 eq npsf) then npsf += 1L
          psf = red_get_psf(npsf, npsf, fwhm, fwhm)
          psf /= total(psf, /double)
          ;; Leave the orignal cmaps alone, we might need them later.
          cmap1r = red_convolve(cmapr, psf)
          cmap1t = red_convolve(cmapt, psf)
        endelse
        
        ;; Apply geometrical transformation from the pinhole calibration to the cavity maps.
        cmap1t = red_apply_camera_alignment(cmap1t, alignment_model, instrument+'-T', amap = amapt)
        cmap1r = red_apply_camera_alignment(cmap1r, alignment_model, instrument+'-R', amap = amapr)

        ;; Clip to the selected FOV
        cmap1r = cmap1r[x0:x1,y0:y1]
        cmap1t = cmap1t[x0:x1,y0:y1]
        
        ;; At this point, the individual cavity maps should be corrected
        ;; for camera misalignments, so they should be aligned with
        ;; respect to the cavity errors on the etalons. So we can sum
        ;; them.
        cmap1 = (cmap1r + cmap1t) / 2.
      endif else undefine, cmap1
    endif else undefine, cmap1



    ;; First get the inverse modulation matrices, make them if
    ;; needed. They are returned in the (size and) orientation of
    ;; the momfbd output.
    self -> inverse_modmatrices, unbprefs[ipref], stokesdir $
       , camr = nbrcamera, immr = immr $
       , camt = nbtcamera, immt = immt $
       , no_ccdtabs = no_ccdtabs


;  ;; The global WB file name should have no tuning info.
;  indx = where(~strmatch(files,'*_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_[+-]*'), Nscans)
;  if Nscans eq 0 then stop
;  wbgfile = files[indx]
;  wbgstate = states[indx]

    hdr = red_readhead(wbgfile)
    red_fitspar_getdates, hdr $
                          , date_beg = date_beg $
                          , date_end = date_end $
                          , date_avg = date_avg $
                          , count_avg = hasdateavg $
                          , comment_avg = comment_avg
    if hasdateavg then begin
      date_avg_split = strsplit(date_avg, 'T', /extract, count = Nsplit)
      ddate = date_avg_split[0]
      if Nsplit gt 1 then ttime = date_avg_split[1] else undefine, ttime
    endif else undefine, ddate, ttime
    ;; Get pointing at center of FOV
    red_wcs_hpl_coords, red_time2double(ttime), metadata_pointing, time_pointing $
                        , hpln, hplt
    
    ;; Get the image scale from the header
    image_scale = float(fxpar(hdr, 'CDELT1'))
    
;  ;; Read wcs extension of wb file to get pointing info
;  fxbopen, wlun, wcfile, 'WCS-TAB', wbdr
;  ttype1 = fxpar(wbdr, 'TTYPE1')
;  fxbread, wlun, wwcs, ttype1
;  fxbclose, wlun
    ;; Get pointing info for the current scan directly
    ;;red_wcs_hpl_coords, t_array[0, *], metadata_pointing, time_pointing, hpln, hplt
;  Or skip having pointing info in the stokes cubes?

;  ;; Find the nb and wb per-tuning files by excluding the global WB image
;  self -> selectfiles, files = files, states = states $
;                       , cam = wbcamera, ustat = '' $
;                       , sel = wbgindx, count = Nscans $
;                       , complement = complement, Ncomplement = Ncomplement
;  ;; We have no special state (or absence of state) to identify
;  ;; the global WB images but we do know that their exposure times
;  ;; are much larger than the ones corresponding to the individual
;  ;; NB states.
;  wbindx = where(states.exposure gt mean(states.exposure)*1.5 $
;                 , Nscans, complement = complement, Ncomplement = Ncomplement)
;  
;  ;; All the per-tuning files and states
;  pertuningfiles = files[complement]
;  pertuningstates = states[complement]

;  ;; Per-tuning files, wb and nb, only for selected scans
;  self -> selectfiles, files = pertuningfiles, states = pertuningstates $
;                       , cam = wbcamera $
;                       , sel = wbindx, count = Nwb
;  wbstates = pertuningstates[wbindx]
;  wbfiles = pertuningfiles[wbindx]
;
;  self -> selectfiles, files = pertuningfiles, states = pertuningstates $
;                       , cam = nbtcamera $
;                       , sel = nbtindx, count = Nnbt
;  nbtstates = pertuningstates[nbtindx]
;  nbtfiles = pertuningfiles[nbtindx]
;  
;  self -> selectfiles, files = pertuningfiles, states = pertuningstates $
;                       , cam = nbrcamera $
;                       , sel = nbrindx, count = Nnbr
;  nbrstates = pertuningstates[nbrindx]
;  nbrfiles = pertuningfiles[nbrindx]
;
    ;; Unique tuning states, sorted by wavelength
    utunindx = uniq(nbrstates.fpi_state, sort(nbrstates.fpi_state))
    Ntuning = n_elements(utunindx)
    sortindx = sort(nbrstates[utunindx].tun_wavelength)
    utuning = nbrstates[utunindx[sortindx]].fpi_state

;    help, prft.wav, prft.pref, nbtstates[utunindx[sortindx]].tun_wavelength*1.d10
    if n_elements(prft.wav) eq 1 then begin
      ;; This is for CHROMIS 3999, Ca II continuum
      nbt_prefilter_curve = prft.pref
      nbr_prefilter_curve = prfr.pref
    endif else begin
      nbt_prefilter_curve = red_intepf(prft.wav, prft.pref, nbtstates[utunindx[sortindx]].tun_wavelength*1.d10)
      nbr_prefilter_curve = red_intepf(prfr.wav, prfr.pref, nbrstates[utunindx[sortindx]].tun_wavelength*1.d10)
    endelse
    
    nbt_rpref = 1.d0/nbt_prefilter_curve
    nbr_rpref = 1.d0/nbr_prefilter_curve

;    ;; Enough LC states to make a Stokes cube?
;    ulc = nbrstates[uniq(nbrstates.lc, sort(nbrstates.lc))].lc
;    Nlc = n_elements(ulc)
;    if Nlc eq 1 then begin      ; or Nlc lt 4
;      print, 'Just a single LC state for scan ' + strtrim(scanno, 2)
;      return
;    endif

    if Nnbt ne Nnbr then stop
    
    Nstokes = 4

    ;; Define the Stokes file names needed for this scan
    these_snames = strarr(Ntuning)    
    for ituning = 0, Ntuning-1 do begin
      these_snames[ituning] = stokesdir $
                              + strjoin(['stokesIQUV' $
                                         , string(scanno, format = '(i05)') $
                                         , prefilter $ ; The WB prefilter
                                         , utuning[ituning] $
                                        ], '_') + '.fits' 
    endfor                      ; ituning

    red_append, snames, these_snames
    
    if keyword_set(redemodulate) then begin
      ;; Delete any stokesIQUV*.fits files that are to be remade.
      red_message, 'Will delete the following Stokes cubes if they exist:'
      print, these_snames, format = '(a0)'
      file_delete, these_snames, /ALLOW_NONEXISTENT
    endif


    ;; Make Stokes cubes for each state (if not done already) 
    todoindx = where(~file_test(these_snames), Ntodo)
    if Ntodo gt 0 then begin
      red_message, 'Will have to make '+strtrim(Ntodo, 2) + ' Stokes cubes for scan '+strtrim(scanno, 2)+'.'

      ;; Get the FOV in the momfbd files.
;    mr = momfbd_read(wbgfile, /names) ; Use /names to avoid reading the data parts
;    mrX01Y01 = mr.roi + mr.margin * [1, -1, 1, -1]
      
      ;; Read the global WB file for this scan.
      wbg = red_readdata(wbgfile)
      
      for ituning = 0, Ntuning-1 do begin

        if strmid(utuning[ituning],0,4) ne unbprefs[ipref] then continue ; Done with another ipref.

        if file_test(these_snames[ituning]) then continue ; Already done

        red_progressbar, ituning, Ntuning, /predict  $
                         , 'Making '+file_basename(these_snames[ituning])

        if n_elements(wbfiles) eq 1 then begin
          ;; selectfiles does not work for a single file? Do this for
          ;; now:
          these_wbindx = [0]
          Nthesewb = 1
          these_nbtindx = [0]
          Nthesenb1 = 1
          these_nbrindx = [0]
          Nthesenbr = 1
        endif else begin
          self -> selectfiles, files = wbfiles, states = wbstates $
                                       , sel = these_wbindx, count = Nthesewb $
                                       , scan = scanno $
                                       , fpi_states = utuning[ituning]
          
          self -> selectfiles, files = nbtfiles, states = nbtstates $
                                       , sel = these_nbtindx, count = Nthesenbt $
                                       , scan = scanno $
                                       , fpi_states = utuning[ituning]
          
          self -> selectfiles, files = nbrfiles, states = nbrstates $
                                       , sel = these_nbrindx, count = Nthesenbr $
                                       , scan = scanno $
                                       , fpi_states = utuning[ituning]
        endelse

        wcs = {wave:dblarr(2,2)   $ ; WCS for this Stokes cube.
               , hplt:dblarr(2,2) $
               , hpln:dblarr(2,2) $
               , time:dblarr(2,2) $
              }
        
        wcs.hpln[0, 0] = hpln - double(image_scale) * (Nx-1)/2.d
        wcs.hpln[1, 0] = hpln + double(image_scale) * (Nx-1)/2.d
        wcs.hpln[0, 1] = hpln - double(image_scale) * (Nx-1)/2.d
        wcs.hpln[1, 1] = hpln + double(image_scale) * (Nx-1)/2.d
        
        wcs.hplt[0, 0] = hplt - double(image_scale) * (Ny-1)/2.d
        wcs.hplt[1, 0] = hplt - double(image_scale) * (Ny-1)/2.d
        wcs.hplt[0, 1] = hplt + double(image_scale) * (Ny-1)/2.d
        wcs.hplt[1, 1] = hplt + double(image_scale) * (Ny-1)/2.d

        wcs.wave = nbtstates[these_nbrindx[0]].tun_wavelength*1d9
        ;; wcs.time = ; Set by demodulate

;        print, ipref, ituning
;        print, snames[ituning]
;        stop

        help, these_nbrindx,  nbrstates, nbtstates, wbstates, nbrstates[these_nbrindx],  nbtstates[these_nbtindx], wbstates[these_wbindx]

        if n_elements(wbstates[these_wbindx]) lt 4 then stop
        
        self -> demodulate, these_snames[ituning], immr, immt $
           , smooth_by_subfield = smooth_by_subfield $ 
           , smooth_by_kernel = smooth_by_kernel $ 
           , clips = clips $
           , cmap = cmap1 $
           , nbrfac = nbr_rpref[ituning] $
           , nbrstates = nbrstates[these_nbrindx] $
           , nbtfac = nbt_rpref[ituning] $
           , nbtstates = nbtstates[these_nbtindx] $
           , noremove_periodic = noremove_periodic $
           , notelmat = notelmat $
           , overwrite = redemodulate $
           , tiles = tiles $
           , units = units $
           , wbg = wbg $
           , wcs = wcs $
           , wbstates = wbstates[these_wbindx] $
           , nthreads = nthreads $
           , nearest = nearest
        
      endfor                    ; ituning

    endif else begin

      print
      red_message, 'No Stokes cubes need to be made for scan '+strtrim(scanno, 2)
      print
      
    endelse
    
  endfor                        ; ipref
  
end
