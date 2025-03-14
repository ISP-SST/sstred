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
                            , redemodulate = redemodulate $
;                              , smooth = smooth $
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
  
  ;; How to smooth the modulation matrices.
  if n_elements(smooth) eq 0 then begin
    ;; Smooth with Gaussian kernel by default
    smooth_by_kernel = 5        ; Default width
    smooth_by_subfield = 0
  endif else begin
    ;; The smooth keyword can either be a number, in which case that
    ;; is the kernel width, or the string "momfbd", in which case we
    ;; do smoothing by subfield using the MOMFBD-estimated PSFs.
    if size(smooth, /tname) eq 'STRING' then begin
      if strlowcase(smooth) eq 'momfbd' then begin
        ;; If the string "momfbd" (or "MOMFBD"), we will smooth by
        ;; subfield. 
        smooth_by_subfield = 1
      endif else begin
        ;; Any string except "momfbd" will result in no smoothing. 
        smooth_by_subfield = 0
        smooth_by_kernel = 0
      endelse
    endif else begin
      ;; Not a string, then hopefully a number
      smooth_by_subfield = 0
      smooth_by_kernel = smooth
    endelse
  endelse
  
  
  ;; Default keywords
  if n_elements(cmap_fwhm) eq 0 then fwhm = 7.0
  if n_elements(tiles) eq 0 or n_elements(clips) eq 0 then begin
    tiles = [12, 16, 32, 64, 72]
    clips = [12, 8,  4,  2 , 1]
  endif


  ;; Camera/detector identification
  self->getdetectors
  wbindx     = where(strmatch(*self.cameras,'Crisp-W'))
  wbcamera   = (*self.cameras)[wbindx[0]]
  wbdetector = (*self.detectors)[wbindx[0]]
  nbtindx     = where(strmatch(*self.cameras,'Crisp-T')) 
  nbtcamera   = (*self.cameras)[nbtindx[0]]
  nbtdetector = (*self.detectors)[nbtindx[0]]
  nbrindx     = where(strmatch(*self.cameras,'Crisp-R')) 
  nbrcamera   = (*self.cameras)[nbrindx[0]]
  nbrdetector = (*self.detectors)[nbrindx[0]]

  ;; We currently do correct for the small scale cavity map in CRISP
  ;; data. (We should get this from earlier meta data!)
  remove_smallscale = 1


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

  wchdr0 = red_readhead(wbgfile)
  datestamp = strtrim(fxpar(wchdr0, 'STARTOBS'), 2)
  timestamp = (strsplit(datestamp, 'T', /extract))[1]
  


  ;; All other files for this scan
  pertuningfiles = file_search(dir + '/*_'+string(scanno, format = '(i05)') $
                               + '_*[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_[+-]*' $
                               + extension $
                               , count = Nfiles)
  self -> extractstates, pertuningfiles, pertuningstates

  ;; Per-tuning files, wb and nb, only for selected scans
  self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                               , cam = wbcamera $
                               , sel = wbindx, count = Nwb
  wbstates = pertuningstates[wbindx]
  wbfiles = pertuningfiles[wbindx]

  self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                               , cam = nbtcamera $
                               , sel = nbtindx, count = Nnbt
  nbtstates = pertuningstates[nbtindx]
  nbtfiles = pertuningfiles[nbtindx]
  
  self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                               , cam = nbrcamera $
                               , sel = nbrindx, count = Nnbr
  nbrstates = pertuningstates[nbrindx]
  nbrfiles = pertuningfiles[nbrindx]

;  stop
  
  ;; Load prefilters
  if ~keyword_set(fitpref_time) then fitpref_time='_'
  ;; Crisp-T

  pfile = self.out_dir + '/prefilter_fits/Crisp-T_'+prefilter+fitpref_time+'prefilter.idlsave'
  if ~file_test(pfile) then begin
    red_message, 'Prefilter file not found: '+pfile
    return
  endif
  restore, pfile                ; Restores variable prf which is a struct
  prft = prf
  
  nbt_units = prft.units
;  nbt_prefilter_curve = prf.pref
;  nbt_prefilter_wav = prf.wav
;  nbt_prefilter_wb = prf.wbint
  
  ;; Crisp-R

  pfile = self.out_dir + '/prefilter_fits/Crisp-R_'+prefilter+fitpref_time+'prefilter.idlsave'
  if ~file_test(pfile) then begin
    red_message, 'Prefilter file not found: '+pfile
    return
  endif
  restore, pfile                ; Restores variable prf which is a struct
  prfr = prf
  
  nbr_units = prfr.units  
;  nbr_prefilter_curve = prf.pref
;  nbr_prefilter_wav = prf.wav
;  nbr_prefilter_wb = prf.wbint
  



  
  if nbr_units ne nbt_units then begin
    red_message, 'Units for Crisp-T and Crisp-R do not match. ' $
                 + 'Please rerun the prefilterfit step for these data.'
    retall
  endif
  units = nbr_units


  if ~keyword_set(nocavitymap) then begin

    ;; Read the original cavity map
    cfile = self.out_dir + 'flats/spectral_flats/' $
            + strjoin([nbtdetector $
                       , prefilter $
                       , 'fit_results.sav'] $
                      , '_')

    if ~file_test(cfile) then begin
      red_message, 'Error, calibration file not found -> '+cfile
      stop
    endif
    restore, cfile                  ; The cavity map is in a struct called "fit". 
    cmapt = reform(fit.pars[1,*,*]) ; Unit is [Angstrom]
    cmapt /= 10.                    ; Make it [nm]
    cmapt = -cmapt                  ; Change sign so lambda_correct = lambda + cmap
    fit = 0B                        ; Don't need the fit struct anymore.
    
    cfile = self.out_dir + 'flats/spectral_flats/' $
            + strjoin([nbrdetector $
                       , prefilter $
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
      npix = 30                 ; Can we get this parameter from earlier headers?
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
    cmap1t = red_apply_camera_alignment(cmap1t, alignment_model, 'Crisp-T', amap = amapt)
    cmap1r = red_apply_camera_alignment(cmap1r, alignment_model, 'Crisp-R', amap = amapr)

    ;; Clip to the selected FOV
    cmap1r = cmap1r[x0:x1,y0:y1]
    cmap1t = cmap1t[x0:x1,y0:y1]
    
    ;; At this point, the individual cavity maps should be corrected
    ;; for camera misalignments, so they should be aligned with
    ;; respect to the cavity errors on the etalons. So we can sum
    ;; them.
    cmap1 = (cmap1r + cmap1t) / 2.
    
  endif


  ;; Store Stokes cubes in separate directories for different smooth
  ;; options:
  stokesdir = dir + '/stokes_sbs'+strtrim(smooth_by_subfield,2) $
              + '_sbk'+strtrim(smooth_by_kernel,2)+'/'

  file_mkdir, stokesdir

  ;; First get the inverse modulation matrices, make them if
  ;; needed. They are returned in the (size and) orientation of
  ;; the momfbd output.
  self -> inverse_modmatrices, prefilter, stokesdir $
     , camr = nbrcamera, immr = immr $
     , camt = nbtcamera, immt = immt $
     , no_ccdtabs = no_ccdtabs

  if Nfiles eq 0 then begin
    red_message, 'No files for scan ' + strtrim(scanno, 2)
    return
  endif
  

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

  
  nbt_prefilter_curve = red_intepf(prft.wav, prft.pref, nbtstates[utunindx[sortindx]].tun_wavelength*1.d10)
  nbr_prefilter_curve = red_intepf(prfr.wav, prfr.pref, nbrstates[utunindx[sortindx]].tun_wavelength*1.d10)

  nbt_rpref = 1.d0/nbt_prefilter_curve
  nbr_rpref = 1.d0/nbr_prefilter_curve

  ;; Enough LC states to make a Stokes cube?
  ulc = nbrstates[uniq(nbrstates.lc, sort(nbrstates.lc))].lc
  Nlc = n_elements(ulc)
  if Nlc eq 1 then begin        ; or Nlc lt 4
    print, 'Just a single LC state for scan ' + strtrim(scanno, 2)
    return
  endif

  if Nnbt ne Nnbr then stop

  Nstokes = 4

  ;; Define the Stokes file names needed for this scan
  snames = strarr(Ntuning)    
  for ituning = 0, Ntuning-1 do begin
    snames[ituning] = stokesdir $
                      + strjoin(['stokesIQUV' $
                                 , string(scanno, format = '(i05)') $
                                 , prefilter $
                                 , utuning[ituning] $
                                ], '_') + '.fits' 
  endfor                        ; ituning

  
  if keyword_set(redemodulate) then begin
    ;; Delete any stokesIQUV*.fits files that are to be remade.
    red_message, 'Will delete the following Stokes cubes if they exist:'
    print, snames, format = '(a0)'
    file_delete, snames, /ALLOW_NONEXISTENT
  endif


  ;; Make Stokes cubes for each state (if not done already) 
  todoindx = where(~file_test(snames), Ntodo)
  if Ntodo gt 0 then begin
    red_message, 'Will have to make '+strtrim(Ntodo, 2) + ' Stokes cubes for scan '+strtrim(scanno, 2)+'.'

    ;; Get the FOV in the momfbd files.
;    mr = momfbd_read(wbgfile, /names) ; Use /names to avoid reading the data parts
;    mrX01Y01 = mr.roi + mr.margin * [1, -1, 1, -1]
    
    ;; Read the global WB file for this scan.
    wbg = red_readdata(wbgfile)
    
    for ituning = 0, Ntuning-1 do begin
      
      if file_test(snames[ituning]) then continue

      red_progressbar, ituning, Ntuning, /predict  $
                       , 'Making '+file_basename(snames[ituning])

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

      self -> demodulate, snames[ituning], immr, immt $
         , smooth_by_subfield = smooth_by_subfield $ 
         , smooth_by_kernel = smooth_by_kernel $ 
         , clips = clips $
         , cmap = cmap1 $
         , nbrfac = nbr_rpref[ituning] $
         , nbrstates = nbrstates[these_nbrindx] $
         , nbtfac = nbt_rpref[ituning] $
         , nbtstates = nbtstates[these_nbtindx] $
         , noremove_periodic = noremove_periodic $
         , overwrite = redemodulate $
         , tiles = tiles $
         , units = units $
         , wbg = wbg $
         , wcs = wcs $
         , wbstates = wbstates[these_wbindx] $
         , nthreads = nthreads $
         , nearest = nearest

    endfor                      ; ituning

  endif else begin

    print
    red_message, 'No Stokes cubes need to be made for scan '+strtrim(scanno, 2)
    print
    
  endelse
  
end
