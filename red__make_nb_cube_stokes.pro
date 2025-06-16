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
                              , clips = clips  $
                              , cmap_fwhm = cmap_fwhm $
                              , fitpref_time = fitpref_time $
                              , integer = integer $
                              , intensitycorrmethod = intensitycorrmethod $ 
                              , nearest = nearest $
                              , noaligncont = noaligncont $  
                              , nocavitymap = nocavitymap $
                              , nocrosstalk = nocrosstalk $
                              , noflipping = noflipping $
                              , nomissing_nans = nomissing_nans $
                              , noremove_periodic = noremove_periodic $
                              , nostretch = nostretch $
                              , notimecor = notimecor $
                              , nthreads = nthreads $
                              , odir = odir $
                              , overwrite = overwrite $
                              , redemodulate = redemodulate $
                              , tiles = tiles $
                              , wbsave = wbsave

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  ;; Deprecated keyword:
  if n_elements(notimecor) gt 0 then begin
    red_message, 'Keyword notimecor is deprecated. Use intensitycorrmethod="none" instead.'
    return
  endif

  instrument = ((typename(self)).tolower())
  
  case instrument of
    'crisp' : begin
      ;; We currently do correct for the small scale cavity map in CRISP
      ;; data. (We should get this from earlier meta data!)
      remove_smallscale = 1
    end
    'chromis' : begin
      remove_smallscale = 0
    end
  endcase

;  if(keyword_set(nearest)) then near = 1 else near = 0
  
  ;; Make prpara
  red_make_prpara, prpara, clips         
  red_make_prpara, prpara, integer
  red_make_prpara, prpara, intensitycorrmethod  
  red_make_prpara, prpara, cmap_fwhm
  red_make_prpara, prpara, nearest
  red_make_prpara, prpara, nocavitymap 
  red_make_prpara, prpara, nocrosstalk 
  red_make_prpara, prpara, nomissing_nans
  red_make_prpara, prpara, nostretch 
  red_make_prpara, prpara, np           
  red_make_prpara, prpara, redemodulate
  red_make_prpara, prpara, tiles        
  red_make_prpara, prpara, wcfile

  if n_elements(nthreads) eq 0 then nthreads = 1 ; Default single thread
  
  ;; Default keywords
  if n_elements(cmap_fwhm) eq 0 then fwhm = 7.0
  if n_elements(tiles) eq 0 or n_elements(clips) eq 0 then begin
    tiles = [8, 16, 32, 64, 128]
    clips = [8, 4,  2,  1,  1  ]
  endif

  if keyword_set(redemodulate) then overwrite = 1

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

  ;; Default for wb cubes without direction parameter
  if n_elements(direction) eq 0 then direction = 0
  
  ;; CRISP has a common prefilter for WB and NB, for CHROMIS this is
  ;; the WB filter.
  prefilter = fxpar(wchead,'FILTER1')

  ;; Read wcs extension of wb file to get pointing info
  fxbopen, wlun, wcfile, 'WCS-TAB', wbdr
  ttype1 = fxpar(wbdr, 'TTYPE1')
  fxbread, wlun, wwcs, ttype1
  fxbclose, wlun

  ;; Prepare for making output file names
  if n_elements(odir) eq 0 then odir = self.out_dir + '/cubes_nb/' 

  ofile = red_strreplace(file_basename(wcfile), 'wb', 'nb')
  ofile = red_strreplace(ofile, 'corrected', 'stokes_corrected')
  filename = odir+ofile

  ;; Already done?
  if file_test(filename) then begin
    if keyword_set(overwrite) then begin
      red_message, 'Overwriting existing data cube: ' + filename
    endif else begin
      red_message, 'This data cube exists already: ' + filename
      return
    endelse
  endif

  file_mkdir, odir


  
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

  case self.filetype of
    'ANA': extension = '.f0'
    'MOMFBD': extension = '.momfbd'
    'FITS': extension = '.fits'
    'MIXED' : extension = '.{fits,momfbd}'
  endcase
  
  ;;extension = (strsplit(wbgfiles[0],'.',/extract))[-1]  
  for jj=0,n_elements(wbgfiles)-1 do begin
    search_dir = file_dirname(wbgfiles[jj])+'/'
    srch = '*_' + string(wbgstates[jj].scannumber, format = '(I05)')+'_*'
    ff = file_search(search_dir + srch + extension) 
    red_append,files,ff
  endfor
  Nfiles = n_elements(files)  

  
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

  ;; Unique nb prefilters
  unbprefs = red_uniquify(pertuningstates[where(~pertuningstates.is_wb)].prefilter, indx = unbprefindx)
  Nnbprefs = n_elements(unbprefs)

  ;; Get the scan selection from wfiles (from the sav file)
  self -> extractstates, wbgfiles, wbgstates
  uscans = wbgstates.scannumber
  Nscans = n_elements(uscans)

  ;; What prefilter fits directory to use?
  if ~keyword_set(fitpref_time) then begin
    fitpref_t='_'
    dt = strtrim(fxpar(wchdr0, 'DATE-AVG'), 2)
    avg_ts = (strsplit(dt, 'T', /extract))[1]
    avg_time = red_time2double(avg_ts)
    pfls = file_search(self.out_dir + '/prefilter_fits/*-T_'+unbprefs+ $
                       '_[0-9][0-9]:[0-9][0-9]:[0-9][0-9]*save', count=Npfls)
;    if Npfls gt 0 then begin
    if Npfls eq Nnbprefs then begin
      tt = dblarr(Npfls)
      ts = strarr(Npfls)
      for ii=0,Npfls-1 do begin
        ts[ii] = (strsplit(file_basename(pfls[ii]),'_',/extract))[2]
        tt[ii] = abs(red_time2double(ts[ii]) - avg_time)
      endfor
      mn = min(tt,jj)
      fitpref_t = '_'+ts[jj]+'_'
    endif
  endif else fitpref_t = '_'+fitpref_time+'_' 

  
  ;; Define the Stokes file names needed for this nb cube, make the
  ;; Stokes cubes if needed.

  for iscan = 0, Nscans-1 do begin
    
    self -> make_stokes_cubes, file_dirname(wbgfiles[iscan]), uscans[iscan] $
       , clips = clips $
       , cmap_fwhm = cmap_fwhm $
       , /nocavitymap $         ; Cavity maps in Stokes cubes aren't used for anything
       , /notelmat $            ; Until Pit has had time to measure them
       , noremove_periodic = noremove_periodic $
       , redemodulate = redemodulate $
       , snames = these_snames $
       , stokesdir = stokesdir $
       , tiles = tiles $
       , nearest = nearest $
       , nthreads = nthreads $
       , fitpref_time = fitpref_t

    if n_elements(snames) eq 0 then begin
      Ntuning = n_elements(these_snames)
      snames = strarr(Nscans, Ntuning)
    endif
    
    snames[iscan, *] = these_snames
    
  endfor                        ; iscan

  self -> extractstates, snames, sstates

  ;; Load prefilters
  self -> get_prefilterfit, unbprefs $
     , units = units_x $
     , preflter_curve = prefilter_curve_x $
     , wave_shifts = wave_shifts_x
  
  for ipref = 0, Nnbprefs-1 do begin

    red_progressbar, ipref, Nnbprefs, 'Check units for '+unbprefs[ipref]
    
    ;; T camera

    pfile = self.out_dir + '/prefilter_fits/'+instrument.capwords()+'-T_'+unbprefs[ipref] $
            +fitpref_t+'prefilter.idlsave'
    if ~file_test(pfile) then begin
      red_message, 'Prefilter file not found: '+pfile
      return
    endif
    restore, pfile              ; Restores variable prf which is a struct
    
    if n_elements(units) eq 0 then begin
      units = prf.units
    endif else begin
      if prf.units ne units then begin
        red_message, 'Units do not match: ' + strtrim(units, 2) + ' '  + strtrim(prf.units, 2)
        red_message, 'Please rerun the prefilterfit step for these data.'
        retall
      endif
    endelse
    
;  nbt_prefilter_curve = prf.pref
;  nbt_prefilter_wav = prf.wav
;  nbt_prefilter_wb = prf.wbint
    
    ;; R camera
    
    pfile = self.out_dir + '/prefilter_fits/'+instrument.capwords()+'-R_'+unbprefs[ipref] $
            +fitpref_t+'prefilter.idlsave'
    if ~file_test(pfile) then begin
      red_message, 'Prefilter file not found: '+pfile
      return
    endif
    restore, pfile              ; Restores variable prf which is a struct
    
    if prf.units ne units then begin
      red_message, 'Units do not match: ' + strtrim(units, 2) + ' '  + strtrim(prf.units, 2)
      red_message, 'Please rerun the prefilterfit step for these data.'
      retall
    endif
    
;  nbr_prefilter_curve = prf.pref
;  nbr_prefilter_wav = prf.wav
;  nbr_prefilter_wb = prf.wbint

  endfor                        ; ipref

  
  if ~keyword_set(nocavitymap) then begin

    ;; There should be one cavity map per prefilter. Read them!

    cavitymaps_pref = fltarr(origNx, origNy, Nnbprefs)
    
    for ipref = 0, Nnbprefs-1 do begin

      red_progressbar, ipref, Nnbprefs, 'Read cavity maps for '+unbprefs[ipref]
      
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
        cmap1t = red_apply_camera_alignment(cmap1t, alignment_model, instrument.capwords()+'-T', amap = amapt)
        cmap1r = red_apply_camera_alignment(cmap1r, alignment_model, instrument.capwords()+'-R', amap = amapr)

        ;; Clip to the selected FOV
        cmap1r = cmap1r[x0:x1,y0:y1]
        cmap1t = cmap1t[x0:x1,y0:y1]

        ;; At this point, the individual cavity maps should be corrected
        ;; for camera misalignments, so they should be aligned with
        ;; respect to the cavity errors on the etalons. So we can sum
        ;; them.
        cmap1 = (cmap1r + cmap1t) / 2.

        ;; Get the orientation right.
        cmap1 = red_rotate(cmap1, direction)

;    Are we doing things in the right order?
;        
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

        
        cavitymaps_pref[*, *, ipref] = cmap1
        
      endif                     ;else undefine, cmap1
    endfor                      ; ipref
  endif                         ;else undefine, cmap1

  
  
  ;; Per-tuning files, wb and nb, only for selected scans
  self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                               , scan = uscans $
                               , cam = wbcamera $
                               , sel = wbindx, count = Nwb
  wbstates = pertuningstates[wbindx]
  wbfiles = pertuningfiles[wbindx]

  self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                               , scan = uscans $
                               , cam = nbtcamera $
                               , sel = nbtindx, count = Nnbt
  nbtstates = pertuningstates[nbtindx]
  nbtfiles = pertuningfiles[nbtindx]
  
  self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                               , scan = uscans $
                               , cam = nbrcamera $
                               , sel = nbrindx, count = Nnbr
  nbrstates = pertuningstates[nbrindx]
  nbrfiles = pertuningfiles[nbrindx]

;  for ipref = 0, Nnbprefs-1 do begin
;
;    ;; Select per-tuning files for the selected scan and the current
;    ;; prefilter.
;
;    self -> selectfiles, files = pertuningfiles, states = pertuningstates $
;                                 , cam = nbtcamera $
;                                 , sel = nbtindx, count = Nnbt $
;                                 , prefilter = unbprefs[ipref]
;    nbtstates = pertuningstates[nbtindx]
;    nbtfiles = pertuningfiles[nbtindx]
;    
;    self -> selectfiles, files = pertuningfiles, states = pertuningstates $
;                                 , cam = nbrcamera $
;                                 , sel = nbrindx, count = Nnbr $
;                                 , prefilter = unbprefs[ipref]
;    nbrstates = pertuningstates[nbrindx]
;    nbrfiles = pertuningfiles[nbrindx]
;
;    ;; For WB, we also filter for the fpi_states from the NBT file
;    ;; selection. This is because filtering for the prefilter does not
;    ;; work for CHROMIS, as selectfiles does not check that NB filters
;    ;; are identical to the one given in the prefilter keyword, only
;    ;; that it "matches" the WB filter. So in this sense, 3634, 3969,
;    ;; and 3999 are all the same.
;    self -> selectfiles, files = pertuningfiles, states = pertuningstates $
;                                 , cam = wbcamera $
;                                 , sel = wbindx, count = Nwb $
;                                 , fpi_state = red_uniquify(nbtstates.fpi_state)
;    wbstates = pertuningstates[wbindx]
;    wbfiles = pertuningfiles[wbindx]
;
;  endfor                        ; ipref
;
;  stop  
  
    
  ;; Unique tuning states, sorted by wavelength
  utunindx = uniq(sstates.fpi_state, sort(sstates.fpi_state))
  Ntuning = n_elements(utunindx)
  sortindx = sort(sstates[utunindx].tun_wavelength)
  utuning = sstates[utunindx[sortindx]].fpi_state
  
  Nstokes = 4
  
  if Nnbt ne Nnbr then stop

;  ;; Do WB correction?
;  wbcor = Nwb eq Nnbt and ~keyword_set(nostretch)


  ;; Load WB image and define the image border
;  tmp = red_readdata(wbgfiles[0])

  ;; If multiple directories, the fov_mask should be the same. Or we
  ;; have to think of something.
  spl = strsplit(wbgfiles[0],'/',/extract)
  cw = where(strmatch(spl,'*cfg*'))
  cfg_dir=strjoin(spl[0:cw],'/')
  if self.filetype eq 'MOMFBD' then begin
    mr = momfbd_read(wbgfiles[0],/nam)
  endif else begin              ; get cropping from cfg file    
    cfg_file = cfg_dir+'/'+'momfbd_reduc_'+wbgstates[0].prefilter+'_'+$
               string(wbgstates[0].scannumber,format='(I05)')+'.cfg'
    cfg = redux_readcfg(cfg_file)
    num_points = long(redux_cfggetkeyword(cfg, 'NUM_POINTS'))
    margin = num_points/8
    sim_xy = redux_cfggetkeyword(cfg, 'SIM_XY', count = cnt)
    if cnt gt 0 then begin
      sim_xy = rdx_str2ints(sim_xy)
      indx = indgen(n_elements(sim_xy)/2)*2
      indy = indx+1
      sim_x = sim_xy[indx]
      sim_y = sim_xy[indy]   
    endif else begin
      sim_x = rdx_str2ints(redux_cfggetkeyword(cfg, 'SIM_X'))
      sim_y = rdx_str2ints(redux_cfggetkeyword(cfg, 'SIM_Y'))
    endelse
    xx0 = min(sim_x) + margin - num_points/2 
    xx1 = max(sim_x) - margin + num_points/2 - 1
    yy0 = min(sim_y) + margin - num_points/2 
    yy1 = max(sim_y) - margin + num_points/2 - 1      
  endelse
  
  if file_test(cfg_dir+'/fov_mask.fits') then begin
    fov_mask = readfits(cfg_dir+'/fov_mask.fits')
    if self.filetype eq 'MOMFBD' then $
       fov_mask = red_crop_as_momfbd(fov_mask, mr) $
    else $
       fov_mask = fov_mask[xx0:xx1,yy0:yy1]
    fov_mask = red_rotate(fov_mask, direction)
  endif
  
  ;; Create cubes for science data and scan-adapted cavity maps.
  cavitymaps = fltarr(Nx, Ny, 1, 1, Nscans)

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
;                                        , *self.cameras,instrument.capwords()+'-T', amap = amapt)
;    cmap1r = red_apply_camera_alignment(cmap1r, alignment_model $
;                                        , *self.cameras,instrument.capwords()+'-R', amap = amapr)
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
  
  ;; Make FITS header for the NB cube
  hdr = wchead                                                ; Start with the WB cube header
  red_headerinfo_deletestep, hdr, /all                        ; Remove make_wb_cube steps 
  red_fitsdelkeyword, hdr, 'VAR_KEYS'                         ; Variable keywords copied later

  red_fitsdelkeyword, hdr, 'STATE'                  ; Not a single state for cube 
  red_fitsdelkeyword, hdr, 'CHECKSUM'               ; Checksum for WB cube
  red_fitsdelkeyword, hdr, 'DATASUM'                  ; Datasum for WB cube
  dindx = where(strmid(hdr, 0, 4) eq 'DATA', Ndata)   ; DATA statistics keywords
  for idata = Ndata-1, 0, -1 do begin
    keyword = strtrim(strmid(hdr[dindx[idata]], 0, 8), 2)
    red_fitsdelkeyword, hdr, keyword
  endfor                        ; idata
  
  red_fitsaddkeyword, hdr, 'BITPIX', -32 ; Floats

  
  anchor = 'DATE'

  ;; Add some keywords
  red_fitsaddkeyword, anchor = anchor, hdr, 'OBS_HDU', 1
  
  ;; Add info to headers
  red_fitsaddkeyword, anchor = anchor, hdr, 'BUNIT', units, 'Units in array'
  red_fitsaddkeyword, anchor = anchor, hdr, 'BTYPE', 'Intensity', 'Type of data in array'

  ;; Copy some info from a Stokes cube header
  shdr = headfits(snames[0])
  red_fitscopykeyword, anchor = anchor, hdr, 'CAMERA',   shdr
  red_fitscopykeyword, anchor = anchor, hdr, 'DETECTOR', shdr
  red_fitscopykeyword, anchor = anchor, hdr, 'DETGAIN',  shdr
  red_fitscopykeyword, anchor = anchor, hdr, 'DETOFFS',  shdr
  red_fitscopykeyword, anchor = anchor, hdr, 'DETMODEL', shdr
  red_fitscopykeyword, anchor = anchor, hdr, 'DETFIRM',  shdr

  ;; Note that exposure times can be different for different tunings.
  ;; Both due to single-frame exposures and due to numbers of frames.
  ;; XPOSURE and TEXPOSUR. 

  ;; WB and NB data come from different cameras.
;  red_fitsaddkeyword, hdr, 'CAMERA', nbtcamera + ',' + nbrcamera
  ;; Get DETGAIN, DETOFFS, DETMODEL, DETFIRM from .fitsheader file,
  ;; i.e., red_readhead(.momfbd file). This has to be handled
  ;; differently for polarimetry data because each Stokes image is a
  ;; mix of two detectors.
;  mhdr = red_readhead(nbrfiles[0]) ; Header of momfbd output file
;  mhdt = red_readhead(nbtfiles[0]) ; Header of momfbd output file
;  red_fitsaddkeyword, hdr, 'DETECTOR', strtrim(fxpar(mhdr,'DETECTOR'), 2) $
;                      + ',' + strtrim(fxpar(mhdt,'DETECTOR'), 2)
;  red_fitsaddkeyword, hdr, 'DETGAIN',  fxpar(mhdr,'DETGAIN') + ',' + fxpar(mhdt,'DETGAIN')
;  red_fitsaddkeyword, hdr, 'DETOFFS',  fxpar(mhdr,'DETOFFS') + ',' + fxpar(mhdt,'DETOFFS')
;  red_fitsaddkeyword, hdr, 'DETMODEL', fxpar(mhdr,'DETMODEL') + ',' + fxpar(mhdt,'DETMODEL')
;  red_fitsaddkeyword, hdr, 'DETFIRM',  fxpar(mhdr,'DETFIRM') + ',' + fxpar(mhdt,'DETFIRM')

  ;; Initialize fits file, set up for writing the data part.
  dims = [Nx, Ny, Ntuning, Nstokes, Nscans] 
  self -> fitscube_initialize, filename, hdr, lun, fileassoc, dims 

  if keyword_set(wbsave) then begin
    wbfilename = red_strreplace(filename, 'nb_', 'wbalign_')
    wbdims = [Nx, Ny, Ntuning, 1, Nscans] 
    wbhdr = hdr
    self -> fitscube_initialize, wbfilename, wbhdr, wblun, wbfileassoc, wbdims 
  endif
  
  ;; Observations metadata variables
  tbeg_array     = dblarr(Ntuning, Nscans) ; Time beginning for state
  tend_array     = dblarr(Ntuning, Nscans) ; Time end for state
  tavg_array     = dblarr(Ntuning, Nscans) ; Time average for state
  date_beg_array = strarr(Ntuning, Nscans) ; DATE-BEG for state
  date_end_array = strarr(Ntuning, Nscans) ; DATE-END for state
  date_avg_array = strarr(Ntuning, Nscans) ; DATE-AVG for state
  exp_array      = fltarr(Ntuning, Nscans) ; Total exposure time
  sexp_array     = fltarr(Ntuning, Nscans) ; Single exposure time
  nsum_array     = lonarr(Ntuning, Nscans) ; Number of summed exposures
  
  wcs = replicate({  wave:dblarr(2,2) $
                     , hplt:dblarr(2,2) $
                     , hpln:dblarr(2,2) $
                     , time:dblarr(2,2) $
                  }, Ntuning, Nscans)

  ;; The narrowband cube is aligned to the wideband cube and all
  ;; narrowband scan positions are aligned to each other. So get hpln
  ;; and hplt from the wideband cube wcs coordinates, this should
  ;; apply to all frames in a scan.
  for iscan = 0L, Nscans-1 do begin
    for ituning = 0, Ntuning-1 do begin
      ;; We rely here on hpln and hplt being the first two tabulated
      ;; coordinates. To make this more general, we should get the
      ;; actual indices from the headers. Maybe later...
      wcs[ituning, iscan].hpln = reform(wwcs[0,*,*,iscan])
      wcs[ituning, iscan].hplt = reform(wwcs[1,*,*,iscan])
    endfor                      ; ituning
  endfor                        ; iscan
  
  iprogress = 0
  Nprogress = Nscans*Ntuning
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

      ;; Load prefilters
      
      ;; T camera

      pfile = self.out_dir + '/prefilter_fits/' + instrument.capwords() + '-T_' $
              + unbprefs[ipref] + fitpref_t + 'prefilter.idlsave'
      if ~file_test(pfile) then begin
        red_message, 'Prefilter file not found: '+pfile
        return
      endif
      restore, pfile            ; Restores variable prf which is a struct

                                ; [m] Shift the wavelengths by this amount
      if n_elements(prf.fitpars) gt 1 then wave_shift = prf.fitpars[1]/10. else wave_shift = 0.
;      nbt_units = prf.units
      if n_elements(prf.wav) eq 1 then begin
        nbt_prefilter_curve = prf.pref
      endif else begin
        nbt_prefilter_curve = red_intepf(prf.wav, prf.pref, nbtstates[utunindx[sortindx]].tun_wavelength*1.d10)
      endelse
                                ;  nbt_prefilter_curve = prf.pref
;  nbt_prefilter_wav = prf.wav
;  nbt_prefilter_wb = prf.wbint
      
      nbt_rpref = 1.d0/nbt_prefilter_curve

      ;; R camera

      pfile = self.out_dir + '/prefilter_fits/'+instrument.capwords() + '-R_' $
              + unbprefs[ipref] + fitpref_t + 'prefilter.idlsave'
      if ~file_test(pfile) then begin
        red_message, 'Prefilter file not found: '+pfile
        return
      endif
      restore, pfile            ; Restores variable prf which is a struct

;      nbr_units = prf.units
      if n_elements(prf.wav) eq 1 then begin
        nbr_prefilter_curve = prf.pref
      endif else begin
        nbr_prefilter_curve = red_intepf(prf.wav, prf.pref, nbrstates[utunindx[sortindx]].tun_wavelength*1.d10)
      endelse
                                ;  nbr_prefilter_curve = prf.pref
;  nbr_prefilter_wav = prf.wav
;  nbr_prefilter_wb = prf.wbint
      
      nbr_rpref = 1.d0/nbr_prefilter_curve

      prefilter_curve = (nbt_prefilter_curve + nbr_prefilter_curve)/2.
;  stop
      



      
      ;; Get some metadata from the Stokes cube header.
      red_fitspar_getdates, stokhdr $
                            , date_beg = date_beg $
                            , date_avg = date_avg $
                            , date_end = date_end 
      
      tbeg_array[ituning, iscan] = red_time2double((strsplit(date_beg,'T',/extract))[1])
      tavg_array[ituning, iscan] = red_time2double((strsplit(date_avg,'T',/extract))[1])
      tend_array[ituning, iscan] = red_time2double((strsplit(date_end,'T',/extract))[1])
      
      date_beg_array[ituning, iscan] = date_beg
      date_avg_array[ituning, iscan] = date_avg
      date_end_array[ituning, iscan] = date_end
      
      ;; Wavelength and time
      wcs[ituning, iscan].wave = red_fitsgetkeyword(stokhdr, 'CRVAL3')
      wcs[ituning, iscan].time = red_fitsgetkeyword(stokhdr, 'CRVAL5')
      ;; The preceding lines should perhaps be replaced with
      ;; properly reading the WCS info?
      
      ;; Exposure time
      exp_array[ituning, iscan]  = red_fitsgetkeyword(stokhdr, 'XPOSURE')
      sexp_array[ituning, iscan] = red_fitsgetkeyword(stokhdr, 'TEXPOSUR')
      nsum_array[ituning, iscan] = red_fitsgetkeyword(stokhdr, 'NSUMEXP')
      
      
;      red_show,wb,w=0
;      red_show,wbim,w=1
;      blink, [0, 1]
      
      ;; Apply derot, align, dewarp based on the output from
      ;; make_wb_cube

      for istokes = 0, Nstokes-1 do begin
        if n_elements(fov_mask) gt 0 then nbim[*, *, istokes] = nbim[*, *, istokes] * fov_mask                    
        
        red_missing, nbim[*, *, istokes] $
                     , nmissing = Nmissing, indx_missing = indx_missing, indx_data = indx_data
        if Nmissing gt 0 then begin
          bg = (nbim[*, *, istokes])[indx_missing[0]]
        endif else begin
          bg = median(nbim[*, *, istokes])
        endelse
        
        frame = red_rotation(nbim[*, *, istokes], ang[iscan], $
                             wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF , $
                             background = bg, nearest = nearest, $
                             stretch_grid = reform(wcGRID[iscan,*,*,*])*sclstr, nthreads=nthreads)
        
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
                               , iscan = iscan, ituning = ituning
      endif
      
      iprogress++               ; update progress counter
      
    endfor                      ; ituning

    if ~keyword_set(nocavitymap) then begin
      if ~keyword_set(nomissing_nans) then bg=!Values.F_NaN
      ;; Apply the same derot, align, dewarp as for the science data
      cmap1 = cavitymaps_pref[*, *, ipref]
      cmap11 = red_rotation(cmap1, ang[iscan] $
                            , wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF $
                            , stretch_grid = reform(wcGRID[iscan,*,*,*])*sclstr $
                            , nthreads=nthreads, background=bg)
      
      cavitymaps[0, 0, 0, 0, iscan] = cmap11
      
    endif 
    
  endfor                        ; iscan
  
  ;; Apply wavelength shift from prefilter fit.
  wcs.wave -= wave_shift
  
  ;; Close fits file.
  free_lun, lun
  whdr = headfits(wcfile)
  csyer_spatial_value = fxpar(whdr, 'CSYER1', comment = csyer_spatial_comment)
  red_fitscube_addwcs, filename, wcs $
                       , csyer_spatial_value = csyer_spatial_value $
                       , csyer_spatial_comment = csyer_spatial_comment $
                       , dimensions = dims
  if keyword_set(wbsave) then begin 
    free_lun, wblun
    red_fitscube_addwcs, wbfilename, wcs $
                         , csyer_spatial_value = csyer_spatial_value $
                         , csyer_spatial_comment = csyer_spatial_comment $
                         , dimensions = dims
  endif 

  ;; Copy the MOMFBD or bypass_momfbd step (or a mix):
  self -> headerinfo_copystep, filename, wcfile, stepnum = 1
  

  ;; Need to copy prpara info from the stokes cube headers to hdr.
  ;;hh = red_readhead(snames[0, 0])
  ;;self -> headerinfo_copystep, filename, hh, prstep = 'DEMODULATION'
  self -> headerinfo_copystep, filename, snames[0, 0] , prstep = 'DEMODULATION'
  ;; Should we examine all the stokes cubes and make sure they are
  ;; done using the same parameters? Or else either stop or note
  ;; somehow in the nb_cube header that they differ. (There is a
  ;; mechanism for that now! Except that potentially this is not per
  ;; scan but per Stokes cube.)


  ;; Add info about this step
  hdr = headfits(filename)
  self -> headerinfo_addstep, hdr $
     , prstep = 'CONCATENATION' $
     , prpara = prpara $
     , prproc = inam $
     , prref = 'Align reference: '+wcfile $
     , comment_prref = 'WB cube file name'

  self -> headerinfo_addstep, hdr $
     , prstep = 'CALIBRATION-INTENSITY-SPECTRAL' $
     , prpara = prpara $
     , prref = ['Hamburg FTS spectral atlas (Neckel 1999)' $
                , 'Calibration data from '+red_timestring(prf.time_avg, n = 0)] $
     , prproc = inam
  red_fitscube_newheader, filename, hdr
  if keyword_set(wbsave) then red_fitscube_newheader, wbfilename, hdr
  
  
  ;; Add cavity maps as WAVE distortions 
  if ~keyword_set(nocavitymap) then red_fitscube_addcmap, filename, cavitymaps

  red_message, 'Add some variable keywords.'

  ;; Add some variable keywords
  self -> fitscube_addvarkeyword, filename, 'DATE-BEG', date_beg_array $
     , anchor = anchor $
     , comment = 'Beginning time of observation' $
     , keyword_method = 'first' $
     , axis_numbers = [3, 5] 
  self -> fitscube_addvarkeyword, filename, 'DATE-END', date_end_array $
     , anchor = anchor $
     , comment = 'End time of observation' $
     , keyword_method = 'last' $
     , axis_numbers = [3, 5] 
  self -> fitscube_addvarkeyword, filename, 'DATE-AVG', date_avg_array $
     , anchor = anchor $
     , comment = 'Average time of observation' $
     , keyword_value = self.isodate + 'T' + red_timestring(mean(tavg_array)) $
     , axis_numbers = [3, 5] 

  ;; Needs prefilter_curve for entire scan
  red_fitscube_addrespappl, filename, prefilter_curve, /tun
  
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
     , axis_numbers = [3, 5] 

  self -> fitscube_addvarkeyword, filename, 'TEXPOSUR', sexp_array $
     , comment = '[s] Single-exposure time' $
     , anchor = anchor $
     , tunit = 's' $
     , keyword_method = 'median' $
     , axis_numbers = [3, 5] 

  self -> fitscube_addvarkeyword, filename, 'NSUMEXP', nsum_array $
     , comment = 'Number of summed exposures' $
     , anchor = anchor $
     , keyword_method = 'median' $
     , axis_numbers = [3, 5]
  
  ;; Correct intensity with respect to solar elevation and exposure
  ;; time.
  self -> fitscube_intensitycorr, filename, intensitycorrmethod = intensitycorrmethod $
     , fitpref_time = fitpref_time 
  
  if ~keyword_set(nocrosstalk) then begin

    print, 'Press any key to make crosstalk correction'
    q=get_kbrd()
    ;; Correct the cube for cross-talk, I --> Q,U,V.
    self -> fitscube_crosstalk, filename
    
  endif                         ;else

  if ~keyword_set(nomissing_nans) then begin
    ;; Set padding pixels to missing-data, i.e., NaN.
    self -> fitscube_missing, filename $
       , /noflip $
       , missing_type = 'nan'
    if keyword_set(wbsave) then begin
      self -> fitscube_missing, wbfilename $
         , /noflip $
         , missing_type = 'nan'
    endif
  endif

  
  if keyword_set(integer) then begin
    self -> fitscube_integer, filename $
       , /delete $
;                              , flip = ~keyword_set(noflipping) $
       , outname = outname $
       , overwrite = overwrite
    filename = outname
  endif
  
  
  if ~keyword_set(noflipping) then $
     red_fitscube_flip, filename, flipfile = flipfile $
                        , overwrite = overwrite
  
  red_message, 'Narrowband cube stored in: ' + filename
  if ~keyword_set(noflipping) then print, flipfile
  
  if keyword_set(wbsave) then begin
    if ~keyword_set(noflipping) then red_fitscube_flip, wbfilename, flipfile = wbflipfile $
       , overwrite = overwrite
    red_message, 'Wideband align cube stored in: ' + wbfilename
    if ~keyword_set(noflipping) then print, wbflipfile
  endif
  
end


a = crispred(/dev)

a -> red::make_nb_cube, 'cubes_wb2/wb_6302_2016-09-19T09:28:36_09:28:36=0,1_corrected_im.fits', /overwrite, odir = 'test/'

stop
a -> red::make_nb_cube, 'cubes_wb2/wb_6302_2016-09-19T09:28:36_09:28:36=0,1_09:30:20=0-4_corrected_im.fits', /overwrite, odir = 'test/'

end

