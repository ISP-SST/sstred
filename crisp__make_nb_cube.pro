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
;     nocrosstalk : in, optional, type=boolean
;
;       Do not correct the (polarimetric) data cube Stokes components
;       for crosstalk from I to Q, U, V.
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
;     verbose : in, optional, type=boolean
;
;       Some extra screen output.
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
;-
pro crisp::make_nb_cube, wcfile $
                         , clips = clips $
                         , cmap_fwhm = cmap_fwhm $
                         , integer = integer $
                         , noaligncont = noaligncont $
                         , nocavitymap = nocavitymap $
                         , nocrosstalk = nocrosstalk $
                         , noflipping = noflipping $
                         , nopolarimetry = nopolarimetry $
                         , nostatistics = nostatistics $
                         , notimecor = notimecor $
;                         , nthreads = nthreads $
                         , overwrite = overwrite $
                         , redemodulate = redemodulate $
                         , scannos = scannos $
                         , smooth = smooth $
                         , tiles = tiles $
                         , verbose = verbose $
                         , wbsave = wbsave

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  ;; Make prpara
  red_make_prpara, prpara, clips         
  red_make_prpara, prpara, integer
  red_make_prpara, prpara, cmap_fwhm
  red_make_prpara, prpara, noaligncont 
  red_make_prpara, prpara, nocavitymap 
  red_make_prpara, prpara, notimecor 
  red_make_prpara, prpara, nopolarimetry 
  red_make_prpara, prpara, np           
  red_make_prpara, prpara, smooth
  red_make_prpara, prpara, redemodulate
  red_make_prpara, prpara, tiles        
  red_make_prpara, prpara, wcfile

  if n_elements(nthreads) eq 0 then nthreads = 1 ; Default single thread
  
;  ;; How to smooth the modulation matrices.
;  if n_elements(smooth) eq 0 then begin
;    ;; Smooth with Gaussian kernel by default
;    smooth_by_kernel = 5        ; Default width
;    smooth_by_subfield = 0
;  endif else begin
;    ;; The smooth keyword can either be a number, in which case that
;    ;; is the kernel width, or the string "momfbd", in which case we
;    ;; do smoothing by subfield using the MOMFBD-estimated PSFs.
;    if size(smooth, /tname) eq 'STRING' then begin
;      if strlowcase(smooth) eq 'momfbd' then begin
;        ;; If the string "momfbd" (or "MOMFBD"), we will smooth by
;        ;; subfield. 
;        smooth_by_subfield = 1
;      endif else begin
;        ;; Any string except "momfbd" will result in no smoothing. 
;        smooth_by_subfield = 0
;        smooth_by_kernel = 0
;      endelse
;    endif else begin
;      ;; Not a string, then hopefully a number
;      smooth_by_subfield = 0
;      smooth_by_kernel = smooth
;    endelse
;  endelse
  
  
  ;; Default keywords
  if n_elements(cmap_fwhm) eq 0 then fwhm = 7.0
  if n_elements(tiles) eq 0 or n_elements(clips) eq 0 then begin
    tiles = [8, 16, 32, 64, 128]
    clips = [8, 4,  2,  1,  1  ]
  endif

  if keyword_set(redemodulate) then overwrite = 1

  ;; Camera/detector identification
  self->getdetectors
  wbindx      = where(strmatch(*self.cameras,'Crisp-W'))
  wbcamera    = (*self.cameras)[wbindx[0]]
  wbdetector  = (*self.detectors)[wbindx[0]]
  nbtindx     = where(strmatch(*self.cameras,'Crisp-T')) 
  nbtcamera   = (*self.cameras)[nbtindx[0]]
  nbtdetector = (*self.detectors)[nbtindx[0]]
  nbrindx     = where(strmatch(*self.cameras,'Crisp-R')) 
  nbrcamera   = (*self.cameras)[nbrindx[0]]
  nbrdetector = (*self.detectors)[nbrindx[0]]

  ;; We currently do correct for the small scale cavity map in CRISP
  ;; data. (We should get this from earlier meta data!)
  remove_smallscale = 1
  
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
            , ['ANG', 'CROP', 'FF', 'GRID', 'ND', 'SHIFT', 'TMEAN', 'X01Y01'] $
            ,   ANG, wcCROP, wcFF, wcGRID, wcND, wcSHIFT, wcTMEAN, wcX01Y01
  ;; Note that the strarr wfiles cannot be read by fxbreadm! Put it in
  ;; wbgfiles (WideBand Global).
  fxbread, bunit, wbgfiles, 'WFILES', 1
  fxbclose, bunit

  ;; CRISP has a common prefilter for WB and NB
  prefilter = fxpar(wchead,'FILTER1')

  ;; Read wcs extension of wb file to get pointing info
  fxbopen, wlun, wcfile, 'WCS-TAB', wbdr
  ttype1 = fxpar(wbdr, 'TTYPE1')
  fxbread, wlun, wwcs, ttype1
  fxbclose, wlun

  x0 = wcX01Y01[0]
  x1 = wcX01Y01[1]
  y0 = wcX01Y01[2]
  y1 = wcX01Y01[3]

  self -> extractstates, wbgfiles, wbgstates
  prefilter = wbgstates[0].prefilter

  
  
  wchdr0 = red_readhead(wbgfiles[0])
  datestamp = strtrim(fxpar(wchdr0, 'STARTOBS'), 2)
  timestamp = (strsplit(datestamp, 'T', /extract))[1]
  
  datadir = file_dirname(wbgfiles[0])+'/'
  extension = (strsplit(wbgfiles[0],'.',/extract))[-1]

  files = file_search(datadir + '*.'+extension, count = Nfiles)      

  ;; If we don't want to make a cube with all scans, this could be
  ;; unnecessarily many files. Takes a long time to do extracstates on
  ;; them.
  undefine,selectfiles 
  Nscans = n_elements(wbgstates.scannumber)
  for iscan = 0, Nscans-1 do begin
    indx = where(strmatch(files, '*_'+string(wbgstates[iscan].scannumber, format = '(i05)')+'_*'), Nmatch)
    if Nmatch ne 0 then red_append, selectfiles, files[indx]
  endfor                        ; iscan
  files = selectfiles
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

;  ufpi_states = pertuningstates[utunindx[sortindx]].fpi_state
;  utunwavelength = pertuningstates[utunindx[sortindx]].tun_wavelength
;  wav = utunwavelength
;  my_prefilters = pertuningstates[utunindx[sortindx]].prefilter
;
;  ;; Unique nb prefilters
;  unbprefindx = uniq(pertuningstates[utunindx].prefilter, sort(pertuningstates[utunindx].prefilter))
;  Nnbprefs = n_elements(unbprefindx)
;  unbprefs = pertuningstates[utunindx[unbprefindx]].prefilter
;  unbprefsref = dblarr(Nnbprefs)
;
;  for inbpref = 0L, Nnbprefs-1 do begin
;    ;; This is the reference point of the fine tuning for this prefilter:
;    unbprefsref[inbpref] = double((strsplit(pertuningstates[utunindx[unbprefindx[inbpref]]].tuning $
;                                            , '_', /extract))[0])
;  endfor                        ; inbpref
;  
;  unbprefsref *= 1e-10          ; [m]

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

  ;; Unique tuning states, sorted by wavelength
  utunindx = uniq(nbrstates.fpi_state, sort(nbrstates.fpi_state))
  Ntuning = n_elements(utunindx)
  sortindx = sort(nbrstates[utunindx].tun_wavelength)
  utuning = nbrstates[utunindx[sortindx]].fpi_state
  
  ;; Make a Stokes cube?
  ulc = nbrstates[uniq(nbrstates.lc, sort(nbrstates.lc))].lc
  Nlc = n_elements(ulc)
  makestokes = (Nlc gt 1) and ~keyword_set(nopolarimetry)

  if makestokes then Nstokes = 4 else Nstokes = 1
  
  if Nnbt ne Nnbr then stop

  ;; Prepare for making output file names
  if n_elements(odir) eq 0 then odir = self.out_dir + '/cubes_nb/' 
  ofile = 'nb_' + prefilter + '_' + datestamp + '_scans=' $ 
          + red_collapserange(uscans, ld = '', rd = '')
  if makestokes then ofile += '_stokes'
  ofile += '_corrected'
;  if keyword_set(integer) then ofile += '_int'
  ofile += '_im.fits'
  filename = odir+ofile

  ;; Already done?
  if file_test(filename) then begin
    if keyword_set(overwrite) then begin
      print, inam + ' : Overwriting existing data cube:'
      print, filename
    endif else begin
      print, inam + ' : This data cube exists already:'
      print, filename
      return
    endelse
  endif

  file_mkdir, odir

  
  ;; Do WB correction?
  wbcor = Nwb eq Nnbt 

  ;; Load prefilters
  
  ;; Crisp-T

  pfile = self.out_dir + '/prefilter_fits/Crisp-T_'+prefilter+'_prefilter.idlsave'
  if ~file_test(pfile) then begin
    print, inam + ' : prefilter file not found: '+pfile
    return
  endif
  restore, pfile                ; Restores variable prf which is a struct

  nbt_units = prf.units
  nbt_prefilter_curve = prf.pref
  nbt_prefilter_wav = prf.wav
  nbt_prefilter_wb = prf.wbint
  
  nbt_rpref = 1.d0/nbt_prefilter_curve

  ;; Crisp-R

  pfile = self.out_dir + '/prefilter_fits/Crisp-R_'+prefilter+'_prefilter.idlsave'
  if ~file_test(pfile) then begin
    print, inam + ' : prefilter file not found: '+pfile
    return
  endif
  restore, pfile                ; Restores variable prf which is a struct

  nbr_units = prf.units  
  nbr_prefilter_curve = prf.pref
  nbr_prefilter_wav = prf.wav
  nbr_prefilter_wb = prf.wbint
  
  nbr_rpref = 1.d0/nbr_prefilter_curve


  
  if nbr_units ne nbt_units then begin
    print, inam + ' : Units for Crisp-T and Crisp-R do not match.'
    print, inam + ' : Please rerun the prefilterfit step for these data.'
    retall
  endif
  units = nbr_units


  ;; Scaling parameter with respect to mean wb intensity (mainly
  ;; varying elevation)
  if keyword_set(notimecor) then begin
    tscl = replicate(1., Nscans)
  endif else begin
    tscl = mean(nbt_prefilter_wb) / wcTMEAN
  endelse

  ;; Load WB image and define the image border
;  tmp = red_readdata(wbgfiles[0])

  ;; Spatial dimensions that match the WB cube
  Nx = wcND[0]
  Ny = wcND[1]

  ;; Create cubes for science data and scan-adapted cavity maps.
  cavitymaps = fltarr(Nx, Ny, Nscans)

  if ~keyword_set(nocavitymap) then begin

    ;; Read the original cavity map
    cfile = self.out_dir + 'flats/spectral_flats/' $
            + strjoin([nbtstates[0].detector $
                       , prefilter $
                       , 'fit_results.sav'] $
                      , '_')

    if ~file_test(cfile) then begin
      print, inam + ' : Error, calibration file not found -> '+cfile
      stop
    endif
    restore, cfile                  ; The cavity map is in a struct called "fit". 
    cmapt = reform(fit.pars[1,*,*]) ; Unit is [Angstrom]
    cmapt /= 10.                    ; Make it [nm]
    cmapt = -cmapt                  ; Change sign so lambda_correct = lambda + cmap
    fit = 0B                        ; Don't need the fit struct anymore.
    
    cfile = self.out_dir + 'flats/spectral_flats/' $
            + strjoin([nbrstates[0].detector $
                       , prefilter $
                       , 'fit_results.sav'] $
                      , '_')

    if ~file_test(cfile) then begin
      print, inam + ' : Error, calibration file not found -> '+cfile
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
    
    ;; Read the output of the pinhole calibrations so we can do the same
    ;; to the cavity maps as was done to the raw data in the momfbd
    ;; step. This output is in a struct "alignments" in the save file
    ;; 'calib/alignments.sav'
    restore, 'calib/alignments.sav'
    ;; Should be based on state1 or state2 in the struct? make_cmaps
    ;; says "just pick one close to continuum (last state?)".
;    indx = where(nbstates[0].prefilter eq alignments.state2.prefilter, Nalign)
    indxt = where(alignments.state2.camera eq 'Crisp-T', Nalignt)
    case Nalignt of
      0    : stop               ; Should not happen!
      1    : amapt = invert(      alignments[indxt].map           )
      else : amapt = invert( mean(alignments[indxt].map, dim = 3) )
    endcase
    indxr = where(alignments.state2.camera eq 'Crisp-R', Nalignr)
    case Nalignr of
      0    : stop               ; Should not happen!
      1    : amapr = invert(      alignments[indxr].map           )
      else : amapr = invert( mean(alignments[indxr].map, dim = 3) )
    endcase

    cmap1r = rdx_img_project(amapr, cmap1r) ; Apply the geometrical mapping
    cmap1r = cmap1r[x0:x1,y0:y1]            ; Clip to the selected FOV
    cmap1t = rdx_img_project(amapt, cmap1t) ; Apply the geometrical mapping
    cmap1t = cmap1t[x0:x1,y0:y1]            ; Clip to the selected FOV

    ;; At this point, the individual cavity maps should be corrected
    ;; for camera misalignments, so they should be aligned with
    ;; respect to the cavity errors on the etalons. So we can sum
    ;; them.
    cmap1 = (cmap1r + cmap1t) / 2.
    
  endif
  
  ;; Make FITS header for the NB cube
  hdr = wchead                      ; Start with the WB cube header
  red_fitsdelkeyword, hdr, 'STATE' ; Not a single state for cube 
  
  
  if makestokes then begin

    ;; Make Stokes cubes
    ;; Define the Stokes file names needed for this nb cube
    snames = strarr(Nscans, Ntuning)    
    for iscan = 0, Nscans-1 do begin
      
      self -> make_stokes_cubes, datadir, uscans[iscan] $
                                 , clips = clips $
                                 , cmap_fwhm = cmap_fwhm $
                                 , nocavitymap = nocavitymap $
                                 , redemodulate = redemodulate $
                                 , smooth = smooth $
                                 , snames = these_snames $
                                 , stokesdir = stokesdir $
                                 , tiles = tiles

      snames[iscan, *] = these_snames
      
    endfor                      ; iscan

    
;    ;;  stokesdir = datadir + '/stokes/'
;
;    ;; Store intermediate Stokes cubes in separate directories for
;    ;; different smooth options: 
;    stokesdir = datadir + '/stokes_sbs'+strtrim(smooth_by_subfield,2) $
;                + '_sbk'+strtrim(smooth_by_kernel,2)+'/'
;
;    file_mkdir, stokesdir
;
;    if keyword_set(redemodulate) then begin
;      ;; We will delete all existing stokesIQUV*.fits files. An
;      ;; alternative would be to delete only the ones involved in the
;      ;; cube to be made.
;
;      dfiles = file_search(stokesdir+'stokesIQUV*.fits', count = Ndelete)
;      if Ndelete gt 0 then begin
;        print, inam + ' : Will delete the following Stokes files:'
;        print, dfiles, format = '(a0)'
;        file_delete, dfiles
;      endif
;    endif
;
;    ;; Define the Stokes file names needed for this nb cube
;    snames = strarr(Nscans, Ntuning)    
;    for iscan = 0, Nscans-1 do begin
;      for ituning = 0, Ntuning-1 do begin
;        snames[iscan, ituning] = stokesdir $
;                                 + strjoin(['stokesIQUV' $
;                                            , string(uscans[iscan], format = '(i05)') $
;                                            , prefilter $
;                                            , utuning[ituning] $
;                                           ], '_') + '.fits' 
;      endfor                    ; ituning
;    endfor                      ; iscan
;
;    ;; Make Stokes cubes for each scan and state (if not done already) 
;    todoindx = where(~file_test(snames), Ntodo)
;    if Ntodo gt 0 then begin
;      print, inam + ' : Will have to make '+strtrim(Ntodo, 2) + ' Stokes cubes.'
;
;;      ;; Get the FOV in the momfbd files.
;;      mr = momfbd_read(wbgfiles[0], /names) ; Use /names to avoid reading the data parts
;;      mrX01Y01 = mr.roi + mr.margin * [1, -1, 1, -1]
;      
;      ;; First get the inverse modulation matrices, make them if
;      ;; needed. They are returned in the (size and) orientation of
;      ;; the momfbd output.
;      self -> inverse_modmatrices, prefilter, stokesdir $
;                                   , camr = nbrcamera, immr = immr $
;                                   , camt = nbtcamera, immt = immt $
;                                   , no_ccdtabs = no_ccdtabs
;
;      swcs = {wave:dblarr(2,2)   $ ; WCS for this Stokes cube.
;              , hplt:dblarr(2,2) $
;              , hpln:dblarr(2,2) $
;              , time:dblarr(2,2) $
;             }
;
;      if Nthreads gt 1 then begin
;
;        ;; Make Nthreads bridges
;        bridges = build_bridges(Nthreads)
;        stop
;        
;        ;; Initialize the bridges
;        for ithread = 0, Nthreads-1 do begin
;
;          red_progressbar, ithread, Nthreads, 'Set up '+strtrim(Nthreads, 2)+' IDL bridges'
;          
;          bridges[ithread] -> execute, 'a=crispred(dev='+strtrim(long(self.developer_mode), 2)+')' ; Set dev!
;
;          ;; Set constant variables
;          bridges[ithread] -> setvar, 'immr', immr
;          bridges[ithread] -> setvar, 'immt', immt  
;          bridges[ithread] -> setvar, 'smooth_by_subfield', smooth_by_subfield
;          bridges[ithread] -> setvar, 'smooth_by_kernel', smooth_by_kernel      
;          bridges[ithread] -> setvar, 'clips', clips
;          bridges[ithread] -> setvar, 'cmap', cmap1
;          bridges[ithread] -> setvar, 'overwrite', keyword_set(redemodulate)
;          bridges[ithread] -> setvar, 'tiles', tiles
;          bridges[ithread] -> setvar, 'units', units
;
;        endfor                  ; ithread
;      endif
;      
;      
;      itodo = 0
;      for iscan = 0, Nscans-1 do begin
;
;        undefine, wbg 
;        
;        for ituning = 0, Ntuning-1 do begin
;
;          if ~file_test(snames[iscan, ituning]) then begin
;
;            ;; Read the global WB file for this scan.
;            if n_elements(wbg) eq 0 then wbg = red_readdata(wbgfiles[iscan])
;            
;            red_progressbar, itodo, Ntodo, /predict  $
;                             , 'Making '+snames[iscan, ituning] 
;
;            self -> selectfiles, files = wbfiles, states = wbstates $
;                                 , sel = these_wbindx, count = Nthesewb $
;                                 , scan = uscans[iscan] $
;                                 , fpi_states = utuning[ituning]
;
;            self -> selectfiles, files = nbtfiles, states = nbtstates $
;                                 , sel = these_nbtindx, count = Nthesenbt $
;                                 , scan = uscans[iscan] $
;                                 , fpi_states = utuning[ituning]
;
;            self -> selectfiles, files = nbrfiles, states = nbrstates $
;                                 , sel = these_nbrindx, count = Nthesenbr $
;                                 , scan = uscans[iscan] $
;                                 , fpi_states = utuning[ituning]
;
;            swcs.hpln = reform(wwcs[0,*,*,iscan])
;            swcs.hplt = reform(wwcs[0,*,*,iscan])
;            swcs.wave = nbtstates[these_nbrindx[0]].tun_wavelength*1d9
;            ;; swcs.time = ; Set by demodulate
;            
;            ;; The demodulate method reads the momfbd output wb, nbt,
;            ;; and wbr images for a particular scan and tuning state
;            ;; and outputs a demodulated Stokes file.
;            if Nthreads gt 1 then begin
;              stop
;              bridge = get_idle_bridge(bridges)
;
;              ;; Needed only for callback:
;;              ud = {i:i,j:j,pout:pout} 
;;              bridge -> setproperty, userdata=ud
;
;
;              ;; Set varying variables
;;              bridge -> setvar, 'in', in[i,j]
;              bridge -> setvar, 'nbrfac', nbr_rpref[ituning]
;              bridge -> setvar, 'nbrstates', nbrstates[these_nbrindx] ; IDL_IDLBRIDGE Error: Unsupported parameter type: IDL_TYP_STRUCT
;
;              bridge -> setvar, 'nbtfac', nbt_rpref[ituning]
;              bridge -> setvar, 'nbtstates', nbtstates[these_nbtindx]
;              bridge -> setvar, 'snames', snames[iscan, ituning]
;              bridge -> setvar, 'wbg', wbg
;              bridge -> setvar, 'wbstates', wbstates[these_wbindx]
;              bridge -> setvar, 'wcs', swcs
;
;
;              ;; Send the demodulation task to the bridge
;;              bridge -> execute, /nowait, 'worker, in, out'
;              bridge -> execute, /nowait $
;                                 , 'self -> demodulate, snames, immr, immt' $
;                                 + ', clips = clips' $                                  
;                                 + ', cmap = cmap' $                                   
;                                 + ', nbrfac = nbrfac' $                       
;                                 + ', nbrstates = nbrstates' $           
;                                 + ', nbtfac = nbtfac' $                       
;                                 + ', nbtstates = nbtstates' $           
;                                 + ', overwrite = overwrite' $                       
;                                 + ', smooth_by_kernel = smooth_by_kernel' $            
;                                 + ', smooth_by_subfield = smooth_by_subfield' $        
;                                 + ', tiles = tiles' $                                  
;                                 + ', units = units' $                                  
;                                 + ', wbg = wbg' $                                      
;                                 + ', wbstates = wbstates' $                
;                                 + ', wcs = wcs'                            
;              
;            endif else begin
;              self -> demodulate, snames[iscan, ituning], immr, immt $
;                                  , smooth_by_subfield = smooth_by_subfield $ 
;                                  , smooth_by_kernel = smooth_by_kernel $ 
;                                  , clips = clips $
;                                  , cmap = cmap1 $
;                                  , nbrfac = nbr_rpref[ituning] $
;                                  , nbrstates = nbrstates[these_nbrindx] $
;                                  , nbtfac = nbt_rpref[ituning] $
;                                  , nbtstates = nbtstates[these_nbtindx] $
;                                  , overwrite = redemodulate $
;                                  , tiles = tiles $
;                                  , units = units $
;                                  , wbg = wbg $
;                                  , wcs = swcs $
;                                  , wbstates = wbstates[these_wbindx]
;              
;            endelse
;            
;            itodo++
;            
;          endif
;          
;        endfor                  ; ituning 
;      endfor                    ; iscan 
;
;      ;; Wait for bridges to finish, then destroy them
;      if Nthreads gt 1 then begin
;        barrier_bridges, bridges
;        burn_bridges, bridges
;      endif
;      
;    endif

    ;; Need to copy some info from the stokes cube headers to hdr.
    ;; E.g. the prpara info.

    
    hh = red_readhead(snames[0, 0])
    self -> headerinfo_copystep, hdr, hh, /last

    ;; At some point we also need to examine all the stokes cubes and
    ;; make sure they are done using the same parameters. Or else
    ;; either stop or note somehow in the nb_cube header that they
    ;; differ.
    
  endif

  
  
  red_fitsaddkeyword, hdr, 'BITPIX', -32
    
  ;; Add info about this step
  self -> headerinfo_addstep, hdr $
                              , prstep = 'Prepare NB science data cube' $
                              , prpara = prpara $
                              , prproc = inam

  anchor = 'DATE'

  ;; Add some keywords
  red_fitsaddkeyword, anchor = anchor, hdr, 'OBS_HDU', 1
  
  ;; Add info to headers
  red_fitsaddkeyword, anchor = anchor, hdr, 'BUNIT', units, 'Units in array'
  red_fitsaddkeyword, anchor = anchor, hdr, 'BTYPE', 'Intensity', 'Type of data in array'

  ;; WB and NB data come from different cameras.
  red_fitsaddkeyword, hdr, 'CAMERA', nbtcamera + ',' + nbrcamera
  ;; Get DETGAIN, DETOFFS, DETMODEL, DETFIRM from .fitsheader file,
  ;; i.e., red_readhead(.momfbd file). This would have to be handled
  ;; differently with CRISP because each Stokes image is a mix of two
  ;; detectors. 
  mhdr = red_readhead(nbrfiles[0]) ; Header of momfbd output file
  mhdt = red_readhead(nbtfiles[0]) ; Header of momfbd output file
  red_fitsaddkeyword, hdr, 'DETECTOR', strtrim(fxpar(mhdr,'DETECTOR'), 2) $
                      + ',' + strtrim(fxpar(mhdt,'DETECTOR'), 2)
;  red_fitsaddkeyword, hdr, 'DETGAIN',  fxpar(mhdr,'DETGAIN') + ',' + fxpar(mhdt,'DETGAIN')
;  red_fitsaddkeyword, hdr, 'DETOFFS',  fxpar(mhdr,'DETOFFS') + ',' + fxpar(mhdt,'DETOFFS')
;  red_fitsaddkeyword, hdr, 'DETMODEL', fxpar(mhdr,'DETMODEL') + ',' + fxpar(mhdt,'DETMODEL')
;  red_fitsaddkeyword, hdr, 'DETFIRM',  fxpar(mhdr,'DETFIRM') + ',' + fxpar(mhdt,'DETFIRM')

  ;; Initialize fits file, set up for writing the data part.
  dims = [Nx, Ny, Ntuning, Nstokes, Nscans] 
  self -> fitscube_initialize, filename, hdr, lun, fileassoc, dims 

  if keyword_set(wbsave) then begin
    wbfilename = strreplace(filename, 'nb_', 'wbalign_')
    wbdims = [Nx, Ny, Ntuning, 1, Nscans] 
    wbhdr = hdr
    self -> fitscube_initialize, wbfilename, wbhdr, wblun, wbfileassoc, wbdims 
  endif
  

  ;; Observations metadata varaibles
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
    wb = (red_readdata(wbgfiles[iscan]))[x0:x1, y0:y1]

    for ituning = 0L, Ntuning - 1 do begin 

      red_progressbar, iprogress, Nprogress $
                       , /predict $
                       , 'Processing scan=' $
                       + strtrim(uscans[iscan], 2) + ' tuning=' + utuning[ituning] 


      if makestokes then begin

        ;; Read the pre-made Stokes cube
        
        nbim = red_readdata(snames[iscan, ituning], head = stokhdr)
        nbim = reform(nbim[x0:x1, y0:y1, 0, *]) * tscl[iscan]

        ;; The Stokes cube is already wbcorrected so we should not
        ;; repeat that here.

        ;; Get some metadata from the Stoeks cube header.
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

      endif else begin
        ;; The NB files in this scan, sorted in tuning wavelength order.
        self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                             , fpi_states = utuning[ituning] $
                             , cam = nbtcamera, scan = uscans[iscan] $
                             , sel = scan_nbtindx, count = count
        scan_nbtfiles = pertuningfiles[scan_nbtindx]
        scan_nbtstates = pertuningstates[scan_nbtindx]
        sortindx = sort(scan_nbtstates.tun_wavelength)
        scan_nbtfiles = scan_nbtfiles[sortindx]
        scan_nbtstates = scan_nbtstates[sortindx]
        self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                             , fpi_states = utuning[ituning] $
                             , cam = nbrcamera, scan = uscans[iscan] $
                             , sel = scan_nbrindx, count = count
        scan_nbrfiles = pertuningfiles[scan_nbrindx]
        scan_nbrstates = pertuningstates[scan_nbrindx]
        sortindx = sort(scan_nbrstates.tun_wavelength)
        scan_nbrfiles = scan_nbrfiles[sortindx]
        scan_nbrstates = scan_nbrstates[sortindx]
;    self -> selectfiles, files = pertuningfiles, states = pertuningstates $
;                         , cam = nbcamera, scan = uscans[iscan] $
;                         , sel = scan_nbindx, count = count
;    scan_nbfiles = pertuningfiles[scan_nbindx]
;    scan_nbstates = pertuningstates[scan_nbindx]

        ;; The WB files in this scan, sorted as the NB files
        self -> selectfiles, files = pertuningfiles, states = pertuningstates $
                             , fpi_states = utuning[ituning] $
                             , cam = wbcamera, scan = uscans[iscan] $
                             , sel = scan_wbindx, count = count
        scan_wbfiles = pertuningfiles[scan_wbindx]
        scan_wbstates = pertuningstates[scan_wbindx]
        ;;match2, scan_nbrstates.fpi_state, scan_wbstates.fpi_state, sortindx
        match2, scan_nbrstates.lc, scan_wbstates.lc, sortindx
        scan_wbfiles = scan_wbfiles[sortindx]
        scan_wbstates = scan_wbstates[sortindx]

        ;; Collect info about this frame here.
        for iim = 0, n_elements(scan_wbindx)-1 do begin
          wbhead = red_readhead(scan_wbfiles[iim])
          red_fitspar_getdates, wbhead $
                                , date_beg = date_beg $
                                , date_avg = date_avg $
                                , date_end = date_end 
          red_append, tbegs, red_time2double((strsplit(date_beg,'T',/extract))[1])
          red_append, tavgs, red_time2double((strsplit(date_avg,'T',/extract))[1])
          red_append, tends, red_time2double((strsplit(date_end,'T',/extract))[1])
          red_append, xps, fxpar(wbhead, 'XPOSURE')
          red_append, texps, fxpar(wbhead, 'TEXPOSUR')
          red_append, nsums, fxpar(wbhead, 'NSUMEXP')
        endfor
        tbeg_array[ituning, iscan] = min(tbegs)
        tavg_array[ituning, iscan] = mean(tavgs)
        tend_array[ituning, iscan] = max(tends)

        date_beg_array[ituning, iscan] = self.isodate + 'T' + red_timestring(tbeg_array[ituning, iscan])
        date_avg_array[ituning, iscan] = self.isodate + 'T' + red_timestring(tavg_array[ituning, iscan])
        date_end_array[ituning, iscan] = self.isodate + 'T' + red_timestring(tend_array[ituning, iscan])

        ;; Wavelength and time
        wcs[ituning, iscan].wave = scan_nbtstates[0].tun_wavelength*1d9
        wcs[ituning, iscan].time = tavg_array[ituning, iscan]

        ;; Exposure time
        exp_array[ituning, iscan]  = total(xps)
        sexp_array[ituning, iscan] = mean(texps)
        nsum_array[ituning, iscan] = round(total(nsums))


;      ;; The NB files in this scan, sorted in tuning wavelength order.
;      self -> selectfiles, files = pertuningfiles, states = pertuningstates $
;                           , fpi_states = utuning[ituning] $
;                           , cam = nbtcamera, scan = uscans[iscan] $
;                           , sel = scan_nbtindx, count = count
;      scan_nbtfiles = pertuningfiles[scan_nbtindx]
;      scan_nbtstates = pertuningstates[scan_nbtindx]
;      self -> selectfiles, files = pertuningfiles, states = pertuningstates $
;                           , fpi_states = utuning[ituning] $
;                           , cam = nbrcamera, scan = uscans[iscan] $
;                           , sel = scan_nbrindx, count = count
;      scan_nbrfiles = pertuningfiles[scan_nbrindx]
;      scan_nbrstates = pertuningstates[scan_nbrindx]
        
        
        ;; Read images, apply prefilter curve and temporal scaling.
        ;; Average all LC states for the two cameras.
        nbim = 0.
        Nim = n_elements(scan_nbrindx) + n_elements(scan_nbtindx)
        for iim = 0, n_elements(scan_nbtindx)-1 do begin

          this_im = (red_readdata(scan_nbrfiles[iim]))[x0:x1, y0:y1] * nbr_rpref[ituning]
          if wbcor then begin
            ;; Get destretch to anchor camera (residual seeing)
            wwi = (red_readdata(scan_wbfiles[iim]))[x0:x1, y0:y1]
            grid1 = red_dsgridnest(wb, wwi, tiles, clips)
            ;; Apply destretch to anchor camera and prefilter correction
            this_im = red_stretch(temporary(this_im), grid1)
          endif
          nbim += this_im

          this_im = (red_readdata(scan_nbtfiles[iim]))[x0:x1, y0:y1] * nbt_rpref[ituning]
          if wbcor then begin
            ;; Apply destretch to anchor camera and prefilter correction
            this_im = red_stretch(temporary(this_im), grid1)
          endif
          nbim += this_im

        endfor
        nbim *= tscl[iscan] / Nim

      endelse
      

      if makestokes then begin
        nbcube = fltarr(Nx, Ny, Nstokes)
        for istokes = 0, Nstokes-1 do begin
          ;; Apply derot, align, dewarp based on the output from
          ;; make_wb_cube
          nbcube[0, 0, istokes] = red_rotation(nbim[*, *, istokes], ang[iscan], $
                                               wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF)
          nbcube[0, 0, istokes] = red_stretch(nbcube[*, *, istokes] $
                                              , reform(wcGRID[iscan,*,*,*]))
          self -> fitscube_addframe, fileassoc, nbcube[*, *, istokes] $
                                     , iscan = iscan, ituning = ituning, istokes = istokes
        endfor                  ; istokes
      endif else begin
        nbim = red_rotation(nbim, ang[iscan] $
                            , wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF)
        nbim = red_stretch(nbim, reform(wcGRID[iscan,*,*,*]))
        self -> fitscube_addframe, fileassoc, nbim $
                                   , iscan = iscan, ituning = ituning
      endelse
      
      if keyword_set(wbsave) then begin
        ;; Same operations as on narrowband image, except for
        ;; "aligncont".
        wbim = wwi * tscl[iscan]
        wbim = red_stretch(temporary(wbim), grid1)
        wbim = red_rotation(temporary(wbim), ang[iscan], $
                            wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF)
        wbim = red_stretch(temporary(wbim), reform(wcGRID[iscan,*,*,*]))
        self -> fitscube_addframe, wbfileassoc, temporary(wbim) $
                                   , iscan = iscan, ituning = ituning
      endif
      
      iprogress++               ; update progress counter

    endfor                      ; ituning

    if ~keyword_set(nocavitymap) then begin
      
      ;; Apply the same derot, align, dewarp as for the science data
      cmap11 = red_rotation(cmap1, ang[iscan] $
                            , wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF)
      cmap11 = red_stretch(temporary(cmap11), reform(wcGRID[iscan,*,*,*]))
      
      cavitymaps[0, 0, iscan] = cmap11

      ;; The following block of code is inactive but we want to keep
      ;; it around in case it is needed later. It does further
      ;; operations on the cavity maps based on what is done to the
      ;; science data, such as stretching based on the extra
      ;; per-tuning wideband objects, as well as blurring based on the
      ;; momfbd psfs. It should probably have a boolean keyword that
      ;; activates it. (This code will have to be updated to the
      ;; current pipeline style before it can be used.)
      ;;
      ;; Note that in this case, for polarimetric data, we might also
      ;; want to read the cavity maps from the single Stokes cubes,
      ;; which are averages of the original cmap, stretched as the
      ;; individual LC state images before combining to form the
      ;; Stokes components.
      if 0 then begin
        wb = (red_readdata(wbf[ss]))[x0:x1,y0:y1]
        for ww = 0L, nw - 1 do begin
          
          iwbf = strjoin([self.camwbtag, (strsplit(file_basename(st.ofiles[ww,ss]),'.',/extract))[1:*]],'.')
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
          
        endfor                  ; ww
      endif
    endif 
    
    if(keyword_set(verbose)) then begin
      print, inam +'scan=',iscan,', max=', max(d1)            
    endif
    
  endfor                        ; iscan

  
  ;; Close fits file.
  self -> fitscube_finish, lun, wcs = wcs
  if keyword_set(wbsave) then self -> fitscube_finish, wblun, wcs = wcs

  ;; Add cavity maps as WAVE distortions 
  if ~keyword_set(nocavitymap) then self -> fitscube_addcmap, filename, cavitymaps

  print, inam + ' : Add some variable keywords.'

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
  
  if makestokes && ~keyword_set(nocrosstalk) then begin

    ;; Correct the cube for cross-talk, I --> Q,U,V.
    self -> fitscube_crosstalk, filename, nostatistics = nostatistics 

  endif else if ~keyword_set(nostatistics) then begin

    ;; Calculate statistics if not done already
    
;    percentiles = [.01, .10, .25, .50, .75, .90, .95, .98, .99]
    
    red_fitscube_statistics, filename, /write $
                             , angles = ang $
                             , full = wcFF $
                             , grid = wcGRID $
                             , origNx = Nxx $
                             , origNy = Nyy $
;                             , percentiles = percentiles $
                             , shifts = wcSHIFT 
    
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
  
  print, inam + ' : Narrowband cube stored in:'
  print, filename
  if ~keyword_set(noflipping) then print, flipfile
  
  if keyword_set(wbsave) then begin
    if ~keyword_set(noflipping) then red_fitscube_flip, wbfilename, flipfile = wbflipfile $
       , overwrite = overwrite
    print, inam + ' : Wideband align cube stored in:'
    print, wbfilename
    if ~keyword_set(noflipping) then print, wbflipfile
  endif
  
end
