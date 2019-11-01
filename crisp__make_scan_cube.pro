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
;     autocrop : in, optional, type=booean
;
;      Try to determine the largest FOV that avoids any bad momfbd
;      subfields along the edges. If this keyword is set, the input
;      value of the crop keyword is ignored and is set to the
;      auto-detected crop parameters.
; 
;     clips : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment.
;
;     cmap_fwhm : in, type=float, default=7
;   
;       FWHM in pixels of kernel used for smoothing the cavity map.
;
;     crop : in, out, optional, type=array, default="[0,0,0,0]"
;
;       The array given here will be used to limit the FOV to
;       [xl+crop[0],xh-crop[1],yl+[crop[3],yh-crop[3]]. If /autocrop,
;       then the auto detected crop is returned in this keyword
;       instead.
;
;     integer : in, optional, type=boolean
;
;       Store as integers instead of floats.
;
;     interactive : in, optional, type=boolean
;
;       Set this keyword to define the data cube FOV by use of the
;       XROI GUI. If autocrop is set, then use the so defined FOV as
;       an initialization of the FOV in the GUI. Otherwise use the
;       crop keyword (or its default).
;
;     limb_data : in, optional, type=boolean
;
;       Set for data where the limb is in the FOV. Disables autocrop.
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
;     nostatistics : in, optional, type=boolean
;  
;       Do not calculate statistics metadata to put in header keywords
;       DATA*. If statistics keywords already exist, then remove them.
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
;     scannos : in, optional, type=string, default="ask"
;
;       Choose scan numbers to include in the sequence by entering a
;       comma-and-dash delimited string, like '2-5,7-20,22-30' or the
;       string '*' to include all.
;
;     smooth : in, optional, type=varies, default=5
;
;       How to smooth the modulation matrices? Set to the string
;       'momfbd' to smooth subfield by subfield using the PSFs
;       estimated by momfbd. Set to a number to smooth by a Gaussian
;       kernel of that width. 
;
;     tiles : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment. 
;
;     tuning_selection : in, optional, type="integer or array"
;
;       The index or indices in the tuning dimension to use for
;       calculating the correction. Should correspond to continuum (or
;       as close to continuum as possible), where the polarimetric
;       signal is minimal. Negative indices are allowed.
;
; 
; :History:
; 
;    2019-03-21 : MGL. First crisp version, based on the chromis one
;                 and code from crisp::make_nb_cube.
; 
;    2019-08-26 : MGL. Do integerization, statistics calculations, and
;                 WB intensity correction by calling subprograms.
; 
;-
pro crisp::make_scan_cube, dir $
                           , autocrop = autocrop $
                           , clips = clips $
                           , cmap_fwhm = cmap_fwhm $
                           , crop = crop $
                           , integer = integer $
                           , interactive = interactive $
                           , limb_data = limb_data $
                           , nocavitymap = nocavitymap $
                           , nopolarimetry = nopolarimetry $
                           , nowbintensitycorr = nowbintensitycorr $
                           , overwrite = overwrite $
                           , redemodulate = redemodulate $
                           , scannos = scannos $
                           , smooth = smooth $
                           , tiles = tiles  $
                           , tuning_selection = tuning_selection
               
  ;; Name of this method
  inam = strlowcase((reverse((scope_traceback(/structure)).routine))[0])

  if n_elements(nowbintensitycorr) eq 0 nowbintensitycorr = 1 ; Temporary default!
  
  ;; Make prpara
  red_make_prpara, prpara, dir
  red_make_prpara, prpara, autocrop
  red_make_prpara, prpara, clips
  red_make_prpara, prpara, crop     
  red_make_prpara, prpara, integer  
  red_make_prpara, prpara, interactive
  red_make_prpara, prpara, limb_data 
  red_make_prpara, prpara, noaligncont 
  red_make_prpara, prpara, nocavitymap    
  red_make_prpara, prpara, overwrite
  red_make_prpara, prpara, smooth
  red_make_prpara, prpara, tiles
  red_make_prpara, prpara, tuning_selection

  ;; It doesn't make sense to make new stokes cubes unless we are also
  ;; prepared to overwrite the scan cube.
  if keyword_set(redemodulate) then overwrite = 1


  ;; Default keywords
  if n_elements(cmap_fwhm) eq 0 then fwhm = 7.0
  if n_elements(clips) eq 0 then clips = [8, 4,  2,  1,  1  ]
  if n_elements(tiles) eq 0 then tiles = [8, 16, 32, 64, 128]

  if keyword_set(limb_data) then autocrop = 0

  ;; Output directory
  if(n_elements(odir) eq 0) then odir = self.out_dir + '/cubes_scan/' 

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
  
  ;; We do currently correct for the small scale cavity map in CRISP
  ;; data. (We should get this from earlier meta data!)
  remove_smallscale = 1
  
  ;; Get metadata from logfiles
  red_logdata, self.isodate, time_r0, r0 = metadata_r0
  red_logdata, self.isodate, time_pointing, diskpos = metadata_pointing, rsun = rsun

  ;; Search for restored images
  case self.filetype of
    'ANA': extension = '.f0'
    'MOMFBD': extension = '.momfbd'
    'FITS': extension = '.fits'
  endcase

  ;; Get the available scan numbers as cheaply as possible.
  wfiles = file_search(dir + wbdetector+'*' + extension, count = Nwfiles) 
  ;; The global WB file name should have no tuning info.
  indx = where(~strmatch(wfiles,'*_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_[+-]*'), Nmatch $
               , complement = complement, Ncomplement = Ncomplement)
  if Nmatch eq 0 then stop
  wbgfiles = wfiles[indx]
  self -> extractstates, wbgfiles, wbgstates
  wfiles = wfiles[complement]
  
  ;; Some info common to all scans
  prefilter = wbgstates[0].prefilter
  datestamp = fxpar(red_readhead(wbgfiles[0]), 'STARTOBS')
  timestamp = (strsplit(datestamp, 'T', /extract))[1]

  ;; Get a subset of the available scans, either through the scannos
  ;; keyword or by a selection dialogue.
  if ~(n_elements(scannos) gt 0 && scannos eq '*') then begin
    if n_elements(scannos) gt 0 then begin
      ;; Selected a subset through the scannos keyword
      uscans = red_expandrange(scannos)
      match2, uscans, wbgstates.scannumber, scanindx
      if max(scanindx eq -1) eq 1 then begin
        print, inam + ' : You asked for scans ' + scannos + '. However, scans ' $
               + red_collapserange(uscans[where(scanindx eq -1)], ld = '', rd = '') $
               + ' are not available.'
        print, 'Please change the scannos keyword to a subset of ' $
               + red_collapserange(wbgstates.scannumber, ld = '', rd = '') $
               + ' and try again.'
        retall
      endif
      Nscans = n_elements(scanindx)
    endif else begin
      ;; Selection dialogue
      selectionlist = strtrim(wbgstates[uniq(wbgstates.scannumber, sort(wbgstates.scannumber))].scannumber, 2)
      tmp = red_select_subset(selectionlist $
                              , qstring = inam + ' : Select scans:' $
                              , count = Nscans, indx = scanindx)
    endelse
    wbgstates = wbgstates[scanindx]
    wbgfiles  = wbgfiles[scanindx]
  endif
  uscans = wbgstates.scannumber
  nscans = n_elements(uscans)
  
  ;; wbgfiles & wbgstates -> selected and existing global WB files in
  ;; directory. Also uscans is the scans to process and Nscans is the
  ;; number of such scans.


  ;; Do we need this?
  ;; Establish the FOV, perhaps based on the selected wb files.
  red_bad_subfield_crop, wbgfiles, crop, autocrop = autocrop,  interactive = interactive
  hdr = red_readhead(wbgfiles[0])
  im_dim = fxpar(hdr, 'NAXIS*')
  x0 = crop[0]
  x1 = im_dim[0]-1 - crop[1]
  y0 = crop[2]
  y1 = im_dim[1]-1 - crop[3]
  Nx = x1 - x0 + 1
  Ny = y1 - y0 + 1
  

  
;  ;; Now let's limit the files & states arrays to only the
;  ;; scans to process.
;  self -> selectfiles, files = files, states = states $
;                       , scan = uscans, sel = pindx
;  files = files[pindx]
;  states = states[pindx]
;
;  ;; files & states -> existing files of selected scans in directory.
;
;  
;  ;; Select the nb and wb "per tuning files" by excluding the global
;  ;; WB images
;  self -> selectfiles, files = files, states = states $
;                       , cam = wbcamera, ustat = '' $
;                       , sel = wbgindx $
;;                       , count = Nscans $
;                       , complement = complement, Ncomplement = Ncomplement
;  ;; We have no special state (or absence of state) to identify
;  ;; the global WB images but we do know that their exposure times
;  ;; are much larger than the ones corresponding to the individual
;  ;; NB states.
;  wbindx = where(states.exposure gt mean(states.exposure)*1.5 $
;                 ,                 , complement = complement, Ncomplement = Ncomplement) 
;  ;; All the per-tuning files and states
;  pertuningfiles = files[complement]
;  pertuningstates = states[complement]
;
;; pertuningfiles & pertuningstates -> existing files of selected
;; scans in directory, excluding the global WB images.

  if ~keyword_set(nocavitymap) then begin

    ;; Read the original cavity map
    cfile = self.out_dir + 'flats/spectral_flats/' $
            + strjoin([nbtdetector $
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
            + strjoin([nbrdetector $
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

  
  for iscan = 0L, Nscans-1 do begin
    
    ;; This is the loop in which the cubes are made.
    
    ;; All momfbd output for this scan
    files = file_search(dir + '*_'+string(uscans[iscan], format = '(i05)')+'_*' + extension, count = Nfiles) 

    ;; Make a Stokes cube?
    red_extractstates, files, lc = lc
    ulc = lc[uniq(lc, sort(lc))]
    Nlc = n_elements(ulc)
    makestokes = (Nlc eq 4) and ~keyword_set(nopolarimetry)

    if makestokes then Nstokes = 4 else Nstokes = 1

    ;; Make output file name
    midpart = prefilter + '_' + datestamp + '_scan=' $ 
              + strtrim(uscans[iscan], 2)
    if makestokes then midpart += '_stokes'
    ofile = 'nb_'+midpart+'_corrected.fits'
    filename = odir+ofile

    ;; Already done?
    if file_test(filename) then begin
      if keyword_set(overwrite) then begin
        print, inam + ' : Overwriting existing data cube:'
        print, filename
      endif else begin
        print, inam + ' : Skipping existing cube, use /overwrite to override:'
        print, filename
        continue
      endelse
    endif

    ;; Make the directory if needed.
    file_mkdir, odir

    self -> selectfiles, files = files, states = states $
                         , cam = wbcamera, sel = windx, count = Nwb
    wfiles  = files[windx]
    wstates = states[windx]
    ;; The global WB file name should have no tuning info.
    indx = where(~strmatch(wfiles,'*_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_[+-]*'), Nmatch $
                 , complement = complement, Ncomplement = Ncomplement)

    ;; Read global WB file to add it as an extension to the cube. And
    ;; possibly use for WB correctionb.
    wbim = (red_readdata(wfiles[indx], head = wbhdr))[x0:x1, y0:y1]

    ;; The non-global WB files
    wfiles  = wfiles[complement]
    wstates = wstates[complement]
    Nwb = Ncomplement
    
    if makestokes then begin
      
      ;; Make Stokes cubes for this scan if needed
      self -> make_stokes_cubes, dir, uscans[iscan] $
                                 , clips = clips $
                                 , cmap_fwhm = cmap_fwhm $
                                 , nocavitymap = nocavitymap $
                                 , redemodulate = redemodulate $
                                 , smooth = smooth $
                                 , snames = snames $
                                 , stokesdir = stokesdir $
                                 , tiles = tiles 

      ;; The names of the Stokes cubes for this scan are returned in
      ;; snames. 
      self -> extractstates, snames, sstates
      Ntuning = n_elements(snames)
      indx = sort(sstates.tun_wavelength)
      snames = snames[indx]
      sstates = sstates[indx]

      ;; No wbstretchcorr is needed for Stokes data, done already when
      ;; demodulating.

      ;; Read a header to use as a starting point.
      nbhdr = red_readhead(snames[0])

      
    endif else begin

      ;; Select files from the two NB cameras

      self -> selectfiles, files = files, states = states $
                           , cam = nbtcamera $
                           , sel = nbtindx, count = Nnbt
      nbtstates = states[nbtindx]
      nbtfiles = files[nbtindx]
      
      self -> selectfiles, files = files, states = states $
                           , cam = nbrcamera $
                           , sel = nbrindx, count = Nnbr
      nbrstates = states[nbrindx]
      nbrfiles = files[nbrindx]
      
      if Nnbt eq Nnbr then Nnb = Nnbt else stop 
      
      ;; Do WB correction?
      if Nwb eq Nnb then wbstretchcorr = 1B else wbstretchcorr = 0B

      ;; Unique tuning states, sorted by wavelength
      utunindx = uniq(nbtstates.fpi_state, sort(nbtstates.fpi_state))
      Ntuning = n_elements(utunindx)
      sortindx = sort(nbtstates[utunindx].tun_wavelength)
      ufpi_states = nbtstates[utunindx[sortindx]].fpi_state
      utunwavelength = nbtstates[utunindx[sortindx]].tun_wavelength


;      ;; Per-tuning files, wb and nb, only for selected scans
;      self -> selectfiles, files = pertuningfiles, states = pertuningstates $
;                           , scan = uscans $
;                           , cam = wbcamera $
;                           , sel = wbindx, count = Nwb
;      wbstates = pertuningstates[wbindx]
;      wbfiles = pertuningfiles[wbindx]




;    ;; The NB files in this scan, sorted in tuning wavelength order.
;    self -> selectfiles, files = pertuningfiles, states = pertuningstates $
;                         , cam = nbcamera, scan = uscans[iscan] $
;                         , sel = scan_nbindx, count = Nnb
;    scan_nbfiles = pertuningfiles[scan_nbindx]
;    scan_nbstates = pertuningstates[scan_nbindx]
;    sortindx = sort(scan_nbstates.tun_wavelength)
;    scan_nbfiles = scan_nbfiles[sortindx]
;    scan_nbstates = scan_nbstates[sortindx]
;
;    ;; The WB files in this scan, sorted as the NB files
;    self -> selectfiles, files = pertuningfiles, states = pertuningstates $
;                         , cam = wbcamera, scan = uscans[iscan] $
;                         , sel = scan_wbindx, count = Nwb
;    scan_wbfiles = pertuningfiles[scan_wbindx]
;    scan_wbstates = pertuningstates[scan_wbindx]
;    match2, scan_nbstates.fpi_state, scan_wbstates.fpi_state, sortindx
;    scan_wbfiles = scan_wbfiles[sortindx]
;    scan_wbstates = scan_wbstates[sortindx]


      nbhdr = red_readhead(nbrfiles[0]) ; Use for main header
      
    
      ;; Read global WB file to use as reference when destretching
      ;; per-tuning wb files and then the corresponding nb files. Plus
      ;; add it as an extension to the cube.
;      wbim = (red_readdata(wfiles[iscan], head = wbhdr))[x0:x1, y0:y1]
      
    endelse
    

    ;; Load prefilters
    
    ;; Crisp-T

    pfile = self.out_dir + '/prefilter_fits/Crisp-T_'+prefilter+'_prefilter.idlsave'
    if ~file_test(pfile) then begin
      print, inam + ' : prefilter file not found: '+pfile
      return
    endif
    restore, pfile              ; Restores variable prf which is a struct

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
    restore, pfile              ; Restores variable prf which is a struct

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

    ;; Make FITS header for the output cube
    hdr = nbhdr
    red_fitsdelkeyword, hdr, 'STATE' ; Not a single state for cube 
    red_fitsaddkeyword, hdr, 'BITPIX', -32
    ;; Add info about this step
    prstep = 'Prepare NB science data cube'

    self -> headerinfo_addstep, hdr $
                                , prstep = prstep $
                                , prpara = prpara $
                                , prproc = inam

    ;; WB and NB data come from different cameras.
    red_fitsaddkeyword, hdr, 'CAMERA', nbtcamera + ',' + nbrcamera
    red_fitsaddkeyword, hdr, 'DETECTOR',  nbtdetector + ',' + nbrdetector

    anchor = 'DATE'



    ;; Add some keywords
    red_fitsaddkeyword, anchor = anchor, hdr, 'OBS_HDU', 1
    
    ;; Add info to headers
    red_fitsaddkeyword, anchor = anchor, hdr, 'BUNIT', units, 'Units in array'
    red_fitsaddkeyword, anchor = anchor, hdr, 'BTYPE', 'Intensity', 'Type of data in array'

    ;; Initialize fits file, set up for writing the data part.
    dims = [Nx, Ny, Ntuning, Nstokes, 1] 
    self -> fitscube_initialize, filename, hdr, lun, fileassoc, dims 


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


    

    for ituning = 0L, Ntuning - 1 do begin 

;      state = ufpi_states[ituning]


      red_progressbar, ituning, Ntuning $
                       , /predict $
                       , 'Assembling the cube'


      if keyword_set(makestokes) then begin
        
        nbims = red_readdata(snames[ituning], head = nbhead)

        ;; Wavelength 
        wcs[ituning].wave = sstates[ituning].tun_wavelength*1d9

        ;; Exposure time
        exp_array[ituning]  = fxpar(nbhead, 'XPOSURE')
        sexp_array[ituning] = fxpar(nbhead, 'TEXPOSUR')
        nsum_array[ituning] = fxpar(nbhead, 'NSUMEXP')

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

        for istokes = 0, Nstokes-1 do begin
          self -> fitscube_addframe, fileassoc, nbims[*, *, 0, istokes] $
                                     , ituning = ituning, istokes = istokes
        endfor                  ; istokes


      endif else begin

        ;; Read and add multiple files per tuning, align them by use
        ;; of extra WB objects, sum them, and add the result to the
        ;; cube.  

        ;; Select files for this tuning, sort them in LC order
        self -> selectfiles, files = nbtfiles, states = nbtstates $
                             , fpi_state = ufpi_states[ituning] $
                             , sel = nbtindx
        tun_nbtfiles  = nbtfiles[nbtindx]
        tun_nbtstates = nbtstates[nbtindx]
        indx = sort(tun_nbtstates.lc)
        tun_nbtfiles  = tun_nbtfiles[indx]
        tun_nbtstates = tun_nbtstates[indx]
        
        self -> selectfiles, files = nbrfiles, states = nbrstates $
                             , fpi_state = ufpi_states[ituning] $
                             , sel = nbrindx
        tun_nbrfiles  = nbrfiles[nbrindx]
        tun_nbrstates = nbrstates[nbrindx]
        indx = sort(tun_nbrstates.lc)
        tun_nbrfiles  = tun_nbrfiles[indx]
        tun_nbrstates = tun_nbrstates[indx]

        self -> selectfiles, files = wfiles, states = wstates $
                             , fpi_state = ufpi_states[ituning] $
                             , sel = wbindx
        tun_wfiles  = wfiles[wbindx]
        tun_wstates = wstates[wbindx]
        indx = sort(tun_wstates.lc)
        tun_wfiles  = tun_wfiles[indx]
        tun_wstates = tun_wstates[indx]

        Nexposures = n_elements(wbindx)
        
        ;; Loop over the exposures
        nbim = 0.0
        undefine, xps, texps, nsums
        for iexposure = 0, Nexposures-1 do begin

          ;; Read images
          nbtim = (red_readdata(tun_nbtfiles[iexposure], head = nbthdr))[x0:x1, y0:y1]
          nbrim = (red_readdata(tun_nbrfiles[iexposure], head = nbrhdr))[x0:x1, y0:y1]

          ;; Apply prefilter curve
          nbtim *= nbt_rpref[ituning]
          nbrim *= nbr_rpref[ituning]
          
          if wbstretchcorr then begin
            wim = (red_readdata(tun_wfiles[iexposure], head = whdr))[x0:x1, y0:y1]
            grid1 = red_dsgridnest(wbim, wim, tiles, clips)
            nbtim = red_stretch(temporary(nbtim), grid1)
            nbrim = red_stretch(temporary(nbrim), grid1)
          endif
          nbim += (nbtim + nbrim) / 2.

          red_fitspar_getdates, nbthdr $
                                , date_beg = date_beg $
                                , date_avg = date_avg $
                                , date_end = date_end 
          red_append, tbegs, red_time2double((strsplit(date_beg,'T',/extract))[1])
          red_append, tavgs, red_time2double((strsplit(date_avg,'T',/extract))[1])
          red_append, tends, red_time2double((strsplit(date_end,'T',/extract))[1])

          ;; These done twice because of the two NB cameras
          red_append, xps,   fxpar(nbthdr, 'XPOSURE')
          red_append, texps, fxpar(nbthdr, 'TEXPOSUR')
          red_append, nsums, fxpar(nbthdr, 'NSUMEXP')
          red_append, xps,   fxpar(nbthdr, 'XPOSURE')
          red_append, texps, fxpar(nbthdr, 'TEXPOSUR')
          red_append, nsums, fxpar(nbthdr, 'NSUMEXP')
          
        endfor                  ; ifile
        nbim /= Nexposures
        
        tbeg_array[ituning] = min(tbegs)
        tavg_array[ituning] = mean(tavgs)
        tend_array[ituning] = max(tends)

        date_beg_array[ituning] = self.isodate + 'T' + red_timestring(tbeg_array[ituning])
        date_avg_array[ituning] = self.isodate + 'T' + red_timestring(tavg_array[ituning])
        date_end_array[ituning] = self.isodate + 'T' + red_timestring(tend_array[ituning])

        ;; Exposure time
        exp_array[ituning]  = total(xps)
        sexp_array[ituning] = mean(texps)
        nsum_array[ituning] = round(total(nsums))

        ;; Wavelength 
        wcs[ituning].wave = tun_nbtstates[0].tun_wavelength*1d9

        ;; Exposure time
        exp_array[ituning]  = total(xps)
        sexp_array[ituning] = mean(texps)
        nsum_array[ituning] = round(total(nsums))

        self -> fitscube_addframe, fileassoc, temporary(nbim) $
                                   , ituning = ituning

      endelse 
      

      ;; Time
      wcs[ituning].time = tavg_array[ituning, 0]

    endfor                      ; ituning

    
    ;; Get pointing at center of FOV for the different tunings. (Or
    ;; get from Stokes cube?)
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
    if ~keyword_set(nocavitymap) then red_fitscube_addcmap, filename, reform(cmap1, Nx, Ny, 1, 1, 1)

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



    if makestokes && ~keyword_set(nocrosstalk) then begin

      ;; Correct the cube for cross-talk, I --> Q,U,V.
      ;; Include tuning_selection in keywords to allow it to be
      ;; specified when calling make_scan_cube and also, if not
      ;; specified, to use the same tuning_selection for multiple
      ;; cubes. Also, including mag_mask will reuse the same
      ;; magnetic-features mask.
      self -> fitscube_crosstalk, filename $
                                  , mag_mask = mag_mask $
                                  , /nostatistics $
                                  , tuning_selection = tuning_selection

    endif

    ;; Include the global WB image as an image extension
    ehdr=wbhdr
    fxaddpar, ehdr, 'XTENSION', 'IMAGE'
    sxdelpar, ehdr, 'SIMPLE'
    check_fits, wbim, ehdr, /update
    fxaddpar, ehdr, 'DATE', red_timestamp(/utc, /iso)
    anchor = 'DATE'
    red_fitsaddkeyword, anchor = anchor, ehdr, 'EXTNAME', 'WBIMAGE', 'Wideband image'
    red_fitsaddkeyword, anchor = anchor, ehdr, 'PCOUNT', 0
    red_fitsaddkeyword, anchor = anchor, ehdr, 'GCOUNT', 1
    writefits, filename, wbim, ehdr, /append

    if ~keyword_set(nowbintensitycorr) then begin
      ;; Correct intensity with respect to solar elevation and
      ;; exposure time.
      self -> fitscube_intensitycorr, filename
    endif

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
      red_fitscube_statistics, filename, /write
    endif
    
    ;; Done with this scan.
    print, inam + ' : Narrowband scan cube stored in:'
    print, filename
    
  endfor                        ; iscan

  
end
