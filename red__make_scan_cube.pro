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
;     circular_fov : in, optional, type=boolean
;
;       Do not make room for corners of the FOV in derotation.
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
;     direction : in, optional, type=integer, default="from class object"
;
;       The relative orientation of reference cameras of different
;       instruments.
;
;     integer : in, optional, type=boolean
;
;       Store as integers instead of floats.
;
;     intensitycorrmethod : in, optional, type="string or boolean", default='fit'
;
;       Indicate whether to do intensity correction based on WB data
;       and with what method. See documentation for red::fitscube_intensitycorr.
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
;     norotation : in, optional, type=boolean
;
;       Do not apply the direction parameter.
;
;     odir : in, optional, type=string, detault='cubes_scan/'
;
;       The output directory.
;
;     overwrite : in, optional, type=boolean
;
;       Don't care if cube is already on disk, overwrite it
;       with a new version.
;
;     padmargin : in, optional, type=integer, default=40
;
;       Amount of rotation padding in pixels.
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
;     tiles : in, optional, type=array
;
;       Used to compute stretch vectors for the wideband alignment. 
;
;     tuning_selection : in, optional, type="integer or array"
;
;       The index or indices in the tuning dimension to use for
;       calculating the polarimetric crosstalk correction. Should
;       correspond to continuum (or as close to continuum as
;       possible), where the polarimetric signal is minimal. Negative
;       indices are allowed.
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
;    2020-01-16 : MGL. New keywords intensitycorrmethod and odir.
; 
;    2020-04-01 : MGL. New keywords direction and norotation.
; 
;    2020-04-27 : MGL. New keyword rotation.
; 
;    2020-07-15 : MGL. Remove keyword smooth.
; 
;    2020-10-28 : MGL. Remove statistics calculations.
;
;    2020-12-15 : JdlCR. There is one element missing in dummy ff
;                 vector. /norotation keyword was not implemented
;                 properly. The cropping info was not applied to the
;                 output image.
;
;    2022-06-24 : JdlCR. Bugfix, the cmap was always rotated
;                 regardless of the /norotation keyword.
;
;    2022-09-04 : MGL. CRISP --> RED. New keyword circular_fov.
;
;    2025-02-20 : MGL. Adapt to new camera alignment model.
; 
;    2025-10-09 : MGL. New keyword no_intensitycorr_timecheck.
;
;-
pro red::make_scan_cube, dir $
                         , autocrop = autocrop $
                         , circular_fov = circular_fov $ 
                         , clips = clips $
                         , cmap_fwhm = cmap_fwhm $
                         , crop = crop $
                         , direction = direction $
                         , filename = filename $
                         , fitpref_time = fitpref_time $
                         , integer = integer $
                         , intensitycorrmethod = intensitycorrmethod $
                         , interactive = interactive $
                         , limb_data = limb_data $
                         , mosaic = mosaic $
                         , nocavitymap = nocavitymap $
                         , nocrosstalk = nocrosstalk $
                         , no_intensitycorr_timecheck = no_intensitycorr_timecheck $
                         , nomissing_nans = nomissing_nans $
                         , nopolarimetry = nopolarimetry $
                         , norotation = norotation $
                         , nthreads = nthreads $
                         , odir = odir $
                         , overwrite = overwrite $
                         , padmargin = padmargin $
                         , redemodulate = redemodulate $
                         , rotation = rotation $
                         , scannos = scannos $
                         , tiles = tiles $
                         , tuning_selection = tuning_selection
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  if n_elements(mosaic) gt 0 and n_elements(scannos) eq 0 then scannos = '0'
  if n_elements(direction) eq 0 then direction = self.direction
  if n_elements(rotation)  eq 0 then rotation  = self.rotation
  if n_elements(padmargin) eq 0 then padmargin  = 40

  ;; Make prpara
  red_make_prpara, prpara, autocrop
  red_make_prpara, prpara, clips
  red_make_prpara, prpara, cmap_fwhm
  red_make_prpara, prpara, crop
  red_make_prpara, prpara, dir
  red_make_prpara, prpara, direction 
  red_make_prpara, prpara, integer  
  red_make_prpara, prpara, intensitycorrmethod
  red_make_prpara, prpara, interactive
  red_make_prpara, prpara, limb_data 
  red_make_prpara, prpara, noaligncont 
  red_make_prpara, prpara, nocavitymap  
  red_make_prpara, prpara, nocrosstalk   
  red_make_prpara, prpara, norotation    
  red_make_prpara, prpara, overwrite
  red_make_prpara, prpara, padmargin 
  red_make_prpara, prpara, rotation    
  red_make_prpara, prpara, tiles
  red_make_prpara, prpara, tuning_selection

  ;; It doesn't make sense to make new stokes cubes unless we are also
  ;; prepared to overwrite the scan cube.
  if keyword_set(redemodulate) then overwrite = 1


  ;; Default keywords
  if n_elements(cmap_fwhm) eq 0 then cmap_fwhm = 7.0
  if n_elements(clips) eq 0 then clips = [8, 4,  2,  1,  1  ]
  if n_elements(tiles) eq 0 then tiles = [8, 16, 32, 64, 128]

  if keyword_set(limb_data) then autocrop = 0

  ;; Output directory
  if(n_elements(odir) eq 0) then odir = self.out_dir + '/cubes_scan/' 

  ;; Camera/detector identification
  self->getdetectors
  wbindx      = where(strmatch(*self.cameras,'*-W'))
  wbcamera    = (*self.cameras)[wbindx[0]]
  wbdetector  = (*self.detectors)[wbindx[0]]
  nbindx      = where(strmatch(*self.cameras,'*-[NTR]')) 
  nbcameras   = (*self.cameras)[nbindx]
  nbdetectors = (*self.detectors)[nbindx]
  Nnbcams     = n_elements(nbcameras)

  instrument = (strsplit(wbcamera, '-', /extract))[0]

  polarimetric_data = self -> polarimetric_data()

;  wbindx     = where(strmatch(*self.cameras,'Crisp-W'))
;  wbcamera   = (*self.cameras)[wbindx[0]]
;  wbdetector = (*self.detectors)[wbindx[0]]
;  nbtindx     = where(strmatch(*self.cameras,'Crisp-T')) 
;  nbtcamera   = (*self.cameras)[nbtindx[0]]
;  nbtdetector = (*self.detectors)[nbtindx[0]]
;  nbrindx     = where(strmatch(*self.cameras,'Crisp-R')) 
;  nbrcamera   = (*self.cameras)[nbrindx[0]]
;  nbrdetector = (*self.detectors)[nbrindx[0]]

  ;; How to handle small scale variations in cavity maps.
  case instrument of
    'Crisp' : begin
      ;; We do correct for the small scale cavity map in CRISP data.
      ;; (We should get this from earlier meta data?)
      remove_smallscale = 1
    end
    'Chromis' : begin
      remove_smallscale = 0
    end
    'Crisp2' : begin
      remove_smallscale = 0
    end
  endcase

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
  wfiles = file_search(dir + '/' + wbdetector+'*' + extension, count = Nwfiles) 
  ;; The global WB file name should have no tuning info.
  indx = where(~strmatch(wfiles,'*_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]_[+-]*'), Nmatch $
               , complement = complement, Ncomplement = Ncomplement)
  if Nmatch eq 0 then stop
  wbgfiles = wfiles[indx]
  self -> extractstates, wbgfiles, wbgstates
  wfiles = wfiles[complement]


  
  ;; Some info common to all scans
  prefilter = wbgstates[0].prefilter
  wbghdr = red_readhead(wbgfiles[0])
  datestamp = fxpar(wbghdr, 'STARTOBS')
  timestamp = (strsplit(datestamp, 'T', /extract))[1]

  ;; Get the alignment mode, projective or polywarp.
  alignment_model = red_align_model(wbgfiles[0])
  
  ;; Get the image scale from the header
  image_scale = float(fxpar(wbghdr, 'CDELT1'))
  
  if ~keyword_set(fitpref_time) then begin
    fitpref_t='_'
    dt = strtrim(fxpar(wbghdr, 'DATE-AVG'), 2)
    avg_ts = (strsplit(dt, 'T', /extract))[1]
    avg_time = red_time2double(avg_ts)
    pfls = file_search(self.out_dir + '/prefilter_fits/*-[NT]_'+prefilter+ $
                       '_[0-9][0-9]:[0-9][0-9]:[0-9][0-9]*save', count=Npfls)
    if Npfls gt 0 then begin
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
  
  red_fitspar_getdates, wbghdr $
                        , date_beg = date_beg $
                        , date_end = date_end $
                        , date_avg = date_avg $
                        , count_avg = hasdateavg $
                        , comment_avg = comment_avg
  
  if hasdateavg then begin
    date_avg_split = strsplit(date_avg, 'T', /extract, count = Nsplit)
    ddate = date_avg_split[0]
    if Nsplit gt 1 then time = date_avg_split[1] else undefine, time
  endif else undefine, ddate, time

  ;; Derotation angle
  
  ang = (red_lp_angles(time, ddate[0], /from_log, offset_angle = rotation))[0]
  
  ;; Get a subset of the available scans, either through the scannos
  ;; keyword or by a selection dialogue.
  if ~(n_elements(scannos) gt 0 && string(scannos) eq '*') then begin
    if n_elements(scannos) gt 0 then begin
      ;; Selected a subset through the scannos keyword
      uscans = red_expandrange(string(scannos))
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
  Nscans = n_elements(uscans)
  
  ;; wbgfiles & wbgstates -> selected and existing global WB files in
  ;; directory. Also uscans is the scans to process and Nscans is the
  ;; number of such scans.


  ;; Do we need this?
  ;; Establish the FOV, perhaps based on the selected wb files.
  x01y01 = red_bad_subfield_crop(wfiles, crop $
                         , autocrop = autocrop  $
                         , direction = direction $
                         , interactive = interactive)
  x0 = x01y01[0] & x1 = x01y01[1] & y0 = x01y01[2] & y1 = x01y01[3]
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

;  if ~keyword_set(nocavitymap) then begin
;
;    ;; Read the original cavity map
;    cfile = self.out_dir + 'flats/spectral_flats/' $
;            + strjoin([nbtdetector $
;                       , prefilter $
;                       , 'fit_results.sav'] $
;                      , '_')
;
;    if ~file_test(cfile) then begin
;      print, inam + ' : Error, calibration file not found -> '+cfile
;      stop
;    endif
;    restore, cfile                  ; The cavity map is in a struct called "fit". 
;    cmapt = reform(fit.pars[1,*,*]) ; Unit is [Angstrom]
;    cmapt /= 10.                    ; Make it [nm]
;    cmapt = -cmapt                  ; Change sign so lambda_correct = lambda + cmap
;    fit = 0B                        ; Don't need the fit struct anymore.
;    
;    cfile = self.out_dir + 'flats/spectral_flats/' $
;            + strjoin([nbrdetector $
;                       , prefilter $
;                       , 'fit_results.sav'] $
;                      , '_')
;
;    if ~file_test(cfile) then begin
;      print, inam + ' : Error, calibration file not found -> '+cfile
;      stop
;    endif
;    restore, cfile                  ; The cavity map is in a struct called "fit". 
;    cmapr = reform(fit.pars[1,*,*]) ; Unit is [Angstrom]
;    cmapr /= 10.                    ; Make it [nm]
;    cmapr = -cmapr                  ; Change sign so lambda_correct = lambda + cmap
;    fit = 0B                        ; Don't need the fit struct anymore.
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
;    cmap1t = red_apply_camera_alignment(cmap1t, alignment_model, instrument+'-T' $
;                                        , pref = prefilter $
;                                        , amap = amapt $
;                                        , /preserve_size)
;    cmap1r = red_apply_camera_alignment(cmap1r, alignment_model, instrument+'-R' $
;                                        , pref = prefilter $
;                                        , amap = amapr $
;                                        , /preserve_size)
;    
;    ;; At this point, the individual cavity maps should be corrected
;    ;; for camera misalignments, so they should be aligned with
;    ;; respect to the cavity errors on the etalons. So we can sum
;    ;; them.
;    cmap1 = (cmap1r + cmap1t) / 2.
;
;    if self.filetype eq 'MOMFBD' then begin
;      ;; Crop the cavity map to the FOV of the momfbd-restored images.
;      mr = momfbd_read(wbgfiles[0],/nam)
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
;  endif

  
  for iscan = 0L, Nscans-1 do begin
    
    ;; This is the loop in which the cubes are made.
    
    ;; All momfbd output for this scan
    files = file_search(dir + '/*_'+string(uscans[iscan], format = '(i05)')+'_*' + extension, count = Nfiles) 

    ;; Make a Stokes cube?
    red_extractstates, files, lc = lc
    ulc = lc[uniq(lc, sort(lc))]
    Nlc = n_elements(ulc)
    makestokes = (Nlc eq 4) and ~keyword_set(nopolarimetry)

    if makestokes then Nstokes = 4 else Nstokes = 1

    ;; Make output file name
    if n_elements(mosaic) ne 0 then $
      midpart = prefilter + '_' + datestamp + '_mos' + string(mosaic, format='(I02)') $
    else $
      midpart = prefilter + '_' + datestamp + '_scan=' + strtrim(uscans[iscan], 2)
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
    ;; possibly use for WB correction.
    wbim = (red_readdata(wfiles[indx], head = wbhdr, direction = direction))[x0:x1, y0:y1]

    ;; The non-global WB files
    wfiles  = wfiles[complement]
    wstates = wstates[complement]
    Nwb = Ncomplement

    spl = strsplit(wfiles[0],'/',/extract)
    cw = where(strmatch(spl,'*cfg*'))
    cfg_dir=strjoin(spl[0:cw],'/')
    if file_test(cfg_dir+'/fov_mask.fits') then begin    
      fov_mask = readfits(cfg_dir+'/fov_mask.fits')
      if self.filetype eq 'MOMFBD' then begin
        mr = momfbd_read(wfiles[0], /nam)
        fov_mask = red_crop_as_momfbd(fov_mask, mr)
      endif else begin          ; get cropping from cfg file      
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
        fov_mask = fov_mask[xx0:xx1,yy0:yy1]
      endelse
      fov_mask = red_rotate(fov_mask, direction)
    endif

    if n_elements(fov_mask) gt 0 then wbim *= fov_mask    
    
    if makestokes then begin
      
      ;; Make Stokes cubes for this scan if needed
      self -> make_stokes_cubes, dir, uscans[iscan] $
         , clips = clips $
         , cmap_fwhm = cmap_fwhm $
         , /nocavitymap $       ; Cavity maps in Stokes cubes aren't used for anything
         , /notelmat $          ; Until Pit has had time to measure them
         , redemodulate = redemodulate $
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

      utunindx = uniq(sstates.fpi_state, sort(sstates.fpi_state))
      Ntuning = n_elements(utunindx)
      sortindx = sort(sstates[utunindx].tun_wavelength)
      ufpi_states = sstates[utunindx[sortindx]].fpi_state
      utunwavelength = sstates[utunindx[sortindx]].tun_wavelength

      ;; Load prefilters
      self -> get_prefilterfit, sstates $ ;unbpertuningstates $
         , units = units $
         , prefilter_curve = prefilter_curve $
         , prf = prf $
         , wave_shifts = wave_shift

    endif else begin

      ;; Select files from the two NB cameras

      self -> selectfiles, files = files, states = states $
                                   , cam = nbcameras $
                                   , sel = nbindx, count = Nnb
      nbstates = states[nbindx]
      nbfiles = files[nbindx]
      
      ;; Do WB correction?
;      if Nwb eq Nnb then wbstretchcorr = 1B else wbstretchcorr = 0B
      wbstretchcorr = 1B
      
      ;; Unique tuning states, sorted by wavelength
      utunindx = uniq(nbstates.fpi_state, sort(nbstates.fpi_state))
      Ntuning = n_elements(utunindx)
      sortindx = sort(nbstates[utunindx].tun_wavelength)
      ufpi_states = nbstates[utunindx[sortindx]].fpi_state
      utunwavelength = nbstates[utunindx[sortindx]].tun_wavelength

      ;; Load prefilters
      self -> get_prefilterfit, nbstates[utunindx[sortindx]] $
         , units = units $
         , prefilter_curve = prefilter_curve $
         , prf = prf $
         , wave_shifts = wave_shift

      rpref = 1.d0/prefilter_curve
      
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


      nbhdr = red_readhead(nbfiles[0]) ; Use for main header
      
    
      ;; Read global WB file to use as reference when destretching
      ;; per-tuning wb files and then the corresponding nb files. Plus
      ;; add it as an extension to the cube.
;      wbim = (red_readdata(wfiles[iscan], head = wbhdr))[x0:x1, y0:y1]

    endelse
    
    if ~keyword_set(norotation) then begin
      if keyword_set(circular_fov) then begin
        ;; Red_rotation.pro only uses keyword full (which is set to ff
        ;; when calling from make_*_cube) if it is an array with at least
        ;; 5 elements. So setting it to -1 is the same as letting it stay
        ;; undefined, but it can still be passed on to make_nb_cube.
        ff = -1    
        ff = [0, -padmargin, padmargin, -padmargin, padmargin, 0]
        ;; Possibly change this to take the shifts into account but not
        ;; the angles. Something like ff = [0, mdx0, mdx1, mdy0, mdy1].
      endif else begin
        ff = [abs(ang),0,0,0,0,ang]
      endelse
      red_missing, wbim $
                   , nmissing = Nmissing, indx_missing = indx_missing, indx_data = indx_data
      if Nmissing gt 0 then begin
        bg = wbim[indx_missing[0]]
      endif else begin
        bg = median(wbim)
      endelse
      wbim_rot = red_rotation(wbim, ang, full = ff $
                              , background = bg $
                              , nthreads=nthreads)
      dims = [size(wbim_rot, /dim), Ntuning, Nstokes, 1] 
    endif else dims = [size(wbim, /dim), Ntuning, Nstokes, 1] 


;    ;; Crisp-T
;
;    pfile = self.out_dir + '/prefilter_fits/Crisp-T_'+prefilter+fitpref_t+'prefilter.idlsave'
;    if ~file_test(pfile) then begin
;      print, inam + ' : prefilter file not found: '+pfile
;      return
;    endif
;    restore, pfile              ; Restores variable prf which is a struct
;
;    wave_shift = prf.fitpars[1]/10. ; [nm] Shift the wavelengths by this amount
;    
;    nbt_units = prf.units
;;    nbt_prefilter_curve = prf.pref
;    nbt_prefilter_curve = red_intepf(prf.wav, prf.pref, utunwavelength*1.d10)
;;    nbt_prefilter_wav = prf.wav
;;    nbt_prefilter_wb = prf.wbint
;
;    nbt_rpref = 1.d0/nbt_prefilter_curve
;
;    ;; Crisp-R
;
;    pfile = self.out_dir + '/prefilter_fits/Crisp-R_'+prefilter+fitpref_t+'prefilter.idlsave'
;    if ~file_test(pfile) then begin
;      print, inam + ' : prefilter file not found: '+pfile
;      return
;    endif
;    restore, pfile              ; Restores variable prf which is a struct
;
;    nbr_units = prf.units  
;;    nbr_prefilter_curve = prf.pref
;    nbr_prefilter_curve = red_intepf(prf.wav, prf.pref, utunwavelength*1.d10)
;    
;;    nbr_prefilter_wav = prf.wav
;;    nbr_prefilter_wb = prf.wbint
;       
;    nbr_rpref = 1.d0/nbr_prefilter_curve
;
;    prefilter_curve = (nbt_prefilter_curve + nbr_prefilter_curve)/2.
;    
;    if nbr_units ne nbt_units then begin
;      print, inam + ' : Units for Crisp-T and Crisp-R do not match.'
;      print, inam + ' : Please rerun the prefilterfit step for these data.'
;      retall
;    endif
;    units = nbr_units

    ;; Make FITS header for the output cube
    hdr = nbhdr
    red_fitsdelkeyword, hdr, 'STATE' ; Not a single state for cube 
    red_fitsaddkeyword, hdr, 'BITPIX', -32
    ;; Add info about this step
    prstep = 'CONCATENATION,SPATIAL-ALIGNMENT,DESTRETCHING'

    self -> headerinfo_addstep, hdr $
                                , prstep = prstep $
                                , prpara = prpara $
                                , prproc = inam
    
    self -> headerinfo_addstep, hdr $
                                , prstep = 'CALIBRATION-INTENSITY-SPECTRAL' $
                                , prpara = prpara $
                                , prref = ['Hamburg FTS spectral atlas (Neckel 1999)' $
                                           , 'Calibration data from '+red_timestring(prf.time_avg, n = 0)] $
                                , prproc = inam

    ;; Now do the cavity maps
    
;    cavitymaps = fltarr(Nx, Ny, 1, 1, 1)
    cavitymaps = fltarr([dims[0:1], 1, 1, 1])

    if ~keyword_set(nocavitymap) then begin

      ;; Read the original cavity map

      if keyword_set(makestokes) then begin
        pindx = where(sstates.prefilter ne '3999') ; No cavity map for the Ca II H continuum
        pindx = pindx[uniq(sstates[pindx].prefilter, sort(sstates[pindx].prefilter))]
        cprefs = sstates[pindx].prefilter
      endif else begin
        pindx = where(nbstates.prefilter ne '3999') ; No cavity map for the Ca II H continuum
        pindx = pindx[uniq(nbstates[pindx].prefilter, sort(nbstates[pindx].prefilter))]
        cprefs = nbstates[pindx].prefilter
      endelse
      Ncprefs = n_elements(cprefs)

      ;; If multiple prefilters during the scan (as for CHROMIS Ca II),
      ;; calculate cmaps independently for them.
      for icprefs = 0, Ncprefs-1 do begin

        ;; Calculate average cmap if multiple cameras
        cmap1 = 0.0
        for icam = 0, Nnbcams-1 do begin

          self -> get_spectral_flats_info, cmap = cmap $
                                                  , detector = nbdetectors[icam] $
                                                  , pref = cprefs[icprefs]
          
          if keyword_set(remove_smallscale) then begin
            ;; If the small scale is already corrected, then include only the
            ;; low-resolution component in the metadata. The blurring kernel
            ;; should match how the low resolution component was removed when
            ;; making flats.
            npix = 30           ; Can we get this parameter from earlier headers?
            cpsf = red_get_psf(npix*2-1,npix*2-1,double(npix),double(npix))
            cpsf /= total(cpsf, /double)
            cmap = red_convolve(temporary(cmap), cpsf)
            this_cmap = cmap
          endif else begin
            ;; If the small scale is not already corrected, then we still want
            ;; to blur the cavity map slightly.
            npsf = round(cmap_fwhm * 7.)
            if((npsf/2)*2 eq npsf) then npsf += 1L
            psf = red_get_psf(npsf, npsf, cmap_fwhm, cmap_fwhm)
            psf /= total(psf, /double)
            ;; Leave the orignal cmap alone, we might need it later.
            this_cmap = red_convolve(cmap, psf)
          endelse

          ;; Apply the geometrical mapping      
          this_cmap = red_apply_camera_alignment(this_cmap $
                                                 , alignment_model, nbcameras[icam] $
                                                 , pref = cprefs[icprefs] $
                                                 , amap = amapt $
                                                 , /preserve_size)

          ;; The orientation should now be the same for multiple cameras
          cmap1 += this_cmap / Nnbcams
          
        endfor                  ; icam
        
        if self.filetype eq 'MOMFBD' then begin
          ;; Crop the cavity map to the FOV of the momfbd-restored images.
          if n_elements(mr) eq 0 then mr = momfbd_read(wfiles[0], /nam)
          cmap1 = red_crop_as_momfbd(cmap1, mr)
        endif else begin ;; Crop with information from the cfg file
          spl = strsplit(wbgfiles[0],'/',/extract)
          cw = where(strmatch(spl,'*cfg*'))
          cfg_dir=strjoin(spl[0:cw],'/')
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
          cmap1 = cmap1[xx0:xx1,yy0:yy1]
        endelse

        ;; Rotate to orientation of the WB camera
        cmap1 = red_rotate(cmap1, direction)

        ;; Clip to the selected FOV
        cmap1 = cmap1[x0:x1,y0:y1]
;        cmap11 = red_rotation(cmap1, ang, full = ff $
;                              , background = bg $
;                              , nthreads=nthreads)
;        
;        cavitymaps[0, 0, 0, 0, 0] = cmap11
;        cavitymaps = reform(cavitymaps,[dims[0:1],1,1,1])
;
        ;; For which tunings were this prefilter used?
;      tindx = where(scan_nbstates.prefilter eq cprefs[icprefs], Nt)
;        tindx = where(unbpertuningstates.prefilter eq cprefs[icprefs], Nt)

;        tindx = indgen(Ntuning)
;        Nt = n_elements(tindx)
;
;        
;        ;; Add cavity maps as WAVE distortions
;        if Nt gt 0 then begin
;          red_fitscube_addcmap, filename, cavitymaps $
;                                , cmap_number = icprefs+1 $
;                                , prefilter = cprefs[icprefs] $
;                                , indx = tindx
;        endif
        
      endfor                    ; icprefs
      
    endif


    
    ;; WB and NB data come from different cameras.
    red_fitsaddkeyword, hdr, 'CAMERA', strjoin(nbcameras,',')
    red_fitsaddkeyword, hdr, 'DETECTOR', strjoin(nbdetectors,',')

    anchor = 'DATE'

    ;; Add some keywords
    red_fitsaddkeyword, anchor = anchor, hdr, 'OBS_HDU', 1
    
    ;; POINT_ID (default - to be set manually later in case of grouping of files)
    date_obs = fxpar(hdr, 'DATE-OBS', count = Ndate_obs)  
    if Ndate_obs ne 0 then red_fitsaddkeyword, anchor = anchor, hdr, 'POINT_ID', date_obs else stop

    ;; Add info to headers
    red_fitsaddkeyword, anchor = anchor, hdr, 'BUNIT', units, 'Units in array'
    red_fitsaddkeyword, anchor = anchor, hdr, 'BTYPE', 'Intensity', 'Type of data in array'
    
    wcs = replicate({  wave:dblarr(2,2) $
                       , hplt:dblarr(2,2) $
                       , hpln:dblarr(2,2) $
                       , time:dblarr(2,2) $
                    }, Ntuning)

    ;; Initialize fits file, set up for writing the data part.
    self -> fitscube_initialize, filename, hdr, lun, fileassoc, dims, wcs = wcs


    ;; Set up for collecting time and wavelength data
    tbeg_array     = dblarr(Ntuning) ; Time beginning for state
    tend_array     = dblarr(Ntuning) ; Time end for state
    tavg_array     = dblarr(Ntuning) ; Time average for state
    date_beg_array = strarr(Ntuning) ; DATE-BEG for state
    date_end_array = strarr(Ntuning) ; DATE-END for state
    date_avg_array = strarr(Ntuning) ; DATE-AVG for state
    exp_array      = fltarr(Ntuning) ; Total exposure time
    sexp_array     = fltarr(Ntuning) ; Single exposure time
    nsum_array     = lonarr(Ntuning) ; Number of summed exposures


;    nbt_tscl = mean(nbt_prefilter_wb)
;    nbr_tscl = mean(nbr_prefilter_wb)
;    tscl = (nbt_tscl+nbr_tscl)/2.
    

    for ituning = 0L, Ntuning - 1 do begin 

;      state = ufpi_states[ituning]


      red_progressbar, ituning, Ntuning $
                       , /predict $
                       , 'Assembling the cube'


      if keyword_set(makestokes) then begin
        
        nbims = red_readdata(snames[ituning], head = nbhead, direction = direction) ;* tscl
        
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

          if n_elements(fov_mask) gt 0 then nbims[*, *, 0, istokes] = nbims[*, *, 0, istokes] * fov_mask

          if(~keyword_set(norotation)) then begin
            ;; bg = median(nbims[x0:x1, y0:y1, 0, istokes])
            red_missing, nbims[x0:x1, y0:y1, 0, istokes] $
                         , nmissing = Nmissing, indx_missing = indx_missing, indx_data = indx_data
            if Nmissing gt 0 then begin
              bg = (nbims[x0:x1, y0:y1, 0, istokes])[indx_missing[0]]
            endif else begin
              bg = median(nbims[x0:x1, y0:y1, 0, istokes])
            endelse
            red_fitscube_addframe, fileassoc $
                                   , red_rotation(nbims[x0:x1, y0:y1, 0, istokes], ang $
                                                  , background = bg $
                                                  , full = ff, nthreads=nthreads) $
                                   , ituning = ituning, istokes = istokes
          endif else begin
            red_fitscube_addframe, fileassoc $
                                   , nbims[x0:x1, y0:y1, 0, istokes] $
                                   , ituning = ituning, istokes = istokes
          endelse
        endfor                  ; istokes


      endif else begin

        ;; Read and add multiple files per tuning, align them by use
        ;; of extra WB objects, sum them, and add the result to the
        ;; cube.  

        ;; Select files for this tuning, sort them in LC order
        self -> selectfiles, files = nbfiles, states = nbstates $
                                     , fpi_state = ufpi_states[ituning] $
                                     , sel = nbindx
        tun_nbfiles  = nbfiles[nbindx]
        tun_nbstates = nbstates[nbindx]
        indx = sort(tun_nbstates.lc)
        tun_nbfiles  = tun_nbfiles[indx]
        tun_nbstates = tun_nbstates[indx]
        
;        self -> selectfiles, files = nbrfiles, states = nbrstates $
;                             , fpi_state = ufpi_states[ituning] $
;                             , sel = nbrindx
;        tun_nbrfiles  = nbrfiles[nbrindx]
;        tun_nbrstates = nbrstates[nbrindx]
;        indx = sort(tun_nbrstates.lc)
;        tun_nbrfiles  = tun_nbrfiles[indx]
;        tun_nbrstates = tun_nbrstates[indx]

        self -> selectfiles, files = wfiles, states = wstates $
                             , fpi_state = ufpi_states[ituning] $
                             , sel = wbindx
        tun_wfiles  = wfiles[wbindx]
        tun_wstates = wstates[wbindx]
        indx = sort(tun_wstates.lc)
        tun_wfiles  = tun_wfiles[indx]
        tun_wstates = tun_wstates[indx]

        Nexposures = n_elements(wbindx)

        ;; Why are we looping over exposures here? Would make more
        ;; sense to loop over LC states. 
        
        ;; Loop over the exposures
        nbim = 0.0
        undefine, xps, texps, nsums
        for iexposure = 0, Nexposures-1 do begin

          ;; Read images
          nbim = 0.0
          for i = 0, n_elements(tun_nbfiles)-1 do $
             nbim += (red_readdata(tun_nbfiles[i], head = nbhdr, direction = direction))[x0:x1, y0:y1] / n_elements(tun_nbfiles)
          
;          nbtim = (red_readdata(tun_nbtfiles[iexposure], head = nbthdr, direction = direction))[x0:x1, y0:y1]
;          nbrim = (red_readdata(tun_nbrfiles[iexposure], head = nbrhdr, direction = direction))[x0:x1, y0:y1]

          ;; Apply prefilter curve
          nbim *= rpref[ituning]
;          nbtim *= nbt_rpref[ituning] ;* nbt_tscl
;          nbrim *= nbr_rpref[ituning] ;* nbr_tscl

;          stop
          if wbstretchcorr then begin
            wim = (red_readdata(tun_wfiles[iexposure], head = whdr, direction = direction))[x0:x1, y0:y1]
            grid1 = red_dsgridnest(wbim, wim, tiles, clips)
;            wims = red_stretch(wim, grid1)
            nbim = red_stretch(temporary(nbim), grid1)
;            nbrim = red_stretch(temporary(nbrim), grid1)
          endif
;          nbim += (nbtim + nbrim) / 2.
          if n_elements(fov_mask) gt 0 then nbim *= fov_mask

          red_fitspar_getdates, nbhdr $
                                , date_beg = date_beg $
                                , date_avg = date_avg $
                                , date_end = date_end 
          red_append, tbegs, red_time2double((strsplit(date_beg,'T',/extract))[1])
          red_append, tavgs, red_time2double((strsplit(date_avg,'T',/extract))[1])
          red_append, tends, red_time2double((strsplit(date_end,'T',/extract))[1])

          ;; These done twice because of the two NB cameras
          red_append, xps,   fxpar(nbhdr, 'XPOSURE')
          red_append, texps, fxpar(nbhdr, 'TEXPOSUR')
          red_append, nsums, fxpar(nbhdr, 'NSUMEXP')
          red_append, xps,   fxpar(nbhdr, 'XPOSURE')
          red_append, texps, fxpar(nbhdr, 'TEXPOSUR')
          red_append, nsums, fxpar(nbhdr, 'NSUMEXP')
          
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
        wcs[ituning].wave = tun_nbstates[0].tun_wavelength*1d9

        ;; Exposure time
        exp_array[ituning]  = total(xps)
        sexp_array[ituning] = mean(texps)
        nsum_array[ituning] = round(total(nsums))

        red_missing, nbim, nmissing = Nmissing, indx_missing = indx_missing, indx_data = indx_data
        if Nmissing gt 0 then begin
          bg = nbim[indx_missing[0]]
        endif else begin
          bg = median(nbim)  
        endelse
        red_fitscube_addframe, fileassoc $
                               , red_rotation(temporary(nbim), ang $
                                              , background = bg $
                                              , full = ff, nthreads=nthreads) $
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
    wcs.hpln[0, 0, *, *] = hpln - double(image_scale) * (Nx-1)/2.d
    wcs.hpln[1, 0, *, *] = hpln + double(image_scale) * (Nx-1)/2.d
    wcs.hpln[0, 1, *, *] = hpln - double(image_scale) * (Nx-1)/2.d
    wcs.hpln[1, 1, *, *] = hpln + double(image_scale) * (Nx-1)/2.d
    
    wcs.hplt[0, 0, *, *] = hplt - double(image_scale) * (Ny-1)/2.d
    wcs.hplt[1, 0, *, *] = hplt - double(image_scale) * (Ny-1)/2.d
    wcs.hplt[0, 1, *, *] = hplt + double(image_scale) * (Ny-1)/2.d
    wcs.hplt[1, 1, *, *] = hplt + double(image_scale) * (Ny-1)/2.d
 
    ;; The wavelength is the tuning wavelength, corrected by the
    ;; prefilterfit.
    for ituning = 0, Ntuning-1 do $
       wcs[ituning].wave = utunwavelength[ituning]*1d9 - wave_shift[ituning]
    
    
    ;; Close fits file 
    self -> fitscube_finish, lun, wcs = wcs, direction = direction

    ;; Add cavity maps as WAVE distortions 
    if ~keyword_set(nocavitymap) then begin 
      if ~keyword_set(norotation) then begin
        red_fitscube_addcmap, filename $
                              , reform(red_rotation(cmap1,ang,full=ff, nthreads=nthreads),[dims[0:1],1,1,1])
      endif else begin
        red_fitscube_addcmap, filename, reform(cmap1,[dims[0:1],1,1,1])
      endelse
    endif 
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
    red_fitscube_addrespappl, filename, prefilter_curve, /tun

    tindx_r0 = where(time_r0 ge min(tavg_array) and time_r0 le max(tavg_array), Nt)
    if Nt gt 0 then begin
      self -> fitscube_addvarkeyword, filename, 'ATMOS_R0' $
         , metadata_r0[*, tindx_r0] $
         , comment = 'Atmospheric coherence length' $
         , tunit = 'm' $
         , extra_coordinate1 = [24, 8] $                                            ; WFS subfield sizes 
         , extra_labels      = ['WFSZ'] $                                           ; Axis labels for metadata_r0
         , extra_names       = ['WFS subfield size'] $                              ; Axis names for metadata_r0
         , extra_units       = ['pix'] $                                            ; Axis units for metadata_r0
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
         , tuning_selection = tuning_selection

    endif

    ;; Include the global WB image as an image extension
    ehdr=wbhdr
    fxaddpar, ehdr, 'XTENSION', 'IMAGE'
    sxdelpar, ehdr, 'SIMPLE'
    if keyword_set(norotation) then begin
      check_fits, wbim, ehdr, /update
    endif else begin
      check_fits, wbim_rot, ehdr, /update
    endelse
    fxaddpar, ehdr, 'DATE', red_timestamp(/utc, /iso)
    anchor = 'DATE'
    red_fitsaddkeyword, anchor = anchor, ehdr, 'EXTNAME', 'WBIMAGE', 'Wideband image'
    red_fitsaddkeyword, anchor = anchor, ehdr, 'PCOUNT', 0
    red_fitsaddkeyword, anchor = anchor, ehdr, 'GCOUNT', 1
    if keyword_set(norotation) then begin
      writefits, filename, wbim, ehdr, /append
    endif else begin
      writefits, filename, wbim_rot, ehdr, /append
    endelse

    ;; Correct intensity with respect to solar elevation and
    ;; exposure time.
    self -> fitscube_intensitycorr, filename $
       , intensitycorrmethod = intensitycorrmethod $
       , notimecheck = no_intensitycorr_timecheck

    if keyword_set(integer) then begin
      ;; Convert to integers
      self -> fitscube_integer, filename $
         , /delete $
         , outname = outname $
         , overwrite = overwrite
      filename = outname
    endif else begin
      if ~keyword_set(nomissing_nans) then begin
        ;; Set padding pixels to missing-data, i.e., NaN.
        self -> fitscube_missing, filename $
                                 , /noflip $
                                 , missing_type = 'nan' 
      endif
    endelse
    
    ;; Done with this scan.
    print, inam + ' : Narrowband scan cube stored in:'
    print, filename
    
  endfor                        ; iscan

  
end
