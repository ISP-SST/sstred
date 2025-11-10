; docformat = 'rst'

;+
; Wrapper method for make_nb_cube with or without demodulation.
; 
; Together with new methods make_nb_cube_stokes and
; make_nb_cube_intensity, based on earlier versions of make_nb_cube
; for CRISP and CHROMIS.
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
; 
; 
; 
; 
; 
; :Keywords:
; 
;    nomissing_nans : in, optional, type=boolean 
;
;      Do not set missing-data padding to NaN. (Set it to the median of
;      each frame instead.)
;
;    nopolarimetry : in, optional, type=boolean
;
;       For a polarimetric dataset, don't make a Stokes cube.
;       Instead combine all LC states for both cameras into a single
;       NB image per tuning, producing a cube similar to that for a
;       data set without polarimetry. (For a nonpolarlimetric dataset,
;       no effect.)
;
;    nostretch : in, optional, type=boolean
;   
;      Compute no intrascan stretch vectors if this is set.
;
;   
;   
; 
; 
; :History:
; 
;   2025-07-24 : MGL. First version as a common method for CRISP,
;                CRISP2, and CHROMIS, both Stokes cubes and cubes
;                without polarimetry. 
;
;   2025-11-02 : MGL. If wcfile undefined, let user select an existing
;                file.
; 
;-
pro red::make_nb_cube, wcfile $
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
                       , nopolarimetry = nopolarimetry $
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
    print, inam + ' : Keyword notimecor is deprecated. Use intensitycorrmethod="none" instead.'
    return
  endif

  ;; Default keywords
  if n_elements(cmap_fwhm) eq 0 then cmap_fwhm = 7.0 ;else fwhm = cmap_fwhm
  if n_elements(tiles) eq 0 or n_elements(clips) eq 0 then begin
    tiles = [8, 16, 32, 64, 84]
    clips = [8, 12,  4,  2,  1]
  endif
  if n_elements(nthreads) eq 0 then nthreads = 1 ; Default single thread

  if n_elements(wcfile) eq 0 then begin
    wcfiles = file_search('cubes_wb/*fits', count = Nfiles)
    if Nfiles eq 0 then return
    finfo = file_info(wcfiles)
    modtimes = finfo.mtime
    wcindx = sort(modtimes)
    wcfiles = wcfiles[wcindx]
    modtimes = modtimes[wcindx]
    timestamps = strarr(Nfiles)
    for i = 0, Nfiles-1 do begin
      get_date,dte,systime(elaps=modtimes[i]),/time
      timestamps[i] = dte
    endfor
    tmp = red_select_subset(timestamps + ' '+wcfiles, default = Nfiles-1 $
                            , indx = sindx $
                            , maxcount = 1 $
                            , qstring = 'Select a WB cube')
    wcfile = wcfiles[sindx[0]]
    print, wcfile
  endif
  
  ;; Make prpara
  red_make_prpara, prpara, clips         
  red_make_prpara, prpara, integer
  red_make_prpara, prpara, intensitycorrmethod
  red_make_prpara, prpara, cmap_fwhm
  red_make_prpara, prpara, nearest
  red_make_prpara, prpara, noaligncont 
  red_make_prpara, prpara, nocavitymap 
  red_make_prpara, prpara, nomissing_nans
  red_make_prpara, prpara, nostretch 
  red_make_prpara, prpara, np           
  red_make_prpara, prpara, overwrite
  red_make_prpara, prpara, tiles        
  red_make_prpara, prpara, wcfile

  
  ;; Camera/detector identification
  self -> cameras $
     , instrument = instrument $
     , Nnbcams = Nnbcams $
     , nb_detectors = nbdetectors $
     , nb_cameras = nbcameras $
     , wb_camera = wbcamera

  polarimetric_data = self -> polarimetric_data()

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

  ;; Default for wb cubes without direction parameter
  if n_elements(direction) eq 0 then direction = 0

  alignment_model = red_align_model(wbgfiles[0])
  
  ;; Don't do any stretching if wcgrid is all zeros.
  nostretch_temporal = total(abs(wcgrid)) eq 0 
  sclstr = 0
  if ~nostretch_temporal then sclstr = 1

  ;; Read wcs extension of wb file to get pointing info
  fxbopen, wlun, wcfile, 'WCS-TAB', wbdr
  ttype1 = fxpar(wbdr, 'TTYPE1')
  fxbread, wlun, wwcs, ttype1
  fxbclose, wlun

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
  
  for jj=0,n_elements(wbgfiles)-1 do begin
    search_dir = file_dirname(wbgfiles[jj])+'/'
    srch = '*_' + (strsplit(wbgfiles[jj],'._',/extract))[-3] +'_*'
    ff = file_search(search_dir + srch + extension) 
    red_append,files,ff
  endfor
  Nfiles = n_elements(files)

  if self.filetype eq 'MOMFBD' then mr = momfbd_read(wbgfiles[0],/nam)

  ;; FOV mask
  ;;  ;; If multiple directories, the fov_mask should be the same. Or we
;;  ;; have to think of something.
  spl = strsplit(wbgfiles[0],'/',/extract)
  cw = where(strmatch(spl,'*cfg*'))
  cfg_dir=strjoin(spl[0:cw],'/')
;;  if self.filetype eq 'MOMFBD' then begin
;;    mr = momfbd_read(wbgfiles[0],/nam)
;;  endif else begin              ; get cropping from cfg file    
;;    cfg_file = cfg_dir+'/'+'momfbd_reduc_'+wbgstates[0].prefilter+'_'+$
;;               string(wbgstates[0].scannumber,format='(I05)')+'.cfg'
;;    cfg = redux_readcfg(cfg_file)
;;    num_points = long(redux_cfggetkeyword(cfg, 'NUM_POINTS'))
;;    margin = num_points/8
;;    sim_xy = redux_cfggetkeyword(cfg, 'SIM_XY', count = cnt)
;;    if cnt gt 0 then begin
;;      sim_xy = rdx_str2ints(sim_xy)
;;      indx = indgen(n_elements(sim_xy)/2)*2
;;      indy = indx+1
;;      sim_x = sim_xy[indx]
;;      sim_y = sim_xy[indy]   
;;    endif else begin
;;      sim_x = rdx_str2ints(redux_cfggetkeyword(cfg, 'SIM_X'))
;;      sim_y = rdx_str2ints(redux_cfggetkeyword(cfg, 'SIM_Y'))
;;    endelse
;;    xx0 = min(sim_x) + margin - num_points/2 
;;    xx1 = max(sim_x) - margin + num_points/2 - 1
;;    yy0 = min(sim_y) + margin - num_points/2 
;;    yy1 = max(sim_y) - margin + num_points/2 - 1      
;;  endelse
  if file_test(cfg_dir+'/fov_mask.fits') then begin
    fov_mask = readfits(cfg_dir+'/fov_mask.fits')
    if self.filetype eq 'MOMFBD' then $
       fov_mask = red_crop_as_momfbd(fov_mask, mr) $
    else $
       fov_mask = fov_mask[xx0:xx1,yy0:yy1]
    fov_mask = red_rotate(fov_mask, direction)
  endif
  

  
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
  pertuningfiles  = files[complement]
  pertuningstates = states[complement]

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
  
  my_prefilters = pertuningstates[xindx].prefilter
  wav = utunwavelength

  ;; Unique nb prefilters
  unbprefs = red_uniquify(unbpertuningstates.prefilter, indx = unbprefindx)
  Nnbprefs = n_elements(unbprefs)
;  unbprefsref = dblarr(Nnbprefs)
;  
;  for inbpref = 0L, Nnbprefs-1 do begin
;    ;; This is the reference point of the fine tuning for this prefilter:
;    unbprefsref[inbpref] = double((strsplit(ufpi_states[0], '_', /extract))[0])
;  endfor                        ; inbpref
  
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
                               , cam = nbcameras $
                               , sel = nbindx, count = Nnb
  nbstates = pertuningstates[nbindx]
  nbfiles = pertuningfiles[nbindx]
  
  ;; Prepare for making output file names
  if(n_elements(odir) eq 0) then odir = self.out_dir + '/cubes_nb/' 
  
  ;; Do we want/can we make a fitscube with Stokes component images?
  case !true of
    keyword_set(nopolarimetry)                 : this_cube_is_stokes = 0
    n_elements(red_uniquify(nbstates.lc)) eq 1 : this_cube_is_stokes = 0
    polarimetric_data                          : this_cube_is_stokes = 1
    else                                       : this_cube_is_stokes = 0 
  endcase
  
  ofile = red_strreplace(file_basename(wcfile), 'wb', 'nb')
  if this_cube_is_stokes then begin
    ofile = red_strreplace(ofile, 'corrected', 'stokes_corrected')
    Nstokes = 4
  endif else begin
    Nstokes = 1
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
      return
    endelse
  endif

  file_mkdir, odir
  
  ;; Load prefilters
  self -> get_prefilterfit, unbpertuningstates $
     , units = units $
     , prefilter_curve = prefilter_curve $
     , prf = prf $
     , wave_shifts = wave_shifts

  ;; Do WB correction?
  wbcor = Nwb eq Nnb/Nnbcams and ~keyword_set(nostretch)

  ;; Spatial dimensions that match the WB cube
  Nx = wcND[0]
  Ny = wcND[1]
  
  ;; Make FITS header for the NB cube
  hdr = wchead                                      ; Start with the WB cube header
  red_headerinfo_deletestep, hdr, /all              ; Remove make_wb_cube steps 
  red_fitsdelkeyword, hdr, 'VAR_KEYS'               ; Variable keywords copied later
  red_fitsdelkeyword, hdr, 'STATE'                  ; Not a single state for cube 
  red_fitsdelkeyword, hdr, 'CHECKSUM'               ; Checksum for WB cube
  red_fitsdelkeyword, hdr, 'DATASUM'                ; Datasum for WB cube
  dindx = where(strmid(hdr, 0, 4) eq 'DATA', Ndata) ; DATA statistics keywords
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

  ;; WB and NB data come from different cameras.
  red_fitsaddkeyword, hdr, 'CAMERA',   nbstates[0].camera
  ;; Get DETGAIN, DETOFFS, DETMODEL, DETFIRM from .fitsheader file,
  ;; i.e., red_readhead(.momfbd file). This would have to be handled
  ;; differently with CRISP because each Stokes image is a mix of two
  ;; detectors. 
  mhd = red_readhead(nbfiles[0]) ; Header of momfbd output file
  red_fitscopykeyword, hdr, 'DETECTOR', mhd
  red_fitscopykeyword, hdr, 'DETGAIN',  mhd
  red_fitscopykeyword, hdr, 'DETOFFS',  mhd
  red_fitscopykeyword, hdr, 'DETMODEL', mhd
  red_fitscopykeyword, hdr, 'DETFIRM',  mhd

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


  ;; Observations metadata variables
  tbeg_array     = dblarr(Ntunings, Nscans) ; Time beginning for state
  tend_array     = dblarr(Ntunings, Nscans) ; Time end for state
  tavg_array     = dblarr(Ntunings, Nscans) ; Time average for state
  date_beg_array = strarr(Ntunings, Nscans) ; DATE-BEG for state
  date_end_array = strarr(Ntunings, Nscans) ; DATE-END for state
  date_avg_array = strarr(Ntunings, Nscans) ; DATE-AVG for state
  exp_array      = fltarr(Ntunings, Nscans) ; Total exposure time
  sexp_array     = fltarr(Ntunings, Nscans) ; Single exposure time
  nsum_array     = lonarr(Ntunings, Nscans) ; Number of summed exposures
  fnumsum_array  = strarr(Ntunings, Nscans) ; Raw frame numbers 



;  if prefilter eq '3950' && ~keyword_set(noaligncont) then $

  help, wbgfiles, utunwavelength, direction

  ;; For CHROMIS Ca II data, get misalignment between wavelengths.
  ;; (Get all zeroes for other types of data.)
  ashifts = self -> get_align_continuum(wbgfiles, utunwavelength, direction $
                                        , noaligncont = noaligncont)


  ;; Initialize fits file, set up for writing the data part.
  dims = [Nx, Ny, Ntunings, Nstokes, Nscans] 
  self -> fitscube_initialize, filename, hdr, lun, fileassoc, dims 
  if keyword_set(wbsave) then begin
    ;; Initialize also the WB version
    wbfilename = red_strreplace(filename, 'nb_', 'wbalign_')
    self -> fitscube_initialize, wbfilename, hdr, wblun, wbfileassoc, dims 
  endif

  ;; Write the actual data to the fitscube file. This is done
  ;; differently depending on whether we want a cube with Stokes data
  ;; or just intensity.
  case this_cube_is_stokes of
    
    !true : self -> make_nb_cube_stokes, wcfile $ 
       , ashifts = ashifts $
       , clips = clips $
       , cshift_mean = cshift_mean $
       , date_avg_array = date_avg_array  $
       , date_beg_array = date_beg_array  $ 
       , date_end_array = date_end_array  $
       , exp_array = exp_array $
       , fileassoc = fileassoc $
       , filename = filename $
       , fnumsum_array  = fnumsum_array $
       , fov_mask = fov_mask $
       , noremove_periodic = noremove_periodic $
       , nsum_array = nsum_array $  
       , nthreads = nthreads $
       , pertuningfiles = pertuningfiles $
       , pertuningstates = pertuningstates $
       , redemodulate = redemodulate $
       , remove_smallscale = remove_smallscale $
       , sexp_array = sexp_array $ 
       , tiles = tiles $
       , wbcor = wbcor $
       , wbfileassoc = wbfileassoc $
       , wbfilename = wbfilename $
       , wbsave = wbsave $
       , wcs = wcs
    
    !false : self -> make_nb_cube_intensity, wcfile $ 
       , ashifts = ashifts $
       , clips = clips $
       , cshift_mean = cshift_mean $
       , date_avg_array = date_avg_array  $
       , date_beg_array = date_beg_array  $ 
       , date_end_array = date_end_array  $
       , exp_array = exp_array $
       , fileassoc = fileassoc $
       , filename = filename $
       , fnumsum_array = fnumsum_array  $
       , fov_mask = fov_mask $
       , nsum_array = nsum_array    $  
       , nthreads = nthreads $
       , pertuningfiles = pertuningfiles $
       , pertuningstates = pertuningstates $
       , prefilter_curve = prefilter_curve $
       , remove_smallscale = remove_smallscale $
       , sexp_array = sexp_array     $ 
       , tiles = tiles $
       , wbcor = wbcor $
       , wbfileassoc = wbfileassoc $
       , wbfilename = wbfilename $
       , wbsave = wbsave $
       , wcs = wcs

    else : stop                 ; Cannot happen!
    
  endcase

  ;; Close fits file(s).
  free_lun, lun
  if keyword_set(wbsave) then free_lun, wblun 

  ;; Process some info from the cube making
  for iscan = 0L, Nscans-1 do begin
    for ituning = 0, Ntunings-1 do begin
      tbeg_array[ituning, iscan] = red_time2double((strsplit(date_beg_array[ituning, iscan],'T',/extract))[1])
      tend_array[ituning, iscan] = red_time2double((strsplit(date_end_array[ituning, iscan],'T',/extract))[1])
      tavg_array[ituning, iscan] = (tbeg_array[ituning, iscan] + tend_array[ituning, iscan]) / 2.
      date_avg_array[ituning, iscan] = self.isodate + 'T' + red_timestring(tavg_array[ituning, iscan])
    endfor                      ; ituning
  endfor                        ; iscan
  

  ;; Set WCS info
  wcs = replicate({  wave:dblarr(2,2) $
                     , hplt:dblarr(2,2) $
                     , hpln:dblarr(2,2) $
                     , time:dblarr(2,2) $
                  }, Ntunings, Nscans)

  ;; The narrowband cube is aligned to the wideband cube and all
  ;; narrowband scan positions are aligned to each other. So get hpln
  ;; and hplt from the wideband cube wcs coordinates, this should
  ;; apply to all frames in a scan.
  for iscan = 0L, Nscans-1 do begin
    for ituning = 0, Ntunings-1 do begin

      ;; We rely here on hpln and hplt being the first two tabulated
      ;; coordinates. To make this more general, we should get the
      ;; actual indices from the headers. Maybe later...
      wcs[ituning, iscan].hpln = reform(wwcs[0,*,*,iscan])
      wcs[ituning, iscan].hplt = reform(wwcs[1,*,*,iscan])
      
      ;; The wavelength is the tuning wavelength, corrected by the
      ;; prefilterfit.
      wcs[ituning, iscan].wave = utunwavelength[ituning]*1d9 - wave_shifts[ituning]
      
      ;; Time info as read from the data
      wcs[ituning, iscan].time = tavg_array[ituning, iscan]

    endfor                      ; ituning
  endfor                        ; iscan


  help, cshift_mean

  
  ;; Write the WCS info to the file(s)
  whdr = headfits(wcfile)
  csyer_spatial_value = fxpar(whdr, 'CSYER1', comment = csyer_spatial_comment)
  red_fitscube_addwcs, filename, wcs $
                       , csyer_spatial_value = csyer_spatial_value $
                       , csyer_spatial_comment = csyer_spatial_comment $
                       , dimensions = dims
  if keyword_set(wbsave) then begin 
    red_fitscube_addwcs, wbfilename, wcs $
                         , csyer_spatial_value = csyer_spatial_value $
                         , csyer_spatial_comment = csyer_spatial_comment $
                         , dimensions = dims
  endif
  
  ;; Copy the MOMFBD or bypass_momfbd step (or a mix):
  self -> headerinfo_copystep, filename, wcfile, stepnum = 1

  if 0 then begin
    ;; Need to copy prpara info from the stokes cube headers to hdr.
    ;;hh = red_readhead(snames[0, 0])
    ;;self -> headerinfo_copystep, filename, hh, prstep = 'DEMODULATION'
    self -> headerinfo_copystep, filename, snames[0, 0] , prstep = 'DEMODULATION'
  endif 

    
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

  
  ;; Now do the cavity maps
  
  cavitymaps = fltarr(Nx, Ny, 1, 1, Nscans)

  if ~keyword_set(nocavitymap) then begin

    ;; Read the original cavity map
    pindx = where(nbstates.prefilter ne '3999') ; No cavity map for the Ca II H continuum
    pindx = pindx[uniq(nbstates[pindx].prefilter, sort(nbstates[pindx].prefilter))]
    cprefs = nbstates[pindx].prefilter
    Ncprefs = n_elements(cprefs)

    ;; If multiple prefilters during the scan (as for CHROMIS Ca II),
    ;; calculate cmaps independently for them.
    for icprefs = 0, Ncprefs-1 do begin

      ;; Calculate average cmap if multiple cameras
      cmap1 = 0.0
      for icam = 0, Nnbcams-1 do begin

;        cfile = self.out_dir + 'flats/spectral_flats/' $
;                + strjoin([nbdetectors[icam] $
;                           , cprefs[icprefs] $
;                           , 'fit_results.sav'] $
;                          , '_')

        cfile = file_search(self.out_dir + 'flats/spectral_flats/' $
                            + nbdetectors[icam] + '_' $
                            + cprefs[icprefs] + '_fit_results.sav', count = Nsearch)

        if Nsearch ne 1 then stop else cfile = cfile[0]
        
        if ~file_test(cfile) then begin
          print, inam + ' : Error, calibration file not found -> '+cfile
          print, 'Please run fitgains for '+cprefs[icprefs]+' or continue without'
          print, 'cavity map for '+cprefs[icprefs]
          stop
        endif
        restore, cfile                 ; The cavity map is in a struct called "fit". 
        cmap = reform(fit.pars[1,*,*]) ; Unit is [Angstrom]
        cmap /= 10.                    ; Make it [nm]
        cmap = -cmap                   ; Change sign so lambda_correct = lambda + cmap
        fit = 0B                       ; Don't need the fit struct anymore.
        
        if keyword_set(remove_smallscale) then begin
          ;; If the small scale is already corrected, then include only the
          ;; low-resolution component in the metadata. The blurring kernel
          ;; should match how the low resolution component was removed when
          ;; making flats.
          npix = 30             ; Can we get this parameter from earlier headers?
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
        
      endfor                    ; icam
      
      if self.filetype eq 'MOMFBD' then begin
        ;; Crop the cavity map to the FOV of the momfbd-restored images.
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

      ;; Now make rotated copies of the cavity map
      for iscan = 0L, Nscans-1 do begin

        if ~keyword_set(nomissing_nans) then bg=!Values.F_NaN
        ;; Apply the same derot, align, dewarp as for the science data
        cmap11 = red_rotation(cmap1, ang[iscan], $
                              wcSHIFT[0,iscan], wcSHIFT[1,iscan], full=wcFF, $
                              stretch_grid=reform(wcGRID[iscan,*,*,*])*sclstr, $
                              nthreads=nthreads, nearest=nearest, $
                              unrotated_shifts = cshift_mean[*,icprefs,iscan], $
                              background=bg)          

        cavitymaps[0, 0, 0, 0, iscan] = cmap11

        ;; The following block of code is inactive but we want to keep
        ;; it around in case it is needed later. It does further
        ;; operations on the cavity maps based on what is done to the
        ;; science data, such as stretching based on the extra
        ;; per-tuning wideband objects, as well as blurring based on the
        ;; momfbd psfs. It should probably have a boolean keyword that
        ;; activates it. (This code will have to be updated to the
        ;; current pipeline style before it can be used.)
        if 0 then begin
          wb = (red_readdata(wbf[ss]))[x0:x1,y0:y1]
          for ww = 0L, nw - 1 do begin
            
            iwbf = strjoin([self.camwbtag, (strsplit(file_basename(st.ofiles[ww,ss]) $
                                                     , '.', /extract))[1:*]],'.')
            iwbf = file_dirname(st.ofiles[ww,ss]) + '/'+iwbf
              
            ;; Load images
            iwb = (red_readdata(iwbf))[x0:x1,y0:y1]
            im = momfbd_read(st.ofiles[ww,ss])
            
            ;; get dewarp from WB
            igrid = rdx_cdsgridnest(wb, iwb, itiles, iclip)
            
            ;; Convolve CMAP and apply wavelength dep. de-warp
            cmap2 = rdx_cstretch((red_mozaic(red_conv_cmap(cmap, im), /crop))[x0:x1, y0:y1], igrid, $
                                 nthreads=nthreads)
            
            ;; Derotate and shift
            cmap2 = red_rotation(temporary(cmap2), ang[ss], total(shift[0,ss]), $
                                 total(shift[1,ss]), full=wcFF, stretch_grid=reform(grid[ss,*,*,*]), $
                                 nearest=nearest, nthreads=nthreads,unrotated_shifts = cshift_mean[*,icprefs,iscan])
            
            ;; Time de-warp
            ;;cmap2 = red_stretch(temporary(cmap2), reform(grid[ss,*,*,*]))
            
          endfor                ; ww
        endif


      endfor                    ; iscan

      ;; For which tunings were this prefilter used?
;      tindx = where(scan_nbstates.prefilter eq cprefs[icprefs], Nt)
      tindx = where(unbpertuningstates.prefilter eq cprefs[icprefs], Nt)
      
      ;; Add cavity maps as WAVE distortions
      if Nt gt 0 then begin
        red_fitscube_addcmap, filename, cavitymaps $
                              , cmap_number = icprefs+1 $
                              , prefilter = cprefs[icprefs] $
                              , indx = tindx
      endif
      
    endfor                      ; icprefs
    
  endif

  red_message, 'Add some variable keywords.'

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

  undefine, fnumsum_total
  for iscan = 0L, Nscans-1 do for ituning = 0L, Ntunings-1 do $
     red_append, fnumsum_total, n_elements(rdx_str2ints(fnumsum_array[ituning,iscan]))
  fnumsum_total = rdx_ints2str([red_uniquify(fnumsum_total)])

  self -> fitscube_addvarkeyword, filename, 'FNUMSUM', fnumsum_array $
     , comment = 'Raw frame numbers' $
     , anchor = anchor $
     , keyword_value = fnumsum_total $
     , axis_numbers = [3, 5]

  ;; Correct intensity with respect to solar elevation and exposure
  ;; time.
  self -> fitscube_intensitycorr, filename $
     , intensitycorrmethod = intensitycorrmethod $
     , fitpref_time = fitpref_time 


  if this_cube_is_stokes and ~keyword_set(nocrosstalk) then begin

    print, 'Press any key to make crosstalk correction'
    q=get_kbrd()
    ;; Correct the cube for cross-talk, I --> Q,U,V.
    self -> fitscube_crosstalk, filename
    
  endif                         ;else


  
  if keyword_set(integer) then begin
    ;; Convert to integers
    self -> fitscube_integer, filename $
       , /delete $
       , outname = outname $
       , overwrite = overwrite
    filename = outname
  endif

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


  if ~keyword_set(noflipping) then $
     red_fitscube_flip, filename, flipfile = flipfile $
                        , overwrite = overwrite

  print, inam + ' : Narrowband cube stored in:'
  print, filename
  if ~keyword_set(noflipping) then print, flipfile
  
  if keyword_set(wbsave) then begin
    if ~keyword_set(noflipping) then red_fitscube_flip, wbfilename, flipfile = wbflipfile $
       , overwrite = overwrite
    red_message, 'Wideband align cube stored in: ' + wbfilename
    if ~keyword_set(noflipping) then print, wbflipfile
  endif


  
end
