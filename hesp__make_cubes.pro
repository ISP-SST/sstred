pro hesp::make_cubes, dir $
                      , direction = direction $
                      , nametag = nametag $
                      , odir = odir $
                      , overwrite = overwrite $
                      , rotation = rotation $
                      , subtract_meanang = subtract_meanang
  
  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  ;; Name of the instrument
  instrument = ((typename(self)).tolower())

  ;; Make prpara
  red_make_prpara, prpara, direction 
  red_make_prpara, prpara, rotation
  red_make_prpara, prpara, subtract_meanang  
  
  ;; Get metadata from logfiles
  red_logdata, self.isodate, time_r0, r0 = metadata_r0, ao_lock = ao_lock
  if keyword_set(use_turret_coordinates) then begin
    red_logdata, self.isodate, time_pointing, turret = metadata_pointing, rsun = rsun
  endif else begin
    red_logdata, self.isodate, time_pointing, diskpos = metadata_pointing, rsun = rsun
  endelse
  red_logdata, self.isodate, time_turret, azel = azel
  

  files = file_search(dir + 'img_fit.*.cube.fits', count = Nfiles)
  
  framenums = lonarr(Nfiles)
  for ifile = 0, Nfiles-1 do framenums[ifile] = long((strsplit(file_basename(files[ifile]), '.', /extract))[1])

  ;; Sort files
  indx = sort(framenums)
  framenums = framenums[indx]
  files = files[indx]
  
  prefilter = '10830'
  datestamp = self.isodate

  
  time = strarr(Nfiles)
  date = strarr(Nfiles)
  tmean = fltarr(Nfiles)

  hdr = red_readhead(files[0])  
  Naxis = fxpar(hdr, 'NAXIS*')
  Nwavelengths = Naxis[0]
  Nstokes = Naxis[1]
  Nx = Naxis[2]
  Ny = Naxis[3]
  
  x01y01 = [19, 72, 17, 82]
  
  x0 = x01y01[0] & x1 = x01y01[1] & y0 = x01y01[2] & y1 = x01y01[3]
  origNx = x1 - x0 + 1
  origNy = y1 - y0 + 1



  ;; Output files
  datatype = 'corrected'
  if n_elements(nametag) eq 0 then begin
    ;;  ofil = 'wb_'+midpart+'_'+datatype+'_im.fits'
    datatags = [datatype, 'im']
  endif else begin
    ;;   ofil = 'wb_'+midpart+'_'+nametag+'_'+datatype+'_im.fits'
    datatags = [nametag, datatype, 'im']
  endelse

  outdir = 'cubes/'
  file_wb = outdir + 'wb_test.fits'
  file_nb = outdir + red_strreplace(file_basename(file_wb), 'wb', 'nb')

  ;; Already done?
  if file_test(file_wb) then begin
    if keyword_set(overwrite) then begin
      print, inam + ' : Overwriting existing data cubes:'
      print, file_wb
      print, file_nb
    endif else begin
      print, inam + ' : These data cubes exist already:'
      print, file_wb
      print, file_nb
      return
    endelse
  endif else begin
    print, inam + ' : Will write these files:'
    print, file_wb
    print, file_nb
    file_mkdir, outdir
  endelse
  
  ;; Observations metadata variables
  tbeg_array     = dblarr(1, Nfiles) ; Time beginning for state
  tavg_array     = dblarr(1, Nfiles) ; Time average for state
  tend_array     = dblarr(1, Nfiles) ; Time end for state
  date_beg_array = strarr(1, Nfiles) ; DATE-BEG for state
  date_avg_array = strarr(1, Nfiles) ; DATE-AVG for state
  date_end_array = strarr(1, Nfiles) ; DATE-END for state
  exp_array      = fltarr(1, Nfiles) ; Total exposure time
  sexp_array     = fltarr(1, Nfiles) ; Single exposure time
  nsum_array     = lonarr(1, Nfiles) ; Number of summed exposures
  date_obs_array = strarr(Nfiles)    ; Datasets for each file

  dims_wb = [origNx, origNy, 1, 1, Nfiles]
  dims_nb = [origNx, origNy, Nwavelengths, Nstokes, Nfiles]
  
  cub_wb = fltarr(dims_wb)
  cub_nb = fltarr(dims_nb)

  ;; WCS info
  wcs_nb = replicate({ wave:dblarr(2,2) $
                       , hplt:dblarr(2,2) $
                       , hpln:dblarr(2,2) $
                       , time:dblarr(2,2) $
                     }, Nwavelengths, Nfiles)
  
  wcs_wb = replicate({ wave:dblarr(2,2) $
                       , hplt:dblarr(2,2) $
                       , hpln:dblarr(2,2) $
                       , time:dblarr(2,2) $
                     }, 1, Nfiles)
  
  wcs_wb[*, *].wave = 10830./10. ; [nm]

  ;; Read headers to get obs_time and load the images into cubes

  for ifile = 0L, Nfiles -1 do begin
    
    red_progressbar, ifile, Nfiles, 'Read headers and load the images into cubes'

    ;; The Stokes cube is the main extension
    tmp = readfits(files[ifile], /silent)
    tmp = tmp[*, *, x0:x1, y0:y1] >(-10) <10 ;                                <------- NOTE!!!
    cub_nb[*, *, *, *, ifile] = transpose(tmp,[2,3,0,1])
    
    ;; The context image is in extension 2
    tmp = readfits(files[ifile], ext = 2, /silent)
    cub_wb[*, *, 0, 0, ifile] = tmp[x0:x1, y0:y1]

    ;; The wavelength coordinates are in extension 1
    if ifile eq 0 then lambda = readfits(files[ifile], ext = 1, /silent) / 10.
    
    hasdateavg = 1
    date_avg = datestamp + 'T' + red_timestring(9*3600.+double(framenums[ifile])/100.)
    date_obs_array[ifile] = datestamp
    tavg_array[0, ifile] = red_time2double((strsplit(date_avg,'T',/extract))[1])

    if hasdateavg then begin
      date_avg_split = strsplit(date_avg, 'T', /extract, count = Nsplit)
      ddate = date_avg_split[0]
      if Nsplit gt 1 then ttime = date_avg_split[1] else undefine, ttime
    endif else undefine, ddate, ttime

    if n_elements(ddate) eq 0 then begin
      print, inam+' : No date and time information for file '+strtrim(ufiles[ifile], 2)
      stop
    endif else begin
      date[ifile] = ddate
      time[ifile] = ttime
    endelse

    ;; Wavelength and time
    for iwave = 0, Nwavelengths-1 do begin
      wcs_nb[iwave, ifile].wave = lambda[iwave]
    endfor                      ; iwave
    
  endfor                        ; ifile



  ;; Intensity correction
  
  ;; Rotation?
  

  timestamps = red_timestring(tavg_array[0])
  filenos_actual = red_collapserange(framenums,r='',l='')
  ;; POINT_ID
  date_obs = datestamp          ;fxpar(hdr, 'DATE-OBS', count = Ndate_obs)  
  if n_elements(point_id) eq 0 then begin
    point_id = fxpar(hdr, 'POINT-ID', count = Npoint_id)
    if Npoint_id eq 0 then point_id = date_obs
  endif


  ;; Make time tabhdu extension with Nfiles rows
  s_array = lonarr(Nfiles)
;  s_array[0] = wstates.filenumber
  t_array = dblarr(1, Nfiles)
  t_array[0] = red_time2double(time) ; In [s] since midnight

  ;; Get pointing at center of FOV
  red_wcs_hpl_coords, t_array[0, *], metadata_pointing, time_pointing $
                      , hpln, hplt
  
  ;; Let's smooth coordinates.
  dt = (t_array[0,-1] - t_array[0,0]) / 60. ; minutes
  if dt le 15. or Nfiles le 3 then fit_expr = 'P[0] + X*P[1]'
  if dt gt 15. and Nfiles gt 3 then fit_expr = 'P[0] + X*P[1] + X*X*P[2]'
  pp = mpfitexpr(fit_expr, t_array[0,*], hpln)
  hpln = red_evalexpr(fit_expr, t_array[0,*], pp)
  pp = mpfitexpr(fit_expr, t_array[0,*], hplt)
  hplt = red_evalexpr(fit_expr, t_array[0,*], pp)
  
  ;; But what we want to tabulate is the pointing in the corners of
  ;; the FOV. Assume hpln and hplt are the coordinates of the center
  ;; of the FOV.
  wcs_wb.hpln[0, 0] = hpln - double(self.image_scale) * (Nx-1)/2.d
  wcs_wb.hpln[1, 0] = hpln + double(self.image_scale) * (Nx-1)/2.d
  wcs_wb.hpln[0, 1] = hpln - double(self.image_scale) * (Nx-1)/2.d
  wcs_wb.hpln[1, 1] = hpln + double(self.image_scale) * (Nx-1)/2.d
  
  wcs_wb.hplt[0, 0] = hplt - double(self.image_scale) * (Ny-1)/2.d
  wcs_wb.hplt[1, 0] = hplt - double(self.image_scale) * (Ny-1)/2.d
  wcs_wb.hplt[0, 1] = hplt + double(self.image_scale) * (Ny-1)/2.d
  wcs_wb.hplt[1, 1] = hplt + double(self.image_scale) * (Ny-1)/2.d
  
  ;; Rebin hpln and hplt to match the nb wcs dimensions
  hpln = rebin(hpln,2,2,Nwavelengths,Nfiles,/samp) 
  hplt = rebin(hplt,2,2,Nwavelengths,Nfiles,/samp) 
  
  for iwave = 0, Nwavelengths-1 do begin
    wcs_nb.hpln[0, 0] = reform(hpln[0, 0, *, *]) - double(self.image_scale) * (Nx-1)/2.d
    wcs_nb.hpln[1, 0] = reform(hpln[1, 0, *, *]) + double(self.image_scale) * (Nx-1)/2.d
    wcs_nb.hpln[0, 1] = reform(hpln[0, 1, *, *]) - double(self.image_scale) * (Nx-1)/2.d
    wcs_nb.hpln[1, 1] = reform(hpln[1, 1, *, *]) + double(self.image_scale) * (Nx-1)/2.d
    
    wcs_nb.hplt[0, 0] = reform(hplt[0, 0, *, *]) - double(self.image_scale) * (Ny-1)/2.d
    wcs_nb.hplt[1, 0] = reform(hplt[1, 0, *, *]) - double(self.image_scale) * (Ny-1)/2.d
    wcs_nb.hplt[0, 1] = reform(hplt[0, 1, *, *]) + double(self.image_scale) * (Ny-1)/2.d
    wcs_nb.hplt[1, 1] = reform(hplt[1, 1, *, *]) + double(self.image_scale) * (Ny-1)/2.d
  endfor                        ; iwave
  
    
  ;; Base cube headers on first input file header

  ;; Delete some keywords that do not make sense for WB cubes.
  red_fitsdelkeyword, hdr, 'FNUMSUM'    ; May add this later
  red_fitsdelkeyword, hdr, 'FRAMENUM'   ; Not applicable here
  red_fitsdelkeyword, hdr, 'STATE'      ; State info is in WCS and FILTER1 keywords

  anchor = 'DATE'

  dateref = self.isodate+'T00:00:00.000000' ; Midnight
  red_fitsaddkeyword, hdr, 'DATEREF', dateref, 'Reference time in ISO-8601', after = 'DATE'

  ;; Add some keywords
  red_fitsaddkeyword, anchor = anchor, hdr, 'OBS_HDU', 1
  red_fitsaddkeyword, anchor = anchor, hdr, 'INSTRUME', 'HeSP'

  red_fitsaddkeyword, anchor = anchor, hdr, 'POINT_ID', point_id
;  red_fitsaddkeyword, anchor = anchor, hdr, 'INFILES', strjoin(files,','), 'Concatenated files'
  red_fitsaddkeyword, anchor = anchor, hdr, 'BUNIT', 'dn', 'Units in array: digital number'
  ;red_fitsaddkeyword, anchor = anchor, hdr, 'BUNIT', 'W m^-2 Hz^-1 sr^-1', 'Units in array'
  red_fitsaddkeyword, anchor = anchor, hdr, 'BTYPE', 'Intensity', 'Type of data in array'

  red_fitsaddkeyword, anchor = anchor, hdr, 'TIMESYS', 'UTC'
  red_fitsaddkeyword, anchor = anchor, hdr, 'WAVEBAND', 'He 10830', ''
  red_fitsaddkeyword, anchor = anchor, hdr, 'WAVEUNIT', -9, 'WAVELNTH in units 10^WAVEUNIT m = nm'
  red_fitsaddkeyword, anchor = anchor, hdr, 'DATE-OBS', date_obs, 'Data set timestamp'
  red_fitsaddkeyword, anchor = anchor, hdr, 'SOLARNET', 0.5, 'Fully SOLARNET-compliant=1.0, partially=0.5'
;  red_fitsaddkeyword, anchor = anchor, hdr, '', '', ''

;  ;; Cadence info.
;  if Nfiles gt 1 then begin
;    dt = float((red_differential(dts))[1:*])
;    red_fitsaddkeyword, anchor = anchor, hdr, 'CADAVG', mean(dt), 'Average of actual cadence'
;    ;; If we can access the planned cadence, it should go into keyword CADENCE.
;    if Nfiles gt 2 then begin
;      red_fitsaddkeyword, anchor = anchor, hdr, 'CADMIN', min(dt),      'Minimum of actual cadence'
;      red_fitsaddkeyword, anchor = anchor, hdr, 'CADMAX', max(dt),      'Maximum of actual cadence'
;      red_fitsaddkeyword, anchor = anchor, hdr, 'CADVAR', stddev(dt)^2, 'Variance of actual cadence'
;    endif
;  endif


  
  hdr_wb = hdr
  hdr_nb = hdr

  check_fits, cub_wb, hdr_wb, /update   ; Get dimensions right
  check_fits, cub_nb, hdr_nb, /update   ; Get dimensions right

  ;; Initialize the output files
  self -> fitscube_initialize, file_wb, hdr_wb, lun_wb, fileassoc_wb, dims_wb
  self -> fitscube_initialize, file_nb, hdr_nb, lun_nb, fileassoc_nb, dims_nb

  for ifile = 0, Nfiles-1 do begin

    ;; Add the wideband frame
    red_fitscube_addframe, fileassoc_wb, cub_wb[*, *, 0, 0, ifile] $
                           , iframe = ifile
    
    for iwave = 0, Nwavelengths-1 do begin
      for istokes = 0, Nstokes-1 do begin
        red_fitscube_addframe, fileassoc_nb, cub_nb[*, *, iwave, istokes, ifile] $
                               , iscan = ifile, ituning = iwave, istokes = istokes
      endfor                    ; istokes
    endfor                      ; iwave

    wcs_wb[*, ifile].time = t_array[ifile]
    wcs_nb[*, ifile].time = t_array[ifile]

  endfor                        ; ifile

  self -> fitscube_finish, lun_nb, wcs = wcs_nb, direction = direction
  self -> fitscube_finish, lun_wb, wcs = wcs_wb, direction = direction



  ;; Add info about this step
  self -> headerinfo_addstep, hdr_wb $
                              , anchor = 'SOLARNET' $
                              , files = files $
                              , cubefile = file_wb $
                              , prstep = 'CONCATENATION' $
                              , prpara = prpara $
                              , prproc = inam
  self -> headerinfo_addstep, hdr_nb $
                              , anchor = 'SOLARNET' $
                              , files = files $
                              , cubefile = file_nb $
                              , prstep = 'CONCATENATION' $
                              , prpara = prpara $
                              , prproc = inam

;  self -> fitscube_missing, file_wb, /noflip, missing_type = 'nan'
;  self -> fitscube_missing, file_nb, /noflip, missing_type = 'nan'


  ;; Add WCS

  ;; Add variable keywords
  

  print, inam + ' : Wrote these files:'
  print, file_wb
  print, file_nb
  
end




a = hespred(/dev)

dir = '/scratch/jleen/tmp/mats/'
a -> make_cubes, dir, /over



end
