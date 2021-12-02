; docformat = 'rst'

;+
; Concatenate multiple fitscubes (nb or scan) into a single fitscube.
;
; Based on assumption that the input fitscubes have the same tuning
; points and polarimetry status. They should also share the same
; (approximate) pointing. They are separated in time and are all
; derotated to Solar N up. It is the user's responsibility to
; make sure this is the case.
;
; Do not use already concatenated files or multi-dataset nb cubes as
; input! 
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
;   infiles : in, type=strarr
; 
;      Array of file names/paths of fitscube files to concatenate into
;      a single file.
; 
;   outfile : in, optional, type=string, default="based on temporally first file"
; 
;      File name/path for the concatenated fitscube.
; 
; :Keywords:
;
;   noflipping : in, optional, type=boolean
;
;      Do not make a spectral cube version.
;
;   overwrite : in, optional, type=boolean
;
;      Overwrite existing outfile.  
;
;   point_id : in, optional, type=string, default="From first file"
;
;      Value for the POINT_ID header keyword. By default the value if
;      POINT_ID in the temporally first file or, if that does not
;      exist, the value of DATE-OBS in the temporally first file.
; 
; :History:
; 
;   2021-11-18 : MGL. First version.
; 
;-
pro red::fitscube_concatenate, infiles, outfile $
                               , noflipping = noflipping $
                               , overwrite = overwrite $
                               , point_id = point_id
  

  ;; Name of this method
  inam = red_subprogram(/low, calling = inam1)

  ;; Make prpara
  red_make_prpara, prpara, point_id 
  
  Nfiles = n_elements(infiles)
  
  date_begs = strarr(Nfiles)
  naxis = lonarr(5, Nfiles)

  for ifile = 0, Nfiles-1 do begin
    if ~file_test(infiles[ifile]) then stop
    h = headfits(infiles[ifile])
    date_begs[ifile] = fxpar(h, 'DATE-BEG')
    naxis[*, ifile] = fxpar(h, 'NAXIS*')
  endfor                        ; ifile
  
  h = headfits(infiles[0])

  date_obs = fxpar(h, 'DATE-OBS')
  isodate = (strsplit(date_obs, 'T', /extract))[0]

  if n_elements(point_id) eq 0 then begin
    point_id = fxpar(h, 'POINT-ID', count = Npoint_id)
    if Npoint_id eq 0 then point_id = date_obs
  endif
  
  ;; Check dimensions
  if max(red_differential(reform(naxis[2, *]))) ne 0 then stop ; Number of tuning points varies!  
  if max(red_differential(reform(naxis[3, *]))) ne 0 then stop   ; Don't mix polarimetric and non-polarimetric data.

  ;; Get dimensions
  Nx = long(max(naxis[0, *]))         ; Make room for largest frames
  Ny = long(max(naxis[1, *]))         ; Make room for largest frames
  Ntuning = naxis[2, 0]               ; 
  Nstokes = naxis[3, 0]               ; 
  Nscans = long(total(naxis[4, *]))   ; Total number of scans

  dimensions = [Nx, Ny, Ntuning, Nstokes, Nscans]
 
  ;; Get metadata from logfiles
  red_logdata, isodate, time_r0, r0 = metadata_r0, ao_lock = ao_lock
  red_logdata, isodate, time_pointing, diskpos = metadata_pointing, rsun = rsun
  red_logdata, isodate, time_turret, azel = azel

  ;; Observations metadata varaibles
  date_beg_array = strarr(Ntuning, Nscans)   ; DATE-BEG for state
  date_end_array = strarr(Ntuning, Nscans)   ; DATE-END for state
  date_avg_array = strarr(Ntuning, Nscans)   ; DATE-AVG for state
  exp_array      = fltarr(Ntuning, Nscans)   ; Total exposure time XPOSURE
  sexp_array     = fltarr(Ntuning, Nscans)   ; Single exposure time TEXPOSUR
  nsum_array     = lonarr(Ntuning, Nscans)   ; Number of summed exposures NSUMEXP
  scannum_array  = lonarr(Nscans)            ; Scan number from which state is copied  
  date_obs_array = strarr(Nscans)            ; Dataset from which state is copied  

  wcs = replicate({ wave:dblarr(2,2) $
                    , hplt:dblarr(2,2) $ 
                    , hpln:dblarr(2,2) $
                    , time:dblarr(2,2) $
                  }, Ntuning, Nscans)
   
  ;; Sort files based on DATE-BEG. This way scans will be concatenated
  ;; in the correct order, assuming the cubes do not overlap in the
  ;; time dimension.
  indx = sort(red_time2double(strmid(date_begs,11)))
  infiles = infiles[indx]
  date_begs = date_begs[indx]  

  self -> extractstates, infiles, instates

  iscan = 0
  for ifile = 0, Nfiles-1 do begin

    Nscans_in = naxis[4, ifile] ; Number of scans in this file

    tmp = red_fitsgetkeyword(infiles[ifile], 'DATE-BEG', variable_values = date_beg_values)  
    tmp = red_fitsgetkeyword(infiles[ifile], 'DATE-END', variable_values = date_end_values)  
    tmp = red_fitsgetkeyword(infiles[ifile], 'DATE-AVG', variable_values = date_avg_values)      
    tmp = red_fitsgetkeyword(infiles[ifile], 'XPOSURE',  variable_values = xposure_values)      
    tmp = red_fitsgetkeyword(infiles[ifile], 'TEXPOSUR', variable_values = texposur_values)      
    tmp = red_fitsgetkeyword(infiles[ifile], 'NSUMEXP',  variable_values = nsumexp_values)

    ;; SCANNUM is not a variable keywords in scancubes
    tmp = red_fitsgetkeyword(infiles[ifile], 'SCANNUM',  variable_values = scannum_values)      
    if n_tags(scannum_values) lt 3 then scannum_values = { values:tmp }
    
    ;; DATE-OBS is not always a variable keyword.
    tmp = red_fitsgetkeyword(infiles[ifile], 'DATE-OBS', variable_values = date_obs_values)    
    if n_tags(date_obs_values) lt 3 then date_obs_values = { values:replicate(tmp, Nscans_in) }
    
    red_fitscube_getwcs, infiles[ifile]  $
                         , coordinates = coordinates $
                         , distortions = distortions

    if naxis[0, ifile] ne Nx and  naxis[1, ifile] ne Ny then begin
      ;; Note that hpln and hplt are corner coordinates so need to be
      ;; scaled to the frame size used here. (Better to use mid point and
      ;; calculate corners from that?) 
      facx = float(Nx)/float(naxis[0, ifile])      
      facy = float(Ny)/float(naxis[1, ifile])
      coordinates.hpln *= facx
      coordinates.hplt *= facy
    endif
    
    red_append, scannos, red_collapserange(scannum_values.values, /nobrackets)
    
    for i = 0L, Nscans_in-1 do begin
      date_beg_array[*, i+iscan] = (reform(date_beg_values.values))[*, i]      
      date_end_array[*, i+iscan] = (reform(date_end_values.values))[*, i]      
      date_avg_array[*, i+iscan] = (reform(date_avg_values.values))[*, i]            
      exp_array[*, i+iscan]      = (reform(xposure_values.values))[*, i]   
      sexp_array[*, i+iscan]     = (reform(texposur_values.values))[*, i] 
      nsum_array[*, i+iscan]     = (reform(nsumexp_values.values))[*, i] 
      scannum_array[i+iscan]     = (reform(scannum_values.values))[i]
      date_obs_array[i+iscan]    = (reform(date_obs_values.values))[i]

      wcs[*, i+iscan] = coordinates[*, i]

      ;; Cavity maps
      if ifile eq 0 and i eq 0 then begin
        Ncprefs = n_elements(distortions) ; # of prefilters for which there are cmaps
        wcs_distortions = fltarr(Nx, Ny, 1, 1, Nscans, Ncprefs)
        cprefs = strarr(Ncprefs)
        wcs_tun_indx = strarr(Ncprefs)
      endif
      for icpref = 0, Ncprefs-1 do begin
        wcs_distortions[*, *, *, *, i+iscan, icpref] = red_centerpic(distortions[icpref].wave[*, *, *, *, i] $
                                                                     , xsize = Nx, ysize = Ny $
                                                                     , z = !Values.F_NaN)
      endfor                    ; icpref

    endfor                      ; i

    for icpref = 0, Ncprefs-1 do begin
      wcs_tun_indx[icpref] = distortions[icpref].tun_index        
      cprefs[icpref] = ''       ; Should be set to the prefilter tag. Not crucial but nicer.         
    endfor                      ; icpref
    
    iscan += Nscans_in
    
  endfor                        ; ifile
  
  ;; Make default outfile name if needed.
  if n_elements(outfile) eq 0 then begin
    timestamps = (stregex(instates[indx].filename,'T([0-2][0-9]:[0-5][0-9]:[0-6][0-9])_',/extra,/sub))[1,*] 
    
    outfile = 'cubes_concatenated/' $
              + red_fitscube_filename('nb' $
                                      , instates[0].prefilter $
                                      , timestamps $                                    
                                      , scannos $                                    
                                      , point_id $
                                      , datatags = ['im'] $
                                     )
  endif
  
  odir = file_dirname(outfile)
  file_mkdir, odir

  print, outfile
  
  ;; Create a header based on the first file's header
  hdr = headfits(infiles[0])
  anchor = 'DATE-OBS'

  ;; Add POINT_ID if not already there, base it on DATE-OBS. Add other
  ;; keywords.
  red_fitsaddkeyword, anchor = anchor, hdr, 'POINT_ID', point_id, 'Identify pointing'
  red_fitsaddkeyword, anchor = anchor, hdr, 'INFILES', strjoin(infiles,','), 'Concatenated files'
  
  ;; Add info about this step
  self -> headerinfo_addstep, hdr $
                              , prstep = 'CONCATENATION' $
                              , prpara = prpara $
                              , prproc = inam $
                              , comment_prref = infiles
  
  
  ;; Initialize the outfile.
  self -> fitscube_initialize, outfile, hdr, lun, fileassoc, dimensions, wcs = wcs
  
  ;; Loop over infiles, copying frames to the outfile. Apply shifts
  ;; and stretch vectors if calculated. 
  iscan = 0L                    ; iscan for output file
  iprogress = 0L
  Nprogress = long(Nscans) * long(Ntuning) * long(Nstokes)
  
  for ifile = 0, Nfiles-1 do begin
    
    Nscans_in = naxis[4, ifile]
    
    for iscan_in = 0L, Nscans_in-1 do begin
      
      for ituning = 0L, Ntuning-1 do begin      
        for istokes = 0L, Nstokes-1 do begin
          
          red_progressbar, iprogress, Nprogress $
                           , infiles[ifile]+' '+strjoin(strtrim([ifile, iscan_in, ituning, istokes], 2), ',')
          
          red_fitscube_getframe, infiles[ifile], frame $
                                 , iscan = iscan_in, ituning = ituning, istokes = istokes

          ;;  Set missing pixels to NaN.
          red_missing, frame, missing_type_wanted = 'nan', /inplace
          
          ;; Pad to largest frame size
          frame = red_centerpic(frame $
                                , xsize = Nx, ysize = Ny $
                                , z = !Values.F_NaN)
          
          red_fitscube_addframe, fileassoc, frame $
                                 , iscan = iscan, ituning = ituning, istokes = istokes

          iprogress++

        endfor                  ; istokes
      endfor                    ; ituning

      iscan++                   ; Advance iscan for output file     

    endfor                      ; iscan_in
    
  endfor                        ; ifile

  ;; Finish the outfile 
  self -> fitscube_finish, lun, wcs=wcs

  ;; Add cavity maps
  for icpref = 0, Ncprefs-1 do begin
    red_fitscube_addcmap, outfile, wcs_distortions[*, *, *, *, *, icpref] $
                          , cmap_number = icpref+1 $
                          , prefilter = cprefs[icpref] $
                          , indx = red_expandrange(wcs_tun_indx[icpref])
  endfor                        ; icpref

  
  
  ;; Add variable keywords.
  self -> fitscube_addvarkeyword, outfile, 'DATE-BEG', date_beg_array $
                                  , anchor = anchor $
                                  , comment = 'Beginning time of observation' $
                                  , keyword_method = 'first' $
                                  , axis_numbers = [3, 5]
  
  self -> fitscube_addvarkeyword, outfile, 'DATE-END', date_end_array $
                                  , anchor = anchor $
                                  , comment = 'End time of observation' $
                                  , keyword_method = 'last' $
                                  , axis_numbers = [3, 5]
  
  self -> fitscube_addvarkeyword, outfile, 'DATE-AVG', date_avg_array $
                                  , anchor = anchor $
                                  , comment = 'Average time of observation' $
                                  , keyword_value = self.isodate + 'T' + red_timestring(mean(wcs.time)) $
                                  , axis_numbers = [3, 5]
  
  self -> fitscube_addvarkeyword, outfile, 'SCANNUM', scannum_array $
                                  , comment = 'Scan number' $
                                  , anchor = anchor $
                                  , keyword_method = 'first' $
                                  , axis_numbers = 5
  
  self -> fitscube_addvarkeyword, outfile, 'DATE-OBS', date_obs_array $
                                  , comment = 'Dataset' $
                                  , anchor = anchor $
                                  , keyword_method = 'first' $
                                  , axis_numbers = 5
  
  self -> fitscube_addvarkeyword, outfile, 'XPOSURE', exp_array $
                                  , comment = 'Summed exposure time' $
                                  , anchor = anchor $
                                  , tunit = 's' $
                                  , keyword_method = 'median' $
                                  , axis_numbers = [3, 5] 

  self -> fitscube_addvarkeyword, outfile, 'TEXPOSUR', sexp_array $
                                  , comment = '[s] Single-exposure time' $
                                  , anchor = anchor $
                                  , tunit = 's' $
                                  , keyword_method = 'median' $
                                  , axis_numbers = [3, 5] 

  self -> fitscube_addvarkeyword, outfile, 'NSUMEXP', nsum_array $
                                  , comment = 'Number of summed exposures' $
                                  , anchor = anchor $
                                  , keyword_method = 'median' $
                                  , axis_numbers = [3, 5]

  
  tindx_r0 = where(time_r0 ge min(wcs.time) and time_r0 le max(wcs.time), Nt)
  if Nt gt 0 then begin
    r0dims = size(metadata_r0, /dim) ; We have not always had two r0 values.
    if r0dims[0] eq 1 then extra_coordinate1 = [24] else extra_coordinate1 = [24, 8]
    self -> fitscube_addvarkeyword, outfile, 'ATMOS_R0' $
                                    , metadata_r0[*, tindx_r0] $
                                    , anchor = anchor $
                                    , comment = 'Atmospheric coherence length' $
                                    , tunit = 'm' $
                                    , extra_coordinate1 = extra_coordinate1 $     ; WFS subfield sizes 
                                    , extra_labels      = ['WFSZ'] $              ; Axis labels for metadata_r0
                                    , extra_names       = ['WFS subfield size'] $ ; Axis names for metadata_r0
                                    , extra_units       = ['pix'] $               ; Axis units for metadata_r0
                                    , keyword_method = 'mean' $
                                    , time_coordinate = time_r0[tindx_r0] $
                                    , time_unit       = 's'

    self -> fitscube_addvarkeyword, outfile, 'AO_LOCK' $
                                    , ao_lock[tindx_r0] $
                                    , anchor = anchor $
                                    , comment = 'Fraction of time the AO was locking, 2s average' $
                                    , keyword_method = 'mean' $
                                    , time_coordinate = time_r0[tindx_r0] $
                                    , time_unit       = 's'

  endif

  tindx_turret = where(time_turret ge min(wcs.time) and time_turret le max(wcs.time), Nt)
  if Nt gt 0 then begin

    self -> fitscube_addvarkeyword, outfile, 'ELEV_ANG' $
                                    , reform(azel[1, tindx_turret]) $
                                    , anchor = anchor $
                                    , comment = 'Elevation angle' $
                                    , keyword_method = 'mean' $
                                    , time_coordinate = time_turret[tindx_turret] $
                                    , time_unit       = 'deg'
    
  end

  ;; Store infiles strarr in outfile. One per frame so we have both
  ;; strarr and scanno for each frame.

  ;; Optionally make spectral cube.
  if ~keyword_set(noflipping) then $
     red_fitscube_flip, outfile, flipfile = flipfile $
                        , overwrite = overwrite
  
  print, inam + ' : Concatenated cube stored in:'
  print, outfile
  if ~keyword_set(noflipping) then print, flipfile
  

end


;cd, '/scratch/olexa/2020-04-21/CHROMIS'
;a = chromisred(/dev)
;infiles = 'cubes_nb/nb_3950_2020-04-21T07:57:08_scans='+['0-10', '11-20']+'_corrected_im.fits'
;
;a -> fitscube_concatenate, infiles, outfile, /overwrite
;
;red_fitscube_getwcs, outfile  $
;                     , coordinates = coordinates $
;                     , distortions = distortions
;
;help, coordinates, distortions
;
;stop

cd, '/scratch/mats/2016.09.19/CRISP-aftersummer'
a = crispred(/dev)

case 2 of
  0 : infiles = 'cubes_scan/nb_6302_2016-09-19T09:30:20_scan='+['0', '4', '2']+'_stokes_corrected.fits'
  1 : infiles = 'cubes_scan/nb_6302_2016-09-19T'+['11:05:36', '11:08:53']+'_scan=0_corrected.fits'
  2 : infiles = 'cubes_nb/' + ['nb_6302_2016-09-19T09:28:36_scans=0,1_stokes_corrected_im.fits', $
                               'nb_6302_2016-09-19T09:30:20_scans=0-2_stokes_corrected_im.fits']
endcase

a -> fitscube_concatenate, infiles, outfile, /noflipping, /overwrite




end
