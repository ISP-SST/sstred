; docformat = 'rst'

;+
; 
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
; :Params:
;
;   work_dir : in, type = string
;
;     A name of the directory to set up the directory tree at and to
;     where the config file and the doit.pro script should be put.
;
;   root_dir : in, type = string
;
;     A name of the directory where the raw data is stored for a date
;     given in isodate.
;
;   cfgfile : in, type = string
;
;     A name of the config file, the default value is 'config.txt'.
;
;   scriptfile : in, type = string
;
;     A name of the set-up script, the default value is 'doit.pro'.
;
;   isodate : in, type = string
;
;     A date given in ISO format (e.g., 2017-08-23) to search raw data
;     for in root_dir.
; 
; :Keywords:
;
;    ampm_cutoff : in, optional, type=string, default='13:00:00'
;
;      Data collected before this time belongs to AM observations,
;      everything else to PM.
;
;    calibrations_only : in, optional, type=boolean
;
;      Set up to process calibration data only.
; 
;    lapalma_setup : in, optional, type=boolean
; 
;      Set up the workdirs to sum all calibration data into
;      timestamped directories, like flats/hh:mm:ss plus make
;      quicklooks.
; 
;    no_observer_metadata : in, optional, type=boolean
;
;      Do not check for OBSERVER metadata in raw data or offer to
;      add it by hand. 
;
;    old_dir : in, optional, type = string
;
;      Copy files from this directory, in particular summed
;      calibration data. Useful if you summed the calibration data in
;      La Palma.
; 
; :History:
; 
;    2017-04-07 : MGL. Revised selection of commands.
; 
;    2017-04-18 : MGL. Improve script file.
; 
;    2017-04-20 : MGL. Do not include left-behind calibrations data
;                 directories in the science data dirs.
; 
;    2017-05-12 : MGL. Use extraclip rather than margin when calling
;                 prepmomfbd, and don't use an empty pref keyword. 
;
;    2017-07-05 : MGL. New keyword calibrations_only.
;
;    2017-08-07 : MGL. Set nthreads in script file and use when
;                 calling the summing methods.
;
;    2017-08-10 : MGL. When in /calibrations_only mode, write summed
;                 output to timestamp subdirectories, with softlinks
;                 to ordinary outdir for darks and flats.
;
;    2017-08-16 : MGL. Stop early if there are no darks and/or flats.
; 
;    2017-07-19 : THI. Change pinholecalib parameter nref defaut value
;                 to 10.
;
;    2018-06-11 : MGL. Do hrz conversion here rather than putting the
;                 commands into the doit script. Split several
;                 commands into one per prefilter.
;
;    2019-07-01 : MGL. New keyword old_dir.
; 
;    2020-03-11 : MGL. Set direction in config file.
; 
;    2020-04-06 : MGL. Set rotation in config file.
; 
;    2020-06-11 : MGL. Add fit_wb_diskcenter step.
; 
;    2021-01-25 : MGL. Check for OBSERVER metadata keyword in raw
;                 data. New keyword no_observer_metadata.
;
;   2025-03-27 : MGL. New keyword lapalma_setup. 
;
;   2025-04-03: MGL. New keyword ampm_cutoff. 
;
;   2025-05-16: MGL. Find any cameras, Chromis-[WDNTR].
; 
;-
pro red_setupworkdir_chromis, work_dir, root_dir, cfgfile, scriptfile, isodate $
                              , ampm_cutoff = ampm_cutoff $
                              , calibrations_only = calibrations_only $
                              , lapalma_setup = lapalma_setup $
                              , no_observer_metadata = no_observer_metadata $
                              , old_dir = old_dir 

  ;; Name of this program
  inam = red_subprogram(/low, calling = inam1)

  instrument = 'CHROMIS'
  
  if n_elements(ampm_cutoff) eq 0 then ampm_cutoff = '13:00:00'
  
  red_metadata_store, fname = work_dir + '/info/metadata.fits' $
                      , [{keyword:'INSTRUME', value:instrument $
                          , comment:'Name of instrument'} $
                         , {keyword:'TELCONFG', value:'Schupmann, imaging table', $
                            comment:'Telescope configuration'}]
  
  ;; Are there copies of/links to already summed calibration data?
  e_darksums   = file_test(work_dir + '/darks')
  e_flatsums   = file_test(work_dir + '/darks')
  e_pinhsums   = file_test(work_dir + '/pinhs')
  e_polcalsums = file_test(work_dir + '/polcal_sums')

  ;; Are there raw darks and flats?
  darksubdirs = red_find_instrumentdirs(root_dir, instrument, '*dark*' $
                                        , count = Ndarkdirs)
  flatsubdirs = red_find_instrumentdirs(root_dir, instrument, '*flat*' $
                                        , count = Nflatdirs)

  if Ndarkdirs eq 0 and ~e_darksums then begin
    print, 'No '+instrument+' darks were found. No setup generated.'
    return
  endif

  if Ndarkdirs gt 0 then begin
    ;; We want the W camera first!
    indx = sort(abs(reform(byte(strmid(file_basename(darksubdirs),0,1,/reverse))) - (byte('W'))[0]))
    darksubdirs = darksubdirs[indx]
  endif
  
  ;; Open two files for writing. Use logical unit Clun for a Config
  ;; file and Slun for a Script file.
  
  openw, Clun, work_dir + cfgfile,    /get_lun
  openw, Slun, work_dir + scriptfile, /get_lun

  ;; Specify the date in the config file, ISO format.
  print, 'Date'
  printf, Clun, '#'
  printf, Clun, '# --- Settings'
  printf, Clun, '#'
  printf, Clun,'isodate = '+isodate
  if isodate gt red_dates(tag = 'CHROMIS Ximea') then begin
    printf, Clun,'image_scale = 0.0355' ; Measured in August 2025.
  endif else begin
    printf, Clun,'image_scale = 0.0379' ; Measured in May 2016.
  endelse
  printf, Clun, 'diversity = 3.35e-3' ; Nominal value for 2016.

  if keyword_set(calibrations_only) then begin
    printf, Slun, 'a = chromisred("'+cfgfile+'",/dev)' 
  endif else begin
    printf, Slun, 'a = chromisred("'+cfgfile+'")' 
  endelse
  printf, Slun, 'root_dir = "' + root_dir + '"'

  ;; Specify default number of threads in script
  Ncpu = !cpu.hw_ncpu
  If Ncpu le 2 then Nthreads = 2 else Nthreads = round(Ncpu*.75) <20
  printf, Slun, 'nthreads='+strtrim(nthreads, 2)
  
  ;; Download SST log files and optionally some other data from the
  ;; web.
  if ~keyword_set(calibrations_only) && ~keyword_set(lapalma_setup) then begin
    print, 'Log files'
    printf, Clun, '#'
    printf, Clun, '# --- Download SST log files'
    printf, Clun, '#'
    printf, Slun
    printf, Slun, 'a -> download ; add ", /all" to get also HMI images and AR maps.'
  endif
  
;  ;; Analyze directories and produce r0 plots (optional)
;  if ~keyword_set(calibrations_only) && ~keyword_set(lapalma_setup) then begin
;    printf, Slun, '; a -> analyze_directories ; Time consuming, do it in a separate IDL session.'
;    printf, Slun, '; red_plot_r0 requires analyze_directories to have been run:'
;    printf, Slun, '; red_plot_r0, /plot8, /mark ; Plot r0 for the whole day.'
;    printf, Slun, '; red_plot_r0, /scan8 ; Plot r0 data for scans. '
;    printf, Slun, '; red_plot_r0, /plotstats ; Plot r0 statistics vs scan number.'
;    printf, Slun, '; a -> plot_pointing ; Plot pointing throughout the day.'
;  endif

  print, 'Cameras'
  printf, Clun, '#'
  printf, Clun, '# --- Cameras'
  printf, Clun, '#'

  if Ndarkdirs gt 0 then begin
    cams = file_basename(darksubdirs)
  endif else begin

    if ~e_darksums then stop

    dnames = file_search(work_dir+'/darks/cam*fits', count = Ndarks)
    if Ndarks eq 0 then stop
    
    cams = red_fitsgetkeyword_multifile(dnames, 'CAMERA', count = cnt)
    if total(cnt) eq 0 then stop

  endelse
  cams = red_uniquify(cams)

  ;; We want the WB camera first!
  indx = sort(abs(reform(byte(strmid(cams,0,1,/reverse)))  - (byte('W'))[0]))
  cams = cams[indx]
  
  is_wb = strmatch(cams, (instrument.tolower()).capwords()+'-[WD]')          
  is_pd = strmatch(cams, (instrument.tolower()).capwords()+'-D')
  
  printf, Clun, 'camera = '+cams, format='(a0)'
  printf, Clun, '#'
;  cams = red_uniquify(file_basename(file_search(root_dir+'/CHROMIS-data/*/Chromis-?')), count = Ncams)
;  cams = [cams[-1],cams[0],cams[1:-2]] ; First W, then D, then the others
;  for icam = 0, Ncams-1 do printf, Clun, 'camera = '+cams[icam]
;  printf, Clun, '#'

  ;; Based on the cameras, we can decide whether these data use the
  ;; polarimetric setup or not. (See also chromis::polarmetry_mode.)
  polarimetric_data = total(strmatch(cams,'*-[TR]')) gt 0
  
  ;; Orientation and rotation of WB camera, see IDL's rotate() and
  ;; offset_angle of red_lp_angles().

  year = strmid(isodate,0,4)
  case 1 of
    total(strmatch(['2016','2017'] $
                        , year)) : begin
      direction = 2    
      rotation = -42.0 
    end
    total(strmatch(['2018', '2019', '2020', '2021', '2022', '2023', '2024'] $
                        , year)) : begin
      direction = 1    
      rotation = -42.0 
    end
    year ge '2025' : begin
      direction = 3 
      rotation = -42.0
    end
    else : stop
  endcase

  if n_elements(direction) gt 0 then printf, Clun, 'direction = '+strtrim(direction, 2)
  if n_elements(rotation)  gt 0 then printf, Clun, 'rotation = '+strtrim(rotation, 2)
;  printf, Clun, 'direction = 2'  ; Orientation of WB camera, see IDL's rotate().
;  printf, Clun, 'rotation = -42.0' ; Rotation of reference camera, se offset_angle of red_lp_angles().
  printf, Clun, '#'
  
  printf, Clun, 'root_dir = ' + root_dir
  printf, Clun, '#'

  print, 'Output'
  printf, Clun, '#'
  printf, Clun, '# --- Output'
  printf, Clun, '#'
  printf, Clun, 'out_dir = ' + work_dir

  print, 'Darks'
  printf, Clun, '#'
  printf, Clun, '# --- Darks'
  printf, Clun, '#'

  if Ndarkdirs gt 0 then begin
    darkdirs = file_dirname(darksubdirs)
    darkdirs = darkdirs[uniq(darkdirs, sort(darkdirs))]
    Ndarkdirs = n_elements(darkdirs)

;    ;; Which darkdirs are collected in the AM?
;    dark_am = file_basename(darkdirs) lt ampm_cutoff

    printf, Slun
    ;; Loop over the darkdirs, write each to the config file and
    ;; to the script file
    for idir = 0, Ndarkdirs-1 do begin

      ;; Config file

      printf, Clun, 'dark_dir = '+red_strreplace(darkdirs[idir] $
                                                 , root_dir, '')

      ;; Script file
      
      if keyword_set(calibrations_only) || keyword_set(lapalma_setup) then begin
        ;; We want to output the summed data in timestamp directories
        ;; so we can handle multiple sets.
        outdir = 'darks/' + file_basename(darkdirs[idir])
        outdir_key = ', outdir="'+outdir+'"'
;        if keyword_set(calibrations_only) then calib_key = ', /softlink' else calib_key = ''
      endif else begin
        outdir_key = ''
;        calib_key = ''
      endelse
      calib_key = ', /softlink'
      printf, Slun, 'a -> sumdark, /sum_in_rdx, /check, dirs=root_dir+"' $
              + red_strreplace(darkdirs[idir], root_dir, '') + '"' $
              + ', nthreads=nthreads' $
              + outdir_key + calib_key
    endfor                      ; idir
  endif                         ; Ndarkdirs


  if Nflatdirs eq 0 && ~e_flatsums then begin
    print, 'No CHROMIS flats were found. Stop after summing darks.'
    free_lun, Clun
    free_lun, Slun
    return
  endif
 
  print, 'Flats'
  if Nflatdirs gt 0 && ~e_flatsums then begin
    ;; There are flats but no already summed flats!

    print, 'Flats'
    printf, Clun, '#'
    printf, Clun, '# --- Flats'
    printf, Clun, '#'

    ;; Directories with camera dirs below:
    flatdirs = file_dirname(flatsubdirs)
    flatdirs = flatdirs[uniq(flatdirs, sort(flatdirs))]
    Nflatdirs = n_elements(flatdirs)

    ;; CHROMIS file names aren't wheel-and-hrz coded from 2022-11-03 
    if isodate lt red_dates(tag = 'CHROMIS tuning metadata') then begin
      ;; Do the linedefs and hrz_calib thing for earlier observations
      red_download_linedefs, isodate, flatdirs, work_dir
      chromis_hrz_zeropoint, work_dir
    endif
    
    ;; Loop over the flatdirs, write each to the config file and
    ;; to the script file
    printf, Slun
    for idir = 0, Nflatdirs-1 do begin

      ;; Config file

      printf, Clun, 'flat_dir = '+red_strreplace(flatdirs[idir] $
                                                 , root_dir, '')

      ;; Script file
      
      ;; Look for wavelengths in those flatsubdirs that match
      ;; flatdirs[idir]! Also collect prefilters.
      indx = where(strmatch(flatsubdirs, flatdirs[idir]+'*'))
      fnames = file_search(flatsubdirs[indx]+'/sst_cam*', count = Nfiles)
      if Nfiles gt 0 then begin
        
        camdirs = strjoin(file_basename(flatsubdirs[indx]), ' ')
;        red_extractstates, fnames, /basename, pref = wls, is_wb = this_is_wb, is_pd = this_is_pd, wav = wav
        
        if isodate lt red_dates(tag = 'CHROMIS tuning metadata') then begin
          ;; Old filenames with wheel and hrz
          red_extractstates, fnames, /basename, pref = wls, is_wb = this_is_wb, is_pd = this_is_pd, wav = wav
        endif else begin
          ;; New filenames with proper tuning but NB files having WB
          ;; prefilter in the names and states.
          red_extractstates, fnames, /basename, wav = wav, fullstate = fullstate
          this_is_wb = strmatch(fnames, '*/Chromis-[WD]/*')          
          this_is_pd = strmatch(fnames, '*/Chromis-D/*')
          if this_is_wb[0] then wls = strmid(fullstate, 0, 4) else wls = strmid(wav, 0, 4)
        endelse

        indx = uniq(wls, sort(wls))
        wls = wls[indx]
        this_is_wb = this_is_wb[indx]
        this_is_pd = this_is_pd[indx]
        wls = wls[where(wls ne '')]
        
        
        if n_elements(wls) gt 0 then begin
          wavelengths = strjoin(wls, ' ')
          
          if keyword_set(calibrations_only) || keyword_set(lapalma_setup) then begin
            ;; We want to output the summed data in timestamp
            ;; directories so we can handle multiple sets.
            outdir = 'flats/' + file_basename(flatdirs[idir])
            outdir_key = ', outdir="'+outdir+'"'
            ;; If there are multiple dark directories, we want to use
            ;; the one that are nearest in time.
            tmp = min(abs(red_time2double(file_basename(flatdirs[idir])) - red_time2double(file_basename(darkdirs))), dindx)
            dark_timestamp_key = ", dark_timestamp = '" + file_basename(darkdirs[dindx]) + "'"
            if keyword_set(calibrations_only) then calib_key = ', /softlink, /store_rawsum' else calib_key = ''
          endif else begin
            outdir_key = ''
            dark_timestamp_key = ''
            calib_key = ''
          endelse
          
          ;; Print to script file
          printf, Slun, 'a -> sumflat, /sum_in_rdx, /check' $
                  + ', dirs=root_dir+"' + red_strreplace(flatdirs[idir], root_dir, '') + '"' $
                  + ', nthreads=nthreads' $
                  + outdir_key + dark_timestamp_key + calib_key $
                  + ' ; ' + camdirs+' ('+wavelengths+')'

          if this_is_wb[0] then red_append, prefilters, wls
          red_append, is_wb, this_is_wb
          red_append, is_pd, this_is_pd
        endif
      endif                     ; Nfiles
    endfor                      ; idir
  endif                         ; Nsubdirs
  
  Nprefilters = n_elements(prefilters)

  if Nprefilters eq 0 then begin
    ;; This can happen if flats were already summed in La Palma and
    ;; then deleted. Look for prefilters in the summed flats directory
    ;; instead.
    if file_test(work_dir + '/flats') then begin
      fdir = work_dir + '/flats/'
      ;;  spawn, 'ls '+work_dir+'/flats/cam*.flat | cut -d. -f2|sort|uniq', prefilters
    endif else if file_test(old_dir + '/flats') then begin
      fdir = old_dir + '/flats/'
      ;;   spawn, 'ls '+old_dir+'/flats/cam*.fits | cut -d. -f2|sort|uniq', prefilters
    endif
    fnames = file_search(fdir+'cam*fits', count = Nfiles)
    if Nfiles eq 0 then stop
    ;; New filenames with proper tuning but NB files having WB
    ;; prefilter in the names and states.
    red_extractstates, fnames, /basename, wav = wav, fullstate = fullstate, pref = prefilters
    Nprefilters = n_elements(prefilters)
    cameras = red_fitsgetkeyword_multifile(fnames, 'CAMERA', counts = cnt)
    is_wb = strmatch(cameras, '*-[WD]')          
    is_pd = strmatch(cameras, '*-D')
  endif 

  if Nprefilters eq 0 then stop

  indx = uniq(prefilters, sort(prefilters))
  is_wb = is_wb[indx]
  is_pd = is_pd[indx]
  prefilters = prefilters[indx]
  Nprefilters = n_elements(prefilters)

  if ~e_pinhsums then begin
    printf, Slun, ''
    print, 'Pinholes'
    printf, Clun, '#'
    printf, Clun, '# --- Pinholes'
    printf, Clun, '#'
    pinhdirs = file_search(root_dir+'/*pinh*/*', count = Ndirs, /fold)
    for idir = 0, Ndirs-1 do begin
      pinhsubdirs = file_search(pinhdirs[idir]+'/chromis*', count = Nsubdirs, /fold)
      if Nsubdirs gt 0 then begin
        printf, Clun, 'pinh_dir = '+red_strreplace(pinhdirs[idir], root_dir, '')
        
        if keyword_set(calibrations_only) || keyword_set(lapalma_setup) then begin
          ;; We want to output the summed data in timestamp directories
          ;; so we can handle multiple sets.
          outdir = 'pinhs/' + file_basename(pinhdirs[idir])
          outdir_key = ', outdir="'+outdir+'"'
          ;; If there are multiple dark directories, we want to use
          ;; the one that are nearest in time.
          tmp = min(abs(red_time2double(file_basename(pinhdirs[idir])) - red_time2double(file_basename(darkdirs))), dindx)
          dark_timestamp_key = ", dark_timestamp = '" + file_basename(darkdirs[dindx]) + "'"
          ;; If there are AM/PM flats directories, we want to take that
          ;; into account.
          flat_am = file_basename(flatdirs) lt ampm_cutoff
          pinh_am = file_basename(pinhdirs[idir]) lt ampm_cutoff
          case pinh_am of

            1 : begin
              ;; AM pinholes
              case total(flat_am) of

                0 : begin
                  ;; There are no AM flats, use PM flats
                  flat_timestamp_key = ", flat_timestamp = ['" $
                                       + strjoin(file_basename(flatdirs), "','") + "']"
                end

                else : begin
                  ;; There are AM flats
                  flat_timestamp_key = ", flat_timestamp = ['" $
                                       + strjoin(file_basename(flatdirs[where(flat_am)]), "','") + "']"
                end
                
              endcase           
            end

            else :  begin
              ;; PM pinholes
              case total(~flat_am) of

                0 : begin
                  ;; There are no PM flats, use AM flats
                  flat_timestamp_key = ", flat_timestamp = ['" $
                                       + strjoin(file_basename(flatdirs), "','") + "']"
                end

                else : begin
                  ;; There are PM flats
                  flat_timestamp_key = ", flat_timestamp = ['" $
                                       + strjoin(file_basename(flatdirs[where(~flat_am)]), "','") + "']"
                end
                
              endcase           
            end
            
          endcase
          
        endif else begin
          outdir_key = ''
          dark_timestamp_key = ''
          flat_timestamp_key = ''
        endelse

        printf, Slun, "a -> sumpinh, /sum_in_rdx, /pinhole_align" $
                + ', nthreads=nthreads' $
                + outdir_key + dark_timestamp_key + flat_timestamp_key $
                + ", dirs=root_dir+'" + red_strreplace(pinhdirs[idir], root_dir, '') + "'" 

      endif else begin

        pinhsubdirs = file_search(pinhdirs[idir]+'/*', count = Nsubdirs)
        for jdir = 0, Nsubdirs-1 do begin
          pinhsubsubdirs = file_search(pinhsubdirs[jdir]+'/chromis*' $
                                       , count = Nsubsubdirs, /fold)
          if Nsubsubdirs gt 0 then begin
            printf, Clun, 'pinh_dir = ' $
                    + red_strreplace(pinhsubdirs[jdir], root_dir, '')
            
            if keyword_set(calibrations_only) || keyword_set(lapalma_setup) then begin
              ;; For /calibrations_only we want to output the summed data in
              ;; timestamp directories so we can handle multiple sets.
              outdir = 'pinhs/' + file_basename(pinhsubdirs[jdir])
              outdir_key = ', outdir="'+outdir+'"'
              ;; If there are multiple dark directories, we want to use
              ;; the one that are nearest in time.
              tmp = min(abs(red_time2double(file_basename(pinhdirs[jdir])) - red_time2double(file_basename(darkdirs))), dindx)
              dark_timestamp_key = ", dark_timestamp = '" + file_basename(darkdirs[dindx]) + "'"
              ;; If there are AM/PM flats directories, we want to take that
              ;; into account.
              flat_am = file_basename(flatdirs) lt ampm_cutoff
              stop              ; Does this case still happen?
            endif else begin
              outdir_key = ''
              dark_timestamp_key = ''
              flat_timestamp_key = ''
            endelse 

            printf, Slun, "a -> sumpinh, /sum_in_rdx, /pinhole_align" $
                    + ', nthreads=nthreads' $
                    + outdir_key + dark_timestamp_key + flat_timestamp_key $
                    + ", dirs=root_dir+'" +  red_strreplace(pinhsubdirs[jdir], root_dir, '') + "'" 
            
          endif
        endfor                  ; jdir
      endelse
    endfor                      ; idir
  endif
  
  if ~keyword_set(calibrations_only) && ~keyword_set(lapalma_setup) then begin  
    printf, Slun, ''
    if isodate gt red_dates(tag = 'CHROMIS Ximea') then begin
      printf, Slun, 'a -> pinholecalib'
    endif else begin
      ;; The new polywarp mechanism does not work with the aspect
      ;; ratio of the old CHROMIS cameras. Yet. Edit this when fixed!
      printf, Slun, 'a -> pinholecalib, /no_polywarp, /verify'
    endelse
;    printf, Slun, 'a -> diversitycalib'
    printf, Slun, ''
  endif

  if polarimetric_data && ~e_polcalsums then begin
    
    print, 'Polcal'
    printf, Clun, '#'
    printf, Clun, '# --- Polcal'
    printf, Clun, '#'
    polcalsubdirs = red_find_instrumentdirs(root_dir, instrument, instrument+'-polc*' $
                                            , count = Npolcalsubdirs)
    
    if Npolcalsubdirs gt 0 then begin
      
      polcaldirs = file_dirname(polcalsubdirs)  
      polcaldirs = polcaldirs[uniq(polcaldirs,sort(polcaldirs))]
      Npolcaldirs = n_elements(polcaldirs)
      polprefs = strarr(Npolcaldirs)

      printf, Slun
      
      for idir = 0, Npolcaldirs-1 do begin
        polcalsubdirs = file_search(polcaldirs[idir]+'/'+instrument+'*' $
                                    , count = Nsubdirs, /fold)
        if Nsubdirs gt 0 then begin
          printf, Clun, 'polcal_dir = ' $
                  + red_strreplace(polcaldirs[idir], root_dir, '')

          if keyword_set(calibrations_only) || keyword_set(lapalma_setup) then begin

            ;; We want to output the summed data in timestamp
            ;; directories so we can handle multiple sets.
            outdir = 'polcal_sums/' + file_basename(polcaldirs[idir])
            outdir_key = ', outdir="'+outdir+'"'
            
          endif else begin
            outdir_key = ''
          endelse

          printf, Slun, 'a -> sumpolcal, /sum_in_rdx, /check, dirs=root_dir+"' $
                  + red_strreplace(polcaldirs[idir], root_dir, '')+'"' $
                  + ', nthreads=nthreads' $
                  + outdir_key
          ;; The prefilter is not part of the path. Try to get it from
          ;; the first data file in the directory.
          files = file_search(polcalsubdirs[0]+'/*', count = Npolfiles)
          if Npolfiles gt 0 then begin
            hh = red_readhead(files[0])
            ;; Unlike for CRISP, the WB filter in FILTER1 is not the same as the NB filter.
            ;;  polprefs[idir] = strtrim(fxpar(hh, 'FILTER1'), 2)
            state = fxpar(hh, 'STATE', count = cnt)
            if cnt eq 0 then stop
            polprefs[idir] = (strsplit(state, '_', /extract))[1]
          endif
        endif
      endfor                    ; idir

      if ~keyword_set(calibrations_only) && ~keyword_set(lapalma_setup) then begin  
        for ipref = 0, n_elements(polprefs)-1 do begin
          printf, Slun, "a -> polcalcube, pref='" + polprefs[ipref] + "'" $
                  + ", nthreads=nthreads"
          printf, Slun, "a -> polcal, pref='" + polprefs[ipref] + "'" $
                  + ", nthreads=nthreads"
        endfor                  ; ipref
      endif
      
    endif else if n_elements(old_dir) gt 0 then begin

      ;; Need to find polprefs and Npolcaldirs from the summed data
      pdir = old_dir + '/polcal_sums/'+instrument+'-T/'
      pfiles = file_search(pdir+'cam*fits', count = Npfiles)
      if Npfiles eq 0 then Npolcaldirs = 0 else begin
        red_extractstates, pfiles, /basename, pref = polprefs
        polprefs = polprefs[uniq(polprefs,sort(polprefs))]
      endelse
      
      for ipref = 0, n_elements(polprefs)-1 do begin
        printf, Slun, "a -> polcalcube, pref='" + polprefs[ipref] + "'" $
                + ", nthreads=nthreads"
        printf, Slun, "a -> polcal, pref='" + polprefs[ipref] + "'" $
                + ", nthreads=nthreads"
      endfor                    ; ipref
    endif

  endif


  if ~keyword_set(calibrations_only) && ~keyword_set(lapalma_setup) then begin  
    ;;if ~is_wb[ipref] then
    printf, Slun
    printf, Slun, "a -> prepflatcubes"
  endif
  
  if ~keyword_set(calibrations_only) && ~keyword_set(lapalma_setup) then begin  
    printf, Slun, ''
    printf, Slun, '; The fitgains step requires the user to look at the fit and determine'
    printf, Slun, '; whether you need to use different keyword settings.'
    printf, Slun, '; Then, if you have already run makegains, rerun it.'
    printf, Slun, "a -> fitgains, rebin=800L, Niter=3L, Nthreads=nthreads, Npar=5L, res=res"
;    printf, Slun, '; If you need per-pixel reflectivities for your analysis'
;    printf, Slun, '; (e.g. for atmospheric inversions) you can set the /fit_reflectivity'
;    printf, Slun, '; keyword:'
;    printf, Slun, '; a -> fitgains, npar = 3, res=res, /fit_reflectivity  '
;    printf, Slun, '; However, running without /fit_reflectivity is safer. In should not'
;    printf, Slun, '; be used for chromospheric lines like 6563 and 8542.'
;    printf, Slun, '; Sometimes you need to add a few spline nodes in order to make the fits work,'
;    printf, Slun, '; particularly just to the red and to the blue of the densely sampled region and'
;    printf, Slun, '; also in blends if they are not propely sampled.'
;    printf, Slun, '; As an example, to do this for 3969 Å, 12.00ms_G10.00 data, do something like'
;    printf, Slun, '; the following:'
;    printf, Slun, "; restore,'flats/spectral_flats/camXXX_12.00ms_G10.00_3969_flats.sav'; Read flats cube"
;    printf, Slun, '; myg = [wav*1.d10, -0.765d0, 0.740d0] ; Add wavelength points'
;    printf, Slun, '; myg=myg[sort(myg)]  ; Sort the wavelength points'
;    printf, Slun, '; a->fitgains, rebin=800L, niter=3L, nthreads=12L, res=res, npar=5L, myg=myg ; Run the fitgain step with the added wavelength points'
;    printf, Slun, '; Then, if you have already run makegains, rerun it.'
    
    printf, Slun, ''

    printf, Slun, "a -> makegains, smooth=3.0, min=0.1, max=4.0, bad=1.0, nthreads = nthreads"

  endif
  

  
  print, 'Prefilter scan'
  printf, Clun, '#'
  printf, Clun, '# --- Prefilter scan'
  printf, Clun, '#'
  Npfs = 0
  pfscandirs = file_search(root_dir+'/*pfscan*/*', count = Ndirs, /fold)
  for i = 0, Ndirs-1 do begin
    pfscansubdirs = file_search(pfscandirs[i]+'/chromis*', count = Nsubdirs, /fold)
    if Nsubdirs gt 0 then begin
      printf, Clun, '# pfscan_dir = '+red_strreplace(pfscandirs[i], root_dir, '')
      Npfs += 1
    endif else begin
      pfscansubdirs = file_search(pfscandirs[i]+'/*', count = Nsubdirs)
      for j = 0, Nsubdirs-1 do begin
        pfscansubsubdirs = file_search(pfscansubdirs[j]+'/chromis*' $
                                       , count = Nsubsubdirs, /fold)
        if Nsubsubdirs gt 0 then begin
          printf, Clun, '# pfscan_dir = ' $
                  + red_strreplace(pfscansubdirs[j], root_dir, '')
          Npfs += 1
        endif                   ; Nsubsubdirs
      endfor                    ; j
    endelse                     ; Nsubdirs
  endfor                        ; i
  ;; If we implement dealing with prefilter scans in the pipeline,
  ;; here is where the command should be written to the script file.

  
  if keyword_set(calibrations_only) then begin
    ;; We don't need science data for the calibrations_only setup. 
    free_lun, Clun
    ;; We'll end the script file with an end statement so it
    ;; can be run with .run or .rnew.
    printf, Slun, 'end'
    free_lun, Slun
    return
  endif

  print, 'Science'
  printf, Clun, '#'
  printf, Clun, '# --- Science data'
  printf, Clun, '# '

  ;; Exclude directories known not to have science data
  red_append, nonsciencedirs, pinhdirs
  red_append, nonsciencedirs, pfscandirs
  red_append, nonsciencedirs, darkdirs
  red_append, nonsciencedirs, flatdirs
  red_append, nonsciencedirs, polcaldirs
  ;; Sometimes the morning calibrations data were not deleted, so they
  ;; have to be excluded too.
  calibsubdirs = red_find_instrumentdirs(root_dir, instrument, '*calib' $
                                         , count = Ncalibdirs)
  if Ncalibdirs gt 0 then begin
    calibdirs = file_dirname(calibsubdirs)
    calibdirs = calibdirs[uniq(calibdirs, sort(calibdirs))]
    red_append, nonsciencedirs, calibdirs
  end
  sciencedirs = file_search(root_dir+'/*/*', count = Ndirs)

  for i = 0, Ndirs-1 do begin

    if total(sciencedirs[i] eq nonsciencedirs) eq 0 then begin
      sciencesubdirs = file_search(sciencedirs[i]+'/chromis*' $
                                   , count = Nsubdirs, /fold)
      if Nsubdirs gt 0 then begin
        red_append, dirarr, red_strreplace(sciencedirs[i], root_dir, '')
      endif else begin
        sciencesubdirs = file_search(sciencedirs[i]+'/*', count = Nsubdirs)
        for j = 0, Nsubdirs-1 do begin
          sciencesubsubdirs = file_search(sciencesubdirs[j]+'/chromis*' $
                                          , count = Nsubsubdirs, /fold)
          if Nsubsubdirs gt 0 then begin
            red_append, dirarr, red_strreplace(sciencesubdirs[j], root_dir, '')
          endif
        endfor                  ; j
      endelse 
    endif
  endfor
  if n_elements(dirarr) gt 0 then printf, Clun, "data_dir = ['"+strjoin(dirarr, "','")+"']"

  printf, Slun
  printf, Slun, ";Assess your science data:"
  printf, Slun, ";a -> summary"
  printf, Slun, "; or perhaps your selection of the following:" 
  printf, Slun, ";a -> summarize_datadir"
  printf, Slun, ";a -> quicklook, /core_and_wings, datasets = '*'        ; NB quicklook"
  printf, Slun, ";a -> quicklook_mosaic, /core_and_wings, datasets = '*' ; NB mosaic quicklook"
  printf, Slun, ";a -> quicklook, cam='"+cams[0]+"', datasets = '*', use_states = '*' ; WB quicklook"

  if keyword_set(lapalma_setup) then begin

    printf, Slun, '; Quicklook will sum its own flats, independent of the ones in the flats/hh:mm:ss/ subdirectories.'
    printf, Slun, 'a -> summarize_datadir'
    
    free_lun, Slun, Clun

    return
    
  endif
  
  printf, Slun, ''
  printf, Slun, 'a -> link_data' 

  printf, Slun, ''
  printf, Slun, ';; -----------------------------------------------------'
  printf, Slun, ';; This is how far we should be able to run unsupervised'
  printf, Slun, 'stop'          
  printf, Slun, ''

  ;;  for ipref = 0, Nprefilters-1 do printf, Slun, "a -> fitprefilter, /mask ;, pref = '"+prefilters[ipref]+"'"
  printf, Slun, "a -> fitprefilter, /mask; , /hints, dir='10:02:45'"
  printf, Slun, ''
  printf, Slun, "a -> fit_wb_diskcenter, tmax='13:00'; for PM data instead tmin='13:00'"

  printf, Slun, ''
  for ipref = 0, Nprefilters-1 do begin
    if is_wb[ipref] then begin
      if polarimetric_data then begin
        printf, Slun, "a -> prepmomfbd" $
;              + ", date_obs='" + isodate + "'" $
                + ", Nremove=1" $
                + ", Nmodes=60" $
                + ", numpoints=128" $
;            + ", margin=5 " $
                + ", global_keywords=['FIT_PLANE']" $
                + ", maxshift=45" $
                + ", /wb_states" $
;              + ", /redux" $
                + ", /unpol" $
                + ", extraclip = [50,50,40,20]" $ ; [75,125,15,15]" $
                + ", pref='" + prefilters[ipref] + "'" $
                + ", dirs=['"+strjoin(file_basename(dirarr), "','")+"']"
      endif else begin
        printf, Slun, "a -> prepmomfbd" $
;              + ", date_obs='" + isodate + "'" $
                + ", Nremove=1" $
                + ", Nmodes=60" $
                + ", numpoints=128" $
;            + ", margin=5 " $
                + ", global_keywords=['FIT_PLANE']" $
                + ", maxshift=45" $
                + ", /wb_states" $
;              + ", /redux" $
                + ", /unpol" $
                + ", extraclip = [50,50,40,20]" $ ; = [75,125,15,15]" $
                + ", pref='" + prefilters[ipref] + "'" $
                + ", dirs=['"+strjoin(file_basename(dirarr), "','")+"']"
      endelse
    endif
  endfor                        ; ipref ;

  printf, Slun, ''
  printf, Slun, ';; Run MOMFBD outside IDL.'
  printf, Slun, ''

  printf, Slun, ';; Post-MOMFBD stuff:' 
  printf, Slun, "a -> align_continuum"
  printf, Slun
  printf, Slun, "a -> make_scan_cube, 'momfbd/.../cfg/results/', /autocrop, scannos = '69', nthreads=nthreads"
  printf, Slun, "a -> fitscube_wcs_improve_spatial, 'cubes_scan/nb....fits' ; If suitable target"
  printf, Slun, "; or "
  printf, Slun, "a -> make_wb_cube, 'momfbd/.../cfg/results/', nthreads=nthreads, /interactive, /autocrop, /align_interactive, outpath = wbpath"
  printf, Slun, "; If not continuing immediately, set something like wbpath='cubes_wb/wb....fits' "
  printf, Slun, "a -> fitscube_wcs_improve_spatial, wbpath ; If suitable target"
  printf, Slun, "; If called without the wbpath, make_nb_cube lets you select one of the existing wb cubes"
  printf, Slun, "a -> make_nb_cube, wbpath, nthreads=nthreads"
  printf, Slun, "; or "
  printf, Slun, "a -> make_mos_cube, 'momfbd/.../cfg/', nthreads=nthreads"

  
  free_lun, Clun
  free_lun, Slun

  ;; Do something about OBSERVER metadata keyword
  if ~keyword_set(no_observer_metadata) then begin
    ;; See if we can find some metadata by looking in the raw data dirs.
        chromis_data_dirs = root_dir + '/' + dirarr + '/Chromis-W/'
        ;; Pick the first file in each.
        chromis_data_files = file_search(chromis_data_dirs+'/*00000_0000000*fits', count = Nfiles) 
        ;; Now look for OBSERVER keywords
        observers = strarr(Nfiles)
        for ifile = 0, Nfiles-1 do begin
      red_progressbar, ifile, Nfiles, 'Looking in CHROMIS data for OBSERVER keyword'
      observers[ifile] = red_fitsgetkeyword(chromis_data_files[ifile], 'OBSERVER')
    endfor
    observers = ['', observers] ; empty default means no OBSERVER keyword in metadata
    indx = uniq(observers, sort(observers))
    print
    if n_elements(indx) gt 1 then begin
      print, inam + ' : Found OBSERVER keyword(s) in the CHROMIS raw data. All is well.'
    endif else begin      
      print, inam + ' : Found no OBSERVER metadata in the CHROMIS raw data.'
      observer = ''
      read, 'Add names for that keyword for the CHROMIS workdir or hit return: ', observer
      ;; Write it to the metadata file
      print
      if observer ne '' then begin
        print, inam+' : Adding to CHROMIS metadata, OBSERVER = '+observer
        red_metadata_store, fname = work_dir + '/info/metadata.fits' $
                            , [{keyword:'OBSERVER', value:observer $
                                , comment:'Observer name(s)'}]
      endif else begin
        print, inam+' : No OBSERVER keyword in CHROMIS metadata.'
      endelse
      print, inam+' : Edit '+work_dir + '/info/metadata.fits if you need to change this.'
    endelse
    print
  endif
  
end
