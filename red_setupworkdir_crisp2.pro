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
;    2022-08-01 : MGL. New version for CRISP2 (and old CRISP with new
;                 cameras) based on the chromis version.
; 
;    2024-07-04 : MGL. If old_dir is given, copy also polcal_sums.
; 
;    2025-01-14 : MGL. No periodic filter needed when we have polcal
;                 flats.
;
;    2025-04-03: MGL. New keywords lapalma_setup, ampm_cutoff. 
; 
;-
pro red_setupworkdir_crisp2, work_dir, root_dir, cfgfile, scriptfile, isodate $
                             , ampm_cutoff = ampm_cutoff $
                             , calibrations_only = calibrations_only $
                             , lapalma_setup = lapalma_setup $
                             , no_observer_metadata = no_observer_metadata $
                             , old_dir = old_dir 

  ;; Name of this program
  inam = red_subprogram(/low, calling = inam1)

  if n_elements(ampm_cutoff) eq 0 then ampm_cutoff = '13:00:00'

  ;; Find out if this is CRISP2 or old CRISP from the data dirs.
  dirs = file_search(root_dir + '/CRISP2-*', count = Ndirs)
  if Ndirs gt 0 then instrument = 'CRISP2' else instrument = 'CRISP'

  ;; Note: "instrument" is in all caps here. In most of the rest of
  ;; the pipeline, only the initial letter is capitalized.

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
  darksubdirs = red_find_instrumentdirs(root_dir, instrument, instrument+'-dark*' $
                                        , count = Ndarkdirs)
  flatsubdirs = red_find_instrumentdirs(root_dir, instrument, instrument+'-flat*' $
                                        , count = Nflatdirs)

  if Ndarkdirs eq 0 and ~e_darksums then begin
    print, 'No '+instrument+' darks were found. No setup generated.'
    return
  endif
  
  ;; We want the W camera first!
  indx = sort(abs(reform(byte(strmid(file_basename(darksubdirs),0,1,/reverse))) - (byte('W'))[0]))
  darksubdirs = darksubdirs[indx]
  
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
  case instrument of
    'CRISP' : begin
      printf, Clun, 'image_scale = 0.044' ; Get from pinhole calibration
      printf, Clun, 'diversity = 4.0e-3'  ; Nominal value for Sep 2023
    end
    'CRISP2' : begin
      printf, Clun, 'image_scale = 0.050' ; Get from pinhole calibration
      printf, Clun, 'diversity = 3.1e-3'  ; Nominal value for Aug 2025
    end
    else : stop
  endcase
  
  if keyword_set(calibrations_only) then begin
    printf, Slun, 'a = crisp2red("'+cfgfile+'",/dev)' 
  endif else begin
    printf, Slun, 'a = crisp2red("'+cfgfile+'")' 
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
  
  ;; Analyze directories and produce r0 plots (optional)
  if ~keyword_set(calibrations_only) && ~keyword_set(lapalma_setup) then begin
    printf, Slun, '; a -> analyze_directories ; Time consuming, do it in a separate IDL session.'
    printf, Slun, '; red_plot_r0 requires analyze_directories to have been run:'
    printf, Slun, '; red_plot_r0, /plot8, /mark ; Plot r0 for the whole day.'
    printf, Slun, '; red_plot_r0, /scan8 ; Plot r0 data for scans. '
    printf, Slun, '; red_plot_r0, /plotstats ; Plot r0 statistics vs scan number.'
    printf, Slun, '; a -> plot_pointing ; Plot pointing throughout the day.'
  endif

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

  ;; Orientation and rotation of WB camera, see IDL's rotate()
  ;; and  offset_angle of red_lp_angles(). 
  case instrument of
    'CRISP' : begin
      direction = 6
      rotation = -42.0 
    end
    'CRISP2' : begin
      direction = 6
      rotation = -42.0 
    end
    else : stop
  endcase
  if n_elements(direction) gt 0 then printf, Clun, 'direction = '+strtrim(direction, 2)
  if n_elements(rotation)  gt 0 then printf, Clun, 'rotation = '+strtrim(rotation, 2)
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

  if Ndarkdirs gt 0 && ~e_darksums then begin
    darkdirs = file_dirname(darksubdirs)
    darkdirs = darkdirs[uniq(darkdirs, sort(darkdirs))]
    Ndarkdirs = n_elements(darkdirs)

    printf, Slun
    ;; Loop over the darkdirs, write each to the config file and
    ;; to the script file
    for idir = 0, Ndarkdirs-1 do begin

      ;; Config file
      printf, Clun, 'dark_dir = '+red_strreplace(darkdirs[idir] $
                                                 , root_dir, '')

      ;; Script file
      
      if keyword_set(calibrations_only) || keyword_set(lapalma_setup) then begin
        ;; For /calibrations_only we want to output the summed data in
        ;; timestamp directories so we can handle multiple sets.
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
  endif                         ; Nsubdirs, ~e_darksums 


  if Nflatdirs eq 0 && ~e_flatsums then begin
    print, 'No '+instrument+' flats were found. Stop after summing darks.'
    free_lun, Clun
    free_lun, Slun
    return
  endif
  
  
  print, 'Flats'
  if Nflatdirs gt 0 && ~e_flatsums then begin
    ;; There are flats but no already summed flats!

    printf, Clun, '#'
    printf, Clun, '# --- Flats'
    printf, Clun, '#'

    ;; Directories with camera dirs below:
    flatdirs = file_dirname(flatsubdirs)
    flatdirs = flatdirs[uniq(flatdirs, sort(flatdirs))]
    Nflatdirs = n_elements(flatdirs)

    ;; Loop over the flatdirs, write each to the config file and
    ;; to the script file
    printf, Slun
    for idir = 0, Nflatdirs-1 do begin

      ;; Config file

      printf, Clun, 'flat_dir = ' + red_strreplace(flatdirs[idir] $
                                                   , root_dir, '')

      ;; Script file
      
      ;; Look for wavelengths in those flatsubdirs that match
      ;; flatdirs[idir]! Also collect prefilters.
      indx = where(strmatch(flatsubdirs, flatdirs[idir]+'*'))
      fnames = file_search(flatsubdirs[indx]+'/sst_cam*', count = Nfiles)
      if Nfiles gt 0 then begin
        
        camdirs = strjoin(file_basename(flatsubdirs[indx]), ' ')
        red_extractstates, fnames, /basename, pref = wls, is_wb = this_is_wb, is_pd = this_is_pd
        indx = uniq(wls, sort(wls))
        wls = wls[indx]
        this_is_wb = this_is_wb[indx]
        this_is_pd = this_is_pd[indx]
        wls = wls[WHERE(wls ne '')]

        if n_elements(wls) gt 0 then begin
          wavelengths = strjoin(wls, ' ')
          
          if keyword_set(calibrations_only) || keyword_set(lapalma_setup) then begin
            ;; We want to output the summed data in timestamp
            ;; directories so we can handle multiple sets.
            outdir = 'flats/' + file_basename(flatdirs[idir])
            outdir_key = ', outdir="'+outdir+'"'
            ;; If there are multiple dark directories, we want to use
            ;; the one that are nearest in time.
            tmp = min(abs(red_time2double(file_basename(flatdirs[idir])) $
                          - red_time2double(file_basename(darkdirs))), dindx)
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

          red_append, prefilters, wls
        endif
      endif                     ; Nfiles
    endfor                      ; idir
  endif                         ; Nsubdirs && ~e_flatsums
  
  Nprefilters = n_elements(prefilters)
  
  if Nprefilters eq 0 then begin
    ;; This can happen if flats were already summed in La Palma and
    ;; then deleted. Look for prefilters in the summed flats directory
    ;; instead.
    if file_test(work_dir + '/flats') then begin
;      fdir = work_dir + '/flats/'
      fnames=file_search(work_dir+'/flats/','*.fits', count = Nsearch)
      if Nsearch gt 0 then begin
        red_extractstates, fnames, /basename, pref = prefilters, cam = cam
        prefilters = red_uniquify(prefilters)
        Nprefilters = n_elements(prefilters)
      endif
    endif else if file_test(old_dir + '/flats') then begin
      ;; Recursively search for flats fits files.
      sss=file_search(old_dir+'/flats/','*.fits', count = Nsearch)
      if Nsearch gt 0 then begin
        fdirs = red_uniquify(file_dirname(sss))
        fnames = file_search(fdirs+'/cam*fits', count = Nfiles)
        red_extractstates, fnames, /basename, pref = prefilters, cam = cam
        prefilters = red_uniquify(prefilters)
        Nprefilters = n_elements(prefilters)
      endif
    endif
  endif
  
  if Nprefilters eq 0 then stop

;  indx = uniq(prefilters, sort(prefilters))
;  prefilters = prefilters[indx]
;  Nprefilters = n_elements(prefilters)

  prefilters = red_uniquify(prefilters, count = Nprefilters)
  
  if ~e_pinhsums then begin
    printf, Slun, ''
    print, 'Pinholes'
    printf, Clun, '#'
    printf, Clun, '# --- Pinholes'
    printf, Clun, '#'
    pinhdirs = file_search(root_dir+'/*pinh*/*', count = Ndirs, /fold)
    for idir = 0, Ndirs-1 do begin
      pinhsubdirs = file_search(pinhdirs[idir]+'/'+instrument+'*', count = Nsubdirs, /fold)
      if Nsubdirs gt 0 then begin
        printf, Slun
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

;;        for ipref = 0, Nprefilters-1 do begin
        printf, Slun, "a -> sumpinh, /sum_in_rdx, /pinhole_align" $
                  + ', nthreads=nthreads' $
                  + outdir_key + dark_timestamp_key + flat_timestamp_key $
                  + ", dirs=root_dir+'" + red_strreplace(pinhdirs[idir], root_dir, '') ;;$
;;                  + "', pref='"+prefilters[ipref]+"'"
;;        endfor                  ; ipref
        endif else begin
        pinhsubdirs = file_search(pinhdirs[idir]+'/*', count = Nsubdirs)
        printf, Slun
        for jdir = 0, Nsubdirs-1 do begin
          pinhsubsubdirs = file_search(pinhsubdirs[jdir]+'/'+instrument+'*' $
                                       , count = Nsubsubdirs, /fold)
          if Nsubsubdirs gt 0 then begin
            printf, Clun, 'pinh_dir = ' $
                    + red_strreplace(pinhsubdirs[jdir], root_dir, '')
            
            if keyword_set(calibrations_only) then begin
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

            ;;for ipref = 0, Nprefilters-1 do begin
            printf, Slun, "a -> sumpinh, /sum_in_rdx, /pinhole_align" $
                    + ', nthreads=nthreads' $
                    + outdir_key + dark_timestamp_key + flat_timestamp_key $
                    + ", dirs=root_dir+'" +  red_strreplace(pinhsubdirs[jdir], root_dir, '') ;;$
            ;;        + "';, pref='" + prefilters[ipref]+"'" 
            ;;endfor              ; ipref
          endif
        endfor                  ; jdir
      endelse
    endfor                      ; idir
  endif
  
  if ~keyword_set(calibrations_only) && ~keyword_set(lapalma_setup) then begin  
    printf, Slun, ''
    printf, Slun, 'a -> pinholecalib'
    printf, Slun, ''
;    printf, Slun, 'a -> diversitycalib'
  endif

;  if e_polcalsums then begin
;    fnames = file_search(work_dir + '/polcal_sums/*.fits', count = cnt)
;    red_extractstates, fnames, /basename, pref = polprefs, cam = cam
;    polprefs = red_uniquify(polprefs)
;    Npolprefs = n_elements(polprefs)
;  endif
  
  if ~e_polcalsums then begin
    
;    fnames = file_search(work_dir+'/polcal_sums/*/*.fits', count = Nsearch)
;    red_extractstates, fnames, /basename, pref = polprefs, cam = cam
;    polprefs = red_uniquify(polprefs)
;    Npolprefs = n_elements(polprefs)
;
;  endif else begin
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
      Npolprefs = n_elements(polprefs)
      
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
            polprefs[idir] = strtrim(fxpar(hh, 'FILTER1'), 2)
          endif
        endif
      endfor                    ; idir
    endif
  endif                         ; ~e_polcalsums

  ;; If there are polcal data, we need to run polcalcube and polcal
  ;; (and possibly make the periodic filter).
  if ~keyword_set(calibrations_only) && ~keyword_set(lapalma_setup) && ~e_polcalsums then begin  
    for ipref = 0, Npolprefs-1 do begin
      printf, Slun, "a -> polcalcube, pref='" + polprefs[ipref] + "'" $
              + ", nthreads=nthreads"
      printf, Slun, "a -> polcal, pref='" + polprefs[ipref] + "'" $
              + ", nthreads=nthreads"
      if isodate ge red_dates(tag = 'polcal flats', explanation = explanation) then begin
        print, explanation
        printf, Slun, "; The periodic filter should not be necessary when we have polcal flats"
        printf, Slun, "; a -> make_periodic_filter,'" + polprefs[ipref] + "'"
      endif else begin
        printf, Slun, "a -> make_periodic_filter,'" + polprefs[ipref] + "'"
      endelse
    endfor                      ; ipref
  endif else begin

    ;; Need to find polprefs and Npolcaldirs from the summed data
    pdir = file_search(work_dir + '/polcal_sums/'+instrument+'-T', /fold)
    pfiles = file_search(pdir+'/cam*fits', count = Npfiles)
    if Npfiles eq 0 then Npolcaldirs = 0 else begin
      red_extractstates, pfiles, /basename, pref = polprefs
      polprefs = red_uniquify(polprefs)
      Npolprefs = n_elements(polprefs)
    endelse
    
    for ipref = 0, Npolprefs-1 do begin
      printf, Slun, "a -> polcalcube, pref='" + polprefs[ipref] + "'" $
              + ", nthreads=nthreads"
      printf, Slun, "a -> polcal, pref='" + polprefs[ipref] + "'" $
              + ", nthreads=nthreads"
      if isodate ge red_dates(tag = 'polcal flats', explanation = explanation) then begin
        print, explanation
        printf, Slun, "; The periodic filter should not be necessary when we have polcal flats"
        printf, Slun, "; a -> make_periodic_filter,'" + polprefs[ipref] + "'"
      endif else begin
        printf, Slun, "a -> make_periodic_filter,'" + polprefs[ipref] + "'"
      endelse
    endfor                      ; ipref

  endelse

  if ~keyword_set(calibrations_only) && ~keyword_set(lapalma_setup) then begin  
    for ipref = 0, Nprefilters-1 do begin
      printf, Slun, "a -> prepflatcubes, pref='"+prefilters[ipref]+"'"
    endfor                      ; ipref
  endif

  if ~keyword_set(calibrations_only) && ~keyword_set(lapalma_setup) then begin  
    printf, Slun, ''
    printf, Slun, '; The fitgains step requires the user to look at the fit and determine'
    printf, Slun, '; whether you need to use different keyword settings.'
    printf, Slun, '; Then, if you have already run makegains, rerun it.'
    for ipref = 0, Nprefilters-1 do begin
      printf, Slun, "a -> fitgains, rebin=800L, Niter=3L, Nthreads=nthreads, Npar=5L, res=res, pref='" $
              + prefilters[ipref] + "'"
    endfor                      ; ipref

    printf, Slun, ''

    case instrument of
      'CRISP'  : printf, Slun, "a -> makegains, smooth=3.0, min=0.1, max=4.0, bad=1.0, nthreads = nthreads"
      'CRISP2' : begin
        printf, Slun, "a -> makegains, smooth=3.0, min=0.1, max=4.0, bad=1.0, nthreads = nthreads"
        ;; CRISP2 WB flats have a patch that makes a "hole" in the
        ;; gains if not taken special care of.
        printf, Slun, "wbflats = file_search('flats/'"+cams[0]+"'_*fits')   ; WB flats special" 
        printf, Slun, "a -> makegains, smooth=3.0, min=0.1, max=50.0, bad=1.0, flatmin=0.01, nthreads = nthreads,files=wbflats"
      end
      else : stop
    endcase
    
  endif


  
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

  sciencesubdirs = red_find_instrumentdirs(root_dir, instrument $
                                           , instrument+['-data*', '-mosaic*'] $
                                           , count = Nsubdirs)
  
  if Nsubdirs gt 0 then begin
    sciencedirs = file_dirname(sciencesubdirs)
    dirarr = red_strreplace(sciencedirs[uniq(sciencedirs, sort(sciencedirs))], root_dir, '')
    printf, Clun, "data_dir = ['"+strjoin(dirarr, "','")+"']"
  endif

  printf, Slun
  printf, Slun, '; Quicklook will sum its own flats if needed, independent of the ones in the flats/hh:mm:ss/ subdirectories.'
  printf, Slun, 'a -> quicklook, /core_and_wings        ; NB quicklook'
  printf, Slun, 'a -> quicklook_mosaic, /core_and_wings ; NB mosaic quicklook'
  printf, Slun, "a -> quicklook, cam='"+cams[0]+"'      ; WB quicklook"

  if keyword_set(lapalma_setup) then begin

    ;; For a La Palma setup, return after quicklook.
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
  printf, Slun, "a -> fit_wb_diskcenter, tmax='13:00'; for PM data instead tmin='13:00'"

  printf, Slun, '; If MOMFBD has problems near the edges, try increasing the margin when calling prepmomfbd.'
  for ipref = 0, Nprefilters-1 do begin
    ;; The number of frames to be discarded by momfbd from every state
    ;; is 2 for data without polarimetry and 1 for data with
    ;; polarimetry. This is because with polarimetry, 1 discarded
    ;; frame from every state means 4 frames are actually discarded
    ;; from each tuning because the LC is changing faster than the
    ;; tuning.
    if instrument eq 'CRISP' then begin
      printf, Slun, "a -> sum_data_intdif, pref = '" + prefilters[ipref] $
              + "', cam = 'Crisp-T', /verbose, /show, /overwrite " $
              + ', nthreads=nthreads' $
              + " ; /all"
      printf, Slun, "a -> sum_data_intdif, pref = '" + prefilters[ipref] $
              + "', cam = 'Crisp-R', /verbose, /show, /overwrite " $
              + ', nthreads=nthreads' $
              + " ; /all"
      printf, Slun, "a -> make_intdif_gains, pref = '" + prefilters[ipref] $
              + "', min=0.1, max=4.0, bad=1.0, smooth=3.0, timeaver=1L, /smallscale ; /all"
    endif
    printf, Slun, "a -> fitprefilter, pref = '"+prefilters[ipref]+"'" $
            + ", /mask ;, /hints, dir='10:02:45'"
    case instrument of
      'CRISP'  : printf, Slun, "a -> prepmomfbd" $
                         + ", /fill_fov" $
                         + ", /wb_states" $
                         + ", date_obs = '" + isodate + "'" $
                         + ", numpoints = 116" $
                         + ", global_keywords=['FIT_PLANE']" $
                         + ", pref = '"+prefilters[ipref]+"'" $
                         + ", margin = 5" $
                         + ", dirs=['"+strjoin(file_basename(dirarr), "','")+"'] "
      'CRISP2' : printf, Slun, "a -> prepmomfbd" $
                         + ", /unpol" $
                         + ", /fill_fov" $
                         + ", /wb_states" $
                         + ", date_obs = '" + isodate + "'" $
                         + ", numpoints = 116" $
                         + ", global_keywords=['FIT_PLANE']" $
                         + ", pref = '"+prefilters[ipref]+"'" $
                         + ", margin = 5" $
                         + ", dirs=['"+strjoin(file_basename(dirarr), "','")+"'] "
      else : stop
    endcase
  
    printf, Slun
  endfor                        ; ipref

  printf, Slun, ''
  printf, Slun, ';; Run MOMFBD outside IDL.'
  printf, Slun, ''

  printf, Slun, ';; Post-MOMFBD stuff:' 
  printf, Slun, "a -> make_scan_cube, 'momfbd/.../cfg/results/', scannos = '69', nthreads=nthreads, /circular_fov"
  printf, Slun, "a -> fitscube_wcs_improve_spatial, 'cubes_scan/nb....fits' ; If suitable target"
  printf, Slun, "; or "
  printf, Slun, "a -> make_wb_cube, 'momfbd/.../cfg/results/', /align_interactive, /circular_fov"
  printf, Slun, "a -> fitscube_wcs_improve_spatial, 'cubes_wb/wb....fits' ; If suitable target"
  printf, Slun, "a -> make_nb_cube, 'cubes_wb/wb....fits', nthreads=nthreads"

  
  free_lun, Clun
  free_lun, Slun

  ;; Do something about OBSERVER metadata keyword
  if ~keyword_set(no_observer_metadata) then begin
    ;; See if we can find some metadata by looking in the raw data dirs.
    data_dirs = root_dir + '/' + dirarr + '/'+instrument.capwords()+'-W/'
    ;; Pick the first file in each.
    data_files = file_search(data_dirs+'/*00000_0000000*fits', count = Nfiles)
    if Nfiles gt 0 then begin
      ;; Now look for OBSERVER keywords
      observers = strarr(Nfiles)
      for ifile = 0, Nfiles-1 do begin
        red_progressbar, ifile, Nfiles, 'Looking in '+instrument+' data for OBSERVER keyword'
        observers[ifile] = red_fitsgetkeyword(data_files[ifile], 'OBSERVER')
      endfor
      observers = ['', observers] ; empty default means no OBSERVER keyword in metadata
      indx = uniq(observers, sort(observers))
      print
      if n_elements(indx) gt 1 then begin
        print, inam + ' : Found OBSERVER keyword(s) in the '+instrument+' raw data. All is well.'
      endif else begin      
        print, inam + ' : Found no OBSERVER metadata in the '+instrument+' raw data.'
        observer = ''
        read, 'Add names for that keyword for the '+instrument+' workdir or hit return: ', observer
        ;; Write it to the metadata file
        print
        if observer ne '' then begin
          print, inam+' : Adding to '+instrument+' metadata, OBSERVER = '+observer
          red_metadata_store, fname = work_dir + '/info/metadata.fits' $
                              , [{keyword:'OBSERVER', value:observer $
                                  , comment:'Observer name(s)'}]
        endif else begin
          print, inam+' : No OBSERVER keyword in '+instrument+' metadata.'
        endelse
        print, inam+' : Edit '+work_dir + '/info/metadata.fits if you need to change this.'
      endelse
      print
    endif else begin
      print, inam+' : Did not find any CHROMIS data to get OBSERVER from.'
    endelse
  endif
  
;  ;; We will now attempt to copy existing sums of calibration data.
;  red_setupworkdir_copy, old_dir, 'darks',       work_dir
;  red_setupworkdir_copy, old_dir, 'flats',       work_dir
;  red_setupworkdir_copy, old_dir, 'pinhs',       work_dir
;  red_setupworkdir_copy, old_dir, 'polcal_sums', work_dir
 
end
