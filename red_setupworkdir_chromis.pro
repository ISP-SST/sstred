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
; 
; :Returns:
; 
; 
; :Params:
; 
; 
; :Keywords:
;
;    calibrations_only : in, optional, type=boolean
;
;      Set up to process calibration data only.
; 
;   
;   
;   
; 
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
; 
;-
pro red_setupworkdir_chromis, work_dir, root_dir, cfgfile, scriptfile, isodate $
                              , calibrations_only = calibrations_only

  red_metadata_store, fname = work_dir + '/info/metadata.fits' $
                      , [{keyword:'INSTRUME', value:'CHROMIS' $
                          , comment:'Name of instrument'} $
                         , {keyword:'TELCONFG', value:'Schupmann, imaging table', $
                            comment:'Telescope configuration'}]


  ;; Are there darks and flats?
  darksubdirs = red_find_instrumentdirs(root_dir, 'chromis', '*dark*' $
                                        , count = Ndarkdirs)
  flatsubdirs = red_find_instrumentdirs(root_dir, 'chromis', '*flat*' $
                                        , count = Nflatdirs)

  if Ndarkdirs eq 0 then begin
    print, 'No CHROMIS darks were found. No setup generated.'
    return
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
  printf, Clun,'image_scale = 0.0379'   ; Measured in May 2016.
  printf, Clun, 'diversity = 3.35e-3'   ; Nominal value for 2016.

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
  if ~keyword_set(calibrations_only) then begin
    print, 'Log files'
    printf, Clun, '#'
    printf, Clun, '# --- Download SST log files'
    printf, Clun, '#'
    printf, Slun
    printf, Slun, 'a -> download ; add ", /all" to get also HMI images and AR maps.'
  endif

  
;  if Nflatdirs gt 0 then begin
;    ;; This step requires flats. And the hrz --> tuning conversion is
;    ;; not needed if we don't have them.
;    printf, Slun
;    printf, Slun, 'a -> hrz_zeropoint' ; Find the reference wavelenght of CHROMIS scans
;  endif
  
  ;; Analyze directories and produce r0 plots (optional)
  if ~keyword_set(calibrations_only) then begin
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
  printf, Clun, 'camera = Chromis-W'
  printf, Clun, 'camera = Chromis-D'
  printf, Clun, 'camera = Chromis-N'
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

    printf, Slun
    ;; Loop over the darkdirs, write each to the config file and
    ;; to the script file
    for idir = 0, Ndarkdirs-1 do begin

      ;; Config file

      printf, Clun, 'dark_dir = '+red_strreplace(darkdirs[idir] $
                                                 , root_dir, '')

      ;; Script file
      
      if keyword_set(calibrations_only) then begin
        ;; For /calibrations_only we want to output the summed data in
        ;; timestamp directories so we can handle multiple sets.
        outdir = 'darks/' + file_basename(darkdirs[idir])
        outdirkey = ', outdir="'+outdir+'", /softlink'
      endif else outdirkey = ''
      
      printf, Slun, 'a -> sumdark, /sum_in_rdx, /check, dirs=root_dir+"' $
              + red_strreplace(darkdirs[idir], root_dir, '') + '"' $
              + ', nthreads=nthreads' $
              + outdirkey 
    endfor                      ; idir
  endif                         ; Nsubdirs


  if Nflatdirs eq 0 then begin
    print, 'No CHROMIS flats were found. Stop after summing darks.'
    free_lun, Clun
    free_lun, Slun
    return
  endif
 
  
  print, 'Flats'
  printf, Clun, '#'
  printf, Clun, '# --- Flats'
  printf, Clun, '#'

  if Nflatdirs gt 0 then begin
    ;; There are CHROMIS flats!

    ;; Directories with camera dirs below:
    flatdirs = file_dirname(flatsubdirs)
    flatdirs = flatdirs[uniq(flatdirs, sort(flatdirs))]
    Nflatdirs = n_elements(flatdirs)

    ;; Do the linedefs and hrz_calib thing
    red_download_linedefs, isodate, flatdirs, work_dir
    chromis_hrz_zeropoint, work_dir

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
        red_extractstates, fnames, /basename, pref = wls, is_wb = this_is_wb, is_pd = this_is_pd
        indx = uniq(wls, sort(wls))
        wls = wls[indx]
        this_is_wb = this_is_wb[indx]
        this_is_pd = this_is_pd[indx]
        wls = wls[WHERE(wls ne '')]

        if n_elements(wls) gt 0 then begin
          wavelengths = strjoin(wls, ' ')
          
          if keyword_set(calibrations_only) then begin
            ;; For /calibrations_only we want to output the summed data
            ;; in timestamp directories so we can handle multiple sets.
            outdir = 'flats/' + file_basename(flatdirs[idir])
            outdirkey = ', outdir="'+outdir+'", /softlink, /store_rawsum'
          endif else outdirkey = ''
        
          ;; Print to script file
          printf, Slun, 'a -> sumflat, /sum_in_rdx, /check' $
                  + ', dirs=root_dir+"' + red_strreplace(flatdirs[idir], root_dir, '') + '"' $
                  + ', nthreads=nthreads' $
                  + outdirkey $
                  + ' ; ' + camdirs+' ('+wavelengths+')'

          red_append, prefilters, wls
          red_append, is_wb, this_is_wb
          red_append, is_pd, this_is_pd
        endif
      endif                     ; Nfiles
    endfor                      ; idir
  endif                         ; Nsubdirs
  
  Nprefilters = n_elements(prefilters)
  if Nprefilters gt 0 then begin
    indx = uniq(prefilters, sort(prefilters))
    is_wb = is_wb[indx]
    is_pd = is_pd[indx]
    prefilters = prefilters[indx]
  endif else begin
    ;; This can happen if flats were already summed in La Palma and
    ;; then deleted. Look for prefilters in the summed flats directory
    ;; instead.
    spawn, 'ls flats/cam*.flat | cut -d. -f2|sort|uniq', prefilters
    ;; Need to set also cameras here.
  endelse
  Nprefilters = n_elements(prefilters)

  if ~keyword_set(calibrations_only) then begin  
    for ipref = 0, Nprefilters-1 do begin
      if ~is_wb[ipref] then printf, Slun, "a -> prepflatcubes, pref='"+prefilters[ipref]+"'"
    endfor                      ; ipref
  endif
  
  if ~keyword_set(calibrations_only) then begin  
    printf, Slun, ''
    printf, Slun, '; The fitgains step requires the user to look at the fit and determine'
    printf, Slun, '; whether npar=3 or npar=4 is needed.'
    for ipref = 0, Nprefilters-1 do begin
      if ~is_wb[ipref] then begin
        printf, Slun, "a -> fitgains, rebin=800L, Niter=3L, Nthreads=12L, Npar=5L, res=res, pref='" $
                + prefilters[ipref] + "'"
      end
    endfor                      ; ipref
    printf, Slun, '; If you need per-pixel reflectivities for your analysis'
    printf, Slun, '; (e.g. for atmospheric inversions) you can set the /fit_reflectivity'
    printf, Slun, '; keyword:'
    printf, Slun, '; a -> fitgains, npar = 3, res=res, /fit_reflectivity  '
    printf, Slun, '; However, running without /fit_reflectivity is safer. In should not'
    printf, Slun, '; be used for chromospheric lines like 6563 and 8542.'
    printf, Slun, '; Sometimes you need to add a few spline nodes in order to make the fits work,'
    printf, Slun, '; particularly just to the red and to the blue of the densely sampled region and'
    printf, Slun, '; also in blends if they are not propely sampled.'
    printf, Slun, '; As an example, to do this for 3969 Å, 12.00ms_G10.00 data, do something like'
    printf, Slun, '; the following:'
    printf, Slun, "; restore,'flats/spectral_flats/camXXX_12.00ms_G10.00_3969_flats.sav'; Read flats cube"
    printf, Slun, '; myg = [wav*1.d10, -0.765d0, 0.740d0] ; Add wavelength points'
    printf, Slun, '; myg=myg[sort(myg)]  ; Sort the wavelength points'
    printf, Slun, '; a->fitgains, rebin=800L, niter=3L, nthreads=12L, res=res, npar=5L, myg=myg ; Run the fitgain step with the added wavelength points'
    printf, Slun, '; Then, if you have already run makegains, rerun it.'

    printf, Slun, ''

    printf, Slun, "a -> makegains, smooth=3.0, min=0.1, max=4.0, bad=1.0"

  endif
    

  printf, Slun, ''
  print, 'Pinholes'
  printf, Clun, '#'
  printf, Clun, '# --- Pinholes'
  printf, Clun, '#'
  pinhdirs = file_search(root_dir+'/*pinh*/*', count = Ndirs, /fold)
  for i = 0, Ndirs-1 do begin
    pinhsubdirs = file_search(pinhdirs[i]+'/chromis*', count = Nsubdirs, /fold)
    if Nsubdirs gt 0 then begin
      printf, Slun
      printf, Clun, 'pinh_dir = '+red_strreplace(pinhdirs[i], root_dir, '')
      
      if keyword_set(calibrations_only) then begin
        ;; For /calibrations_only we want to output the summed data in
        ;; timestamp directories so we can handle multiple sets.
        outdir = 'pinhs/' + file_basename(pinhdirs[i])
        outdirkey = ', outdir="'+outdir+'"'
      endif else outdirkey = ''

      for ipref = 0, Nprefilters-1 do begin
        printf, Slun, "a -> sumpinh, /sum_in_rdx, /pinhole_align" $
                + ', nthreads=nthreads' $
                + outdirkey $
                + ", dirs=root_dir+'" + red_strreplace(pinhdirs[i], root_dir, '') $
                + "', pref='"+prefilters[ipref]+"'"
      endfor                    ; ipref
    endif else begin
      pinhsubdirs = file_search(pinhdirs[i]+'/*', count = Nsubdirs)
      printf, Slun
      for j = 0, Nsubdirs-1 do begin
        pinhsubsubdirs = file_search(pinhsubdirs[j]+'/chromis*' $
                                     , count = Nsubsubdirs, /fold)
        if Nsubsubdirs gt 0 then begin
          printf, Clun, 'pinh_dir = ' $
                  + red_strreplace(pinhsubdirs[j], root_dir, '')
          
          if keyword_set(calibrations_only) then begin
            ;; For /calibrations_only we want to output the summed data in
            ;; timestamp directories so we can handle multiple sets.
            outdir = 'pinhs/' + file_basename(pinhsubdirs[j])
            outdirkey = ', outdir="'+outdir+'"'
          endif else outdirkey = ''

          for ipref = 0, Nprefilters-1 do begin
            printf, Slun, "a -> sumpinh, /sum_in_rdx, /pinhole_align" $
                    + ', nthreads=nthreads' $
                    + outdirkey $
                    + ", dirs=root_dir+'" +  red_strreplace(pinhsubdirs[j], root_dir, '') $
                    + "';, pref='" + prefilters[ipref]+"'" 
          endfor                ; ipref
        endif
      endfor                    ; j
    endelse
  endfor                        ; i

  if ~keyword_set(calibrations_only) then begin  
    printf, Slun, ''
    printf, Slun, 'a -> pinholecalib, nref=10'
;    printf, Slun, 'a -> diversitycalib'
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

  ;; Exclude directories know not to have science data
  red_append, nonsciencedirs, pinhdirs
  red_append, nonsciencedirs, pfscandirs
  red_append, nonsciencedirs, darkdirs
  red_append, nonsciencedirs, flatdirs
  ;; Sometimes the morning calibrations data were not deleted, so they
  ;; have to be excluded too.
  calibsubdirs = red_find_instrumentdirs(root_dir, 'chromis', '*calib' $
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

  printf, Slun, ''
  printf, Slun, 'a -> link_data' 
  
  printf, Slun, ''
  printf, Slun, ';; -----------------------------------------------------'
  printf, Slun, ';; This is how far we should be able to run unsupervised'
  printf, Slun, 'stop'          
  printf, Slun, ''

  ;;  for ipref = 0, Nprefilters-1 do printf, Slun, "a -> fitprefilter, /mask ;, pref = '"+prefilters[ipref]+"'"
  printf, Slun, "a -> fitprefilter, /mask"

  printf, Slun, ''
  for ipref = 0, Nprefilters-1 do begin
    if is_wb[ipref] then begin
      printf, Slun, "a -> prepmomfbd" $
;              + ", date_obs='" + isodate + "'" $
              + ", Nremove=2" $
              + ", Nmodes=60" $
              + ", numpoints=128" $
;            + ", margin=5 " $
              + ", global_keywords=['FIT_PLANE']" $
              + ", maxshift=45" $
              + ", /wb_states" $
;              + ", /redux" $
              + ", /unpol" $
              + ", extraclip = [75,125,15,15]" $
              + ", pref='" + prefilters[ipref] + "'" $
              + ";, dirs=['09:03:24','09:24:21']" 
    endif
  endfor                        ; ipref ;

  printf, Slun, ''
  printf, Slun, ';; Run MOMFBD outside IDL.'
  printf, Slun, ''

  printf, Slun, ';; Post-MOMFBD stuff:' 
  printf, Slun, "a -> align_continuum"
  printf, Slun
  printf, Slun, "a -> make_wb_cube, 'momfbd/.../cfg/results/', /interactive, /autocrop"
  printf, Slun, "a -> make_nb_cube, 'cubes_wb/wb....fits'"
  printf, Slun, "; or "
  printf, Slun, "a -> make_scan_cube, 'momfbd/.../cfg/results/', /autocrop, scannos = '69'"

;  
;  printf, Slun, "a->polish_tseries" $
;;            + ", /full" $
;;            + ", /fitsoutput" $
;          + ", xbd=1280, ybd=1024" $
;          + ", np=5"
;
;  printf, Slun, "a->make_crispex, /float, /aligncont"
;; a -> make_crispex, /noflat, /scans_only, /float, /aligncont, /wbwrite

  
  free_lun, Clun
  free_lun, Slun

  
end
