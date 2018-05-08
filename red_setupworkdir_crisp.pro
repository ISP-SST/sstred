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
; 
; 
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
; 
; 
;-
pro red_setupworkdir_crisp, work_dir, root_dir, cfgfile, scriptfile, isodate $
                            , calibrations_only = calibrations_only

  red_metadata_store, fname = work_dir + '/info/metadata.fits' $
                      , [{keyword:'INSTRUME', value:'CRISP' $
                          , comment:'Name of instrument'} $
                         , {keyword:'TELCONFG', value:'Schupmann, imaging', $
                            comment:'Telescope configuration'}]

  
  ;; Are there darks and flats?
  darksubdirs = red_find_instrumentdirs(root_dir, 'crisp', 'dark' $
                                        , count = Ndarkdirs)
  flatsubdirs = red_find_instrumentdirs(root_dir, 'crisp', 'flat' $
                                        , count = Nflatdirs)
  
  if Ndarkdirs eq 0 then begin
    print, 'No CRISP darks were found. No setup generated.'
    return
  endif

  ;; Open two files for writing. Use logical unit Clun for a Config
  ;; file and Slun for a Script file.

  openw, Clun, work_dir + cfgfile, /get_lun
  openw, Slun, work_dir + scriptfile , /get_lun

  ;; Specify the date in the config file, ISO format.
  print, 'Date'
  printf, Clun, '#'
  printf, Clun, '# --- Date'
  printf, Clun, '#'
  printf, Clun,'isodate = '+isodate

  ;; printf, Slun, '.r crispred'
  if keyword_set(calibrations_only) then begin
    printf, Slun, 'a = crispred("'+cfgfile+'",/dev)' 
  endif else begin
    printf, Slun, 'a = crispred("'+cfgfile+'")' 
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
    printf, Slun, 'a -> download ; add ", /all" to get also HMI images and AR maps.'
  endif
  
  print, 'Cameras'
  printf, Clun, '#'
  printf, Clun, '# --- Cameras'
  printf, Clun, '#'
  printf, Clun, 'camera = Crisp-T'
  printf, Clun, 'camera = Crisp-R'
  printf, Clun, 'camera = Crisp-W'
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
    ;; There are CRISP darks!

    ;; Directories with camera dirs below:    
    darkdirs = file_dirname(darksubdirs)
    darkdirs = darkdirs[uniq(darkdirs, sort(darkdirs))]
    Ndarkdirs = n_elements(darkdirs)

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

      ;; Print to script file
      printf, Slun, 'a -> sumdark, /sum_in_rdx, /check, dirs=root_dir+"' $
              + red_strreplace(darkdirs[idir], root_dir, '') + '"' $
              + ', nthreads=nthreads' $
              + outdirkey 
      
    endfor                      ; idir
  endif                         ; Nsubdirs
  

  if Nflatdirs eq 0 then begin
    print, 'No CRISP flats were found. Stop after summing darks.'
    free_lun, Clun
    free_lun, Slun
    return
  endif
 
  
  print, 'Flats'
  printf, Clun, '#'
  printf, Clun, '# --- Flats'
  printf, Clun, '#'

  if Nflatdirs gt 0 then begin
    ;; There are CRISP flats!

    ;; Directories with camera dirs below:
    flatdirs = file_dirname(flatsubdirs)
    flatdirs = flatdirs[uniq(flatdirs, sort(flatdirs))]
    Nflatdirs = n_elements(flatdirs)

    ;; Loop over the flatdirs, write each to the config file and
    ;; to the script file
    for idir = 0, Nflatdirs-1 do begin

      ;; Config file

      printf, Clun, 'flat_dir = '+red_strreplace(flatdirs[idir] $
                                                 , root_dir, '')

      ;; Script file
      
      ;; Look for wavelengths in those flatsubdirs that match
      ;; flatdirs[idir]! Also collect prefilters.
      indx = where(strmatch(flatsubdirs, flatdirs[idir]+'*'))
      fnames = file_search(flatsubdirs[indx]+'/cam*', count = Nfiles)
      if Nfiles gt 0 then begin
        
        camdirs = strjoin(file_basename(flatsubdirs[indx]), ' ')

        red_extractstates, fnames, /basename, pref = wls
        wls = wls[uniq(wls, sort(wls))]
        wls = wls[WHERE(wls ne '')]
        wavelengths = strjoin(wls, ' ')

        
        if keyword_set(calibrations_only) then begin
          ;; For /calibrations_only we want to output the summed data
          ;; in timestamp directories so we can handle multiple sets.
          outdir = 'flats/' + file_basename(flatdirs[idir])
          outdirkey = ', outdir="'+outdir+'", /softlink, /store_rawsum'
        endif else outdirkey = ''

        ;; Print to script file
        printf, Slun, 'a -> sumflat, /sum_in_rdx, /check, dirs=root_dir+"' $
                + red_strreplace(flatdirs[idir], root_dir, '')+ '"' $
                + ', nthreads=nthreads' $
                + outdirkey $
                + '  ; ' + camdirs+' ('+wavelengths+')'

        red_append, prefilters, wls

      endif                     ; Nfiles
    endfor                      ; idir
  endif                         ; Nsubdirs

  
  if n_elements(prefilters) gt 0 then $
     prefilters = prefilters[uniq(prefilters, sort(prefilters))]
  Nprefilters = n_elements(prefilters)

  if Nprefilters eq 0 then begin
    ;; This can happen if flats were already summed in La Palma. Look
    ;; for prefilters in the summed flats directory instead.
    spawn, 'ls flats/cam*.flat | cut -d. -f2|sort|uniq', prefilters
    Nprefilters = n_elements(prefilters)
  endif
  
  ;; For the 7772 Å prefilter a /no_descatter keyword may be needed in
  ;; some of the method calls, so add it commented out. (This if
  ;; because for some years we don't have properly prepared
  ;; backgains and psfs for the relevant cameras.)
  if ~keyword_set(calibrations_only) then begin  
    maybe_nodescatter = strarr(Nprefilters)
    indx7772 = where(prefilters eq '7772', N7772)
    if N7772 gt 0 then maybe_nodescatter[indx7772] = '; , /no_descatter'
  endif
  
;  if ~keyword_set(calibrations_only) then begin  
;    for ipref = 0, Nprefilters-1 do begin
;      printf, Slun, "a -> makegains, pref='" + prefilters[ipref] $
;              + "' " + maybe_nodescatter[ipref]
;    endfor
;  endif

  print, 'Pinholes'
  printf, Clun, '#'
  printf, Clun, '# --- Pinholes'
  printf, Clun, '#'
  pinhdirs = file_search(root_dir+'/pinh*/*', count = Ndirs, /fold)
  for i = 0, Ndirs-1 do begin
    pinhsubdirs = file_search(pinhdirs[i]+'/crisp*', count = Nsubdirs, /fold)
    if Nsubdirs gt 0 then begin
      printf, Clun, 'pinh_dir = '+red_strreplace(pinhdirs[i], root_dir, '')
      
      if keyword_set(calibrations_only) then begin
        ;; For /calibrations_only we want to output the summed data in
        ;; timestamp directories so we can handle multiple sets.
        outdir = 'pinhs/' + file_basename(pinhdirs[i])
        outdirkey = ', outdir="'+outdir+'"'
      endif else outdirkey = ''

      printf, Slun, 'a -> sumpinh, /sum_in_rdx, /pinhole_align, dirs=root_dir+"' $
              + red_strreplace(pinhdirs[i], root_dir, '') + '"' $
              + ', nthreads=nthreads' $
              + outdirkey 
    endif else begin
      pinhsubdirs = file_search(pinhdirs[i]+'/*', count = Nsubdirs)
      for j = 0, Nsubdirs-1 do begin
        pinhsubsubdirs = file_search(pinhsubdirs[j]+'/crisp*' $
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

          printf, Slun, 'a -> sumpinh, /sum_in_rdx, /pinhole_align, dirs=root_dir+"' $
                  + red_strreplace(pinhsubsubdirs[j], root_dir, '')  + '"' $
                  + ', nthreads=nthreads' $
                  + outdirkey 
        endif                   ; Nsubsubdirs
      endfor                    ; j
    endelse                     ; Nsubdirs
  endfor                        ; i


  print, 'Polcal'
  printf, Clun, '#'
  printf, Clun, '# --- Polcal'
  printf, Clun, '#'
;  Npol = 0
  polcaldirs = file_search(root_dir+'/polc*/*', count = Npol, /fold)
  if Npol gt 0 then begin
    polprefs = file_basename(polcaldirs)
    for i = 0, Npol-1 do begin
      polcalsubdirs = file_search(polcaldirs[i]+'/crisp*' $
                                  , count = Nsubdirs, /fold)
      if Nsubdirs gt 0 then begin
        printf, Clun, 'polcal_dir = ' $
                + red_strreplace(polcaldirs[i], root_dir, '')
;        Npol += 1
;          printf, Slun, 'a -> setpolcaldir, root_dir+"' $
;                  + red_strreplace(polcaldirs[i], root_dir, '')+'"'
        
        if keyword_set(calibrations_only) then begin
          ;; For /calibrations_only we want to output the summed data in
          ;; timestamp directories so we can handle multiple sets.
          outdir = 'polcal_sums/' + file_basename(polcaldirs[i])
          outdirkey = ', outdir="'+outdir+'"'
        endif else outdirkey = ''

        printf, Slun, 'a -> sumpolcal, /sum_in_rdx, /check, dirs=root_dir+"' $
                + red_strreplace(polcaldirs[i], root_dir, '')+'"' $
                + ', nthreads=nthreads' $
                + outdirkey 
      endif else begin
        polcalsubdirs = file_search(polcaldirs[i]+'/*', count = Nsubdirs)
        for j = 0, Nsubdirs-1 do begin
          polcalsubsubdirs = file_search(polcalsubdirs[j]+'/crisp*' $
                                         , count = Nsubsubdirs, /fold)
          if Nsubsubdirs gt 0 then begin
            printf, Clun, 'polcal_dir = ' $
                    + red_strreplace(polcalsubdirs[j], root_dir, '')
;              Npol += 1
;              printf, Slun, 'a -> setpolcaldir, root_dir+"' $
;                      + red_strreplace(polcalsubdirs[j], root_dir, '')+'"'
            
            if keyword_set(calibrations_only) then begin
              ;; For /calibrations_only we want to output the summed data in
              ;; timestamp directories so we can handle multiple sets.
              outdir = 'polcal_sums/' + file_basename(polcalsubdirs[j])
              outdirkey = ', outdir="'+outdir+'"'
            endif else outdirkey = ''
            
            printf, Slun, 'a -> sumpolcal, /sum_in_rdx, /check, dirs=root_dir+"' $
                    + red_strreplace(polcalsubdirs[j], root_dir, '')+'"' $
                    + ', nthreads=nthreads' $
                    + outdirkey 
          endif
        endfor                  ; j
      endelse
    endfor                      ; i

    if ~keyword_set(calibrations_only) then begin  
      for ipref = 0, Npol-1 do begin
        printf, Slun, "a -> polcalcube, pref='" + polprefs[ipref] + "' " $
                + ", nthreads=nthreads" $
                + maybe_nodescatter[ipref] 
        printf, Slun, "a -> polcal, pref='" + polprefs[ipref] + "' " $
                + ", nthreads=nthreads"
      endfor                    ; ipref
    endif
      
  endif else begin
    polprefs = ''
  endelse                       ; Npol
  
  
  print, 'Prefilter scan'
  printf, Clun, '#'
  printf, Clun, '# --- Prefilter scan'
  printf, Clun, '#'
  Npfs = 0
  pfscandirs = file_search(root_dir+'/pfscan*/*', count = Ndirs, /fold)
  for i = 0, Ndirs-1 do begin
    pfscansubdirs = file_search(pfscandirs[i]+'/crisp*', count = Nsubdirs, /fold)
    if Nsubdirs gt 0 then begin
      printf, Clun, '# pfscan_dir = '+red_strreplace(pfscandirs[i], root_dir, '')
      Npfs += 1
    endif else begin
      pfscansubdirs = file_search(pfscandirs[i]+'/*', count = Nsubdirs)
      for j = 0, Nsubdirs-1 do begin
        pfscansubsubdirs = file_search(pfscansubdirs[j]+'/crisp*' $
                                       , count = Nsubsubdirs, /fold)
        if Nsubsubdirs gt 0 then begin
          printf, Clun, '# pfscan_dir = ' $
                  + red_strreplace(pfscansubdirs[j], root_dir, '')
          Npfs += 1
        endif
      endfor
    endelse
  endfor
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

  ;;  sciencedirs = file_search(root_dir+'/sci*/*', count = Ndirs, /fold)
  red_append, nonsciencedirs, pinhdirs
  red_append, nonsciencedirs, polcaldirs
  red_append, nonsciencedirs, pfscandirs
  red_append, nonsciencedirs, darkdirs
  red_append, nonsciencedirs, flatdirs
  sciencedirs = file_search(root_dir+'/*/*', count = Ndirs)

  for i = 0, Ndirs-1 do begin

    if total(sciencedirs[i] eq nonsciencedirs) eq 0 then begin
      sciencesubdirs = file_search(sciencedirs[i]+'/crisp*' $
                                   , count = Nsubdirs, /fold)
      if Nsubdirs gt 0 then begin
        red_append, dirarr, red_strreplace(sciencedirs[i], root_dir, '')
      endif else begin
        sciencesubdirs = file_search(sciencedirs[i]+'/*', count = Nsubdirs)
        for j = 0, Nsubdirs-1 do begin
          sciencesubsubdirs = file_search(sciencesubdirs[j]+'/crisp*' $
                                          , count = Nsubsubdirs, /fold)
          if Nsubsubdirs gt 0 then begin
            red_append, dirarr, red_strreplace(sciencesubdirs[j], root_dir, '')
          endif                 ; Nsubsubdirs
        endfor                  ; j
      endelse                   ; Nsubdirs
    endif
  endfor                        ; i
  if n_elements(dirarr) gt 0 then printf, Clun, "data_dir = ['"+strjoin(dirarr, "','")+"']"

  printf, Slun, 'a -> link_data' 
  
  for ipref = 0, Nprefilters-1 do begin
;    if total(prefilters[ipref] eq polprefs) gt 0 then begin
      printf, Slun, "a -> prepflatcubes, pref='"+prefilters[ipref]+"'" $
              + maybe_nodescatter[ipref]
;    endif else begin
;      printf, Slun, "a -> prepflatcubes_lc4, pref='"+prefilters[ipref]+"'" $
;              + maybe_nodescatter[ipref]
;    endelse
  endfor                        ; ipref

  
;  printf, Slun, ''
;  printf, Slun, 'a -> getalignclips_new' 
;  printf, Slun, 'a -> getoffsets' 
  
  printf, Slun, ''
  printf, Slun, 'a -> pinholecalib, nref=10'
  
  printf, Slun, ''
  printf, Slun, ';; -----------------------------------------------------'
  printf, Slun, ';; This is how far we should be able to run unsupervised'
  printf, Slun, 'stop'          
  printf, Slun, ''

  printf, Slun, '; The fitgains step requires the user to look at the fit and determine'
  printf, Slun, '; whether npar=3 or npar=4 is needed.'
  printf, Slun, 'a -> fitgains, npar = 2, res=res' 
  printf, Slun, '; If you need per-pixel reflectivities for your analysis'
  printf, Slun, '; (e.g. for atmospheric inversions) you can set the /fit_reflectivity'
  printf, Slun, '; keyword:'
  printf, Slun, '; a -> fitgains, npar = 3, res=res, /fit_reflectivity  '
  printf, Slun, '; However, running without /fit_reflectivity is safer. In should not'
  printf, Slun, '; be used for chromospheric lines like 6563 and 8542.'
  printf, Slun, ''
  printf, Slun, "a -> makegains, smooth=3.0, min=0.1, max=4.0, bad=1.0"
  printf, Slun, ''

  printf, Slun, '; If MOMFBD has problems near the edges, try to increase the margin in the call the prepmomfbd.'
  for ipref = 0, Nprefilters-1 do begin
    printf, Slun, "a -> sum_data_intdif, pref = '" + prefilters[ipref] $
            + "', cam = 'Crisp-T', /verbose, /show, /overwrite " + maybe_nodescatter[ipref] + " ; /all"
    printf, Slun, "a -> sum_data_intdif, pref = '" + prefilters[ipref] $
            + "', cam = 'Crisp-R', /verbose, /show, /overwrite " + maybe_nodescatter[ipref] + " ; /all"
    printf, Slun, "a -> make_intdif_gains3, pref = '" + prefilters[ipref] $
            + "', min=0.1, max=4.0, bad=1.0, smooth=3.0, timeaver=1L, /smallscale ; /all"
    if strmid(prefilters[ipref], 0, 2) eq '63' then begin
      printf, Slun, "a -> fitprefilter, fixcav = 2.0d, pref = '"+prefilters[ipref]+"', shift=-0.5"
    endif else begin
      printf, Slun, "a -> fitprefilter, fixcav = 2.0d, pref = '"+prefilters[ipref]+"'"
    endelse
    printf, Slun, "a -> prepmomfbd, /wb_states, date_obs = '" + isodate $
            + "', numpoints = 88, pref = '"+prefilters[ipref]+"', margin = 5 " $
            + maybe_nodescatter[ipref]
  endfor


  printf, Slun, ''
  printf, Slun, ';; Run MOMFBD outside IDL.'
  printf, Slun, ''

  printf, Slun, ';; Post-MOMFBD stuff:' 
  printf, Slun, 'a -> make_unpol_crispex, /noflat [, /scans_only,/wbwrite]        ; For unpolarized data'
  if Npol gt 0 then begin
    printf, Slun, 'pol = a->polarim(/new)' 
    printf, Slun, 'for i = 0, n_elements(pol)-1, 1 do pol[i]->demodulate,/noflat' 
    printf, Slun, 'a -> make_pol_crispex [, /scans_only,/wbwrite]          ; For polarized data'
  endif
  printf, Slun, 'a -> polish_tseries, np = 3 [, /negangle, xbd =, ybd =, tstep = ...]'
  
  free_lun, Clun
  free_lun, Slun


end
