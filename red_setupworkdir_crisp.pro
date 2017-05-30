pro red_setupworkdir_crisp, work_dir, root_dir, cfgfile, scriptfile, isodate

  red_metadata_store, fname = work_dir + '/info/metadata.fits' $
                      , [{keyword:'INSTRUME', value:'CRISP' $
                          , comment:'Name of instrument'} $
                         , {keyword:'TELCONFG', value:'Schupmann, imaging', $
                            comment:'Telescope configuration'}]

  
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
  printf, Slun, 'a = crispred("'+cfgfile+'")' 
  printf, Slun, 'root_dir = "' + root_dir + '"'

  ;; Download SST log files and optionally some other data from the web.
  print, 'Log files'
  printf, Clun, '#'
  printf, Clun, '# --- Download SST log files'
  printf, Clun, '#'
  printf, Slun, 'a -> download ; add ", /all" to get also HMI images and AR maps.'

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
  
  darksubdirs = red_find_instrumentdirs(root_dir, 'crisp', 'dark' $
                                        , count = Nsubdirs)
  if Nsubdirs gt 0 then begin
    darkdirs = file_dirname(darksubdirs)
    darkdirs = darkdirs[uniq(darkdirs, sort(darkdirs))]
    for idir = 0, n_elements(darkdirs)-1 do begin
      printf, Clun, 'dark_dir = '+red_strreplace(darkdirs[idir] $
                                                 , root_dir, '')
      printf, Slun, 'a -> sumdark, /check, dirs=root_dir+"' $
              + red_strreplace(darkdirs[idir], root_dir, '') + '"'
    endfor                      ; idir
  endif                         ; Nsubdirs
  
  print, 'Flats'
  printf, Clun, '#'
  printf, Clun, '# --- Flats'
  printf, Clun, '#'

  flatsubdirs = red_find_instrumentdirs(root_dir, 'crisp', 'flat' $
                                        , count = Nsubdirs)
  if Nsubdirs gt 0 then begin
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
        ;; Print to script file
        printf, Slun, 'a -> sumflat, /check, dirs=root_dir+"' $
                + red_strreplace(flatdirs[idir], root_dir, '') $
                + '"  ; ' + camdirs+' ('+wavelengths+')'

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
  maybe_nodescatter = strarr(Nprefilters)
  indx7772 = where(prefilters eq '7772', N7772)
  if N7772 gt 0 then maybe_nodescatter[indx7772] = '; , /no_descatter'

  for ipref = 0, Nprefilters-1 do begin
    printf, Slun, "a -> makegains, pref='" + prefilters[ipref] $
            + "' " + maybe_nodescatter[ipref]
  endfor

  print, 'Pinholes'
  printf, Clun, '#'
  printf, Clun, '# --- Pinholes'
  printf, Clun, '#'
  pinhdirs = file_search(root_dir+'/pinh*/*', count = Ndirs, /fold)
  for i = 0, Ndirs-1 do begin
    pinhsubdirs = file_search(pinhdirs[i]+'/crisp*', count = Nsubdirs, /fold)
    if Nsubdirs gt 0 then begin
      printf, Clun, 'pinh_dir = '+red_strreplace(pinhdirs[i], root_dir, '')
      printf, Slun, 'a -> setpinhdir, root_dir+"' $
              + red_strreplace(pinhdirs[i], root_dir, '')+'"'
;        printf, Slun, 'a -> sumpinh_new'
      for ipref = 0, Nprefilters-1 do begin
        printf, Slun, "a -> sumpinh, /pinhole_align, pref='" $
                + prefilters[ipref]+"'" $
                + maybe_nodescatter[ipref]
      endfor
    endif else begin
      pinhsubdirs = file_search(pinhdirs[i]+'/*', count = Nsubdirs)
      for j = 0, Nsubdirs-1 do begin
        pinhsubsubdirs = file_search(pinhsubdirs[j]+'/crisp*' $
                                     , count = Nsubsubdirs, /fold)
        if Nsubsubdirs gt 0 then begin
          printf, Clun, 'pinh_dir = ' $
                  + red_strreplace(pinhsubdirs[j], root_dir, '')
          printf, Slun, 'a -> setpinhdir, root_dir+"' $
                  + red_strreplace(pinhsubdirs[j], root_dir, '')+'"'
;              printf, Slun, 'a -> sumpinh_new'
          for ipref = 0, Nprefilters-1 do begin
            printf, Slun, "a -> sumpinh, /pinhole_align, pref='" $
                    + prefilters[ipref] + "'" $
                    + maybe_nodescatter[ipref]
          endfor                ; ipref
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
        printf, Slun, 'a -> setpolcaldir, root_dir+"' $
                + red_strreplace(polcaldirs[i], root_dir, '')+'"'
        printf, Slun, 'a -> sumpolcal, /check'
      endif else begin
        polcalsubdirs = file_search(polcaldirs[i]+'/*', count = Nsubdirs)
        for j = 0, Nsubdirs-1 do begin
          polcalsubsubdirs = file_search(polcalsubdirs[j]+'/crisp*' $
                                         , count = Nsubsubdirs, /fold)
          if Nsubsubdirs gt 0 then begin
            printf, Clun, 'polcal_dir = ' $
                    + red_strreplace(polcalsubdirs[j], root_dir, '')
;              Npol += 1
            printf, Slun, 'a -> setpolcaldir, root_dir+"' $
                    + red_strreplace(polcalsubdirs[j], root_dir, '')+'"'
            printf, Slun, 'a -> sumpolcal, /check' 
          endif
        endfor                  ; j
      endelse
    endfor                      ; i

    for ipref = 0, Npol-1 do begin
      printf, Slun, "a -> polcalcube, pref='"+polprefs[ipref]+"' " $
              + maybe_nodescatter[ipref]
      printf, Slun, "a -> polcal, pref='"+polprefs[ipref]+"', nthreads=" $
              + strtrim(Nthreads, 2)
    endfor                      ; ipref
    
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
    if total(prefilters[ipref] eq polprefs) gt 0 then begin
      printf, Slun, "a -> prepflatcubes, pref='"+prefilters[ipref]+"'" $
              + maybe_nodescatter[ipref]
    endif else begin
      printf, Slun, "a -> prepflatcubes_lc4, pref='"+prefilters[ipref]+"'" $
              + maybe_nodescatter[ipref]
    endelse
  endfor                        ; ipref

  
  printf, Slun, ''
  printf, Slun, 'a -> getalignclips_new' 
  printf, Slun, 'a -> getoffsets' 
  
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