; docformat = 'rst'

;+
; Makes a pipeline config file.
; 
; :Categories:
;
;    SST observations
; 
; 
; :Author:
; 
;    Mats Löfdahl, 2013-07-08
; 
; 
; :Params:
; 
;    root_dir : in, string
; 
;      The top directory of your saved data (or a regular expression
;      that matches it). If this directory name does not contain a
;      date, an attempt will be made to get the date from out_dir.
; 
; :Keywords:
; 
;    cfgfile : in, optional, type=string, default='config.txt'
; 
;      The name of the generated config file.
; 
;    download_all : in, optional, type=boolean
;
;      Set this to download auxiliary data, like SDO/HMI images and AR
;      maps. Otherwise download SST log file only.
;
;    scriptfile : in, optional, type=string, default='doit.pro'
; 
;      The name of the generated script file. The script file can be
;      run in an idl session with "@doit.pro" (assuming the default
;      name). It will perform the basic things, like co-adding of
;      darks, flats, etc. Later commands, that involve human
;      interaction are present in the file but commented out.
; 
;    date : in, optional, type=string
; 
;      The date (in iso format) the data was collected.
; 
;    out_dir : in, optional, type=string, default='current directory'
; 
;      The output directory to be used by crispred.
; 
;    lapalma : in, optional, type=boolean
; 
;       If this is set, will search for root_dir's date in
;       "/data/disk?/*/" and "/data/camera?/*/" (where data is usually
;       found in La Palma)'
; 
;    stockholm : in, optional, type=boolean
; 
;       If this is set, will search for out_dir's date in
;       "/mnt/sand??" (where data is usually found in Stockholm). If
;       no other information is given about where to look for data,
;       stockholm is assumed.
; 
; 
; :History:
; 
;    2013-07-10 : MGL. Will now get a date from out_dir if none is
;                 present in root_dir.
; 
;    2013-08-29 : MGL. Any subdirectory that is not a known
;                 calibration data directory (darks, flats, etc.) is a
;                 potential science data directory.
; 
;    2013-08-30 : MGL. Take care of prefilter scans.
; 
;    2013-11-25 : MGL. Renamed from original name sst_makeconfig and
;                 moved to Crispred repository. You can now use
;                 regular expressions for root_dir. New keywords
;                 "date", "lapalma" and "stockholm". Made root_dir a
;                 keyword. 
; 
;    2013-11-26 : MGL. Various improvements. The "new" flat field
;                 procedure. Find out which cameras and wavelengths
;                 are present in the raw flats directories.
; 
;    2013-11-29 : MGL. Changed the order of some commands in the
;                 doit.pro file. Add "/descatter" keyword to
;                 sum_data_intdif call only for wavelengths > 7700. 
; 
;    2013-12-02 : MGL. Bugfixes. Add slash at end of root_dir. Can now
;                 deal with empty raw flats directories. Gets the
;                 prefilters from the summed flats instead.
;
;    2013-12-02 : MGL. Move prepflatcubes[_lc4] to the unsupervised
;                 part. Default parameters for fitgains[_ng]. Deal
;                 with data sets with more than a single polcal
;                 directory and prefilter. Set number of threads to
;                 use based on the number of CPUs.
;
;    2013-12-09 : Pit. Allow root_dir to be the actual date directory.
;
;    2013-12-19 : MGL. Download SST logfiles and some other data from
;                 the web. 
;
;    2013-12-20 : MGL. Change calls from fitgains or fitgains_ng to
;                 fitgains or fitgains,/fit_reflectivity and with
;                 npar=2.
;
;    2013-12-22 : MGL. New keyword: download_all. Make downloading log
;                 files only the default.
;
;    2014-01-08 : MGL. Don't do downloading of log files
;                 directly, put the command to do it in the script
;                 file. Add isodate to the config file.
;
;    2014-01-09 : MGL. Bugfix isodate in config file.
;
;    2014-01-10 : MGL. The download command is now a method, write it
;                 in that form in the script.; 
;
;    2014-01-22 : MGL. Adapt to string functions moved to the str_
;                 namespace. 
;
;    2014-01-23 : MGL. No need to give /lapalma keyword, we'll
;                 know by examining the host name.
;
;    2014-09-08 : MGL. Changed the wording of comments on the fitgains
;                 method written to doit.pro.
;
;    2015-04-07 : MGL. Changed the default path for data in Stockholm
;                 to "/mnt/sand??/" (and not its subdirectory
;                 "Incoming/".
;
;    2015-08-12 : THI. Use demodulate rather than (recently renamed)
;                 demodulate2. 
;
;
;-
pro red_setupworkdir, root_dir = root_dir $
                      , out_dir = out_dir $
                      , cfgfile = cfgfile $
                      , scriptfile = scriptfile $
                      , download_all = download_all $
                      , sand = sand $
                      , date = date $
                      , stockholm = stockholm $
                      , lapalma = lapalma

  if n_elements(cfgfile) eq 0 then cfgfile = 'config.txt'
  if n_elements(scriptfile) eq 0 then scriptfile = 'doit.pro'

  if n_elements(out_dir) eq 0 then out_dir = getenv('PWD')  
  if ~strmatch(out_dir,'*/') then out_dir += '/'
  
  if n_elements(date) eq 0 then begin
     ;; Date not specified. Does root_dir include the date?
     if n_elements(root_dir) eq 0 then begin
        date_known = 0
     endif else begin
        pos = stregex(root_dir,'/[0-9][0-9][0-9][0-9][.-][0-9][0-9][.-][0-9][0-9]')
        if pos ne -1 then begin
           ;; Date present in root_dir
           date = strmid(root_dir, pos+1, 10)
           date_known = 1
        endif else begin
           date_known = 0
        endelse
     endelse
  endif

  if n_elements(date) eq 0 then begin
     ;; Get the date from out_dir?
     pos = stregex(out_dir,'/[0-9][0-9][0-9][0-9][.-][0-9][0-9][.-][0-9][0-9]')
     if pos eq -1 then begin
        print, 'sst_makeconfig : No date in either root_dir or out_dir.'
        retall
     endif
     date = strmid(out_dir, pos+1, 10)
     date = red_strreplace(date, '-', '.', n = 2)
;     root_dir = root_dir+date+'/'
  endif

  isodate = red_strreplace(date, '.', '-', n = 2)

  date_momfbd = isodate
  date = red_strreplace(isodate, '-', '.', n = 2)

  ;; Where to look for data?
  known_site = keyword_set(lapalma) $
               or keyword_set(stockholm)
  
  if ~known_site and n_elements(root_dir) eq 0 then begin
     hostname = getenv('HOSTNAME')
     if strpos(hostname,'royac.iac.es') ne -1 then begin
        lapalma = 1
     endif else begin
        stockholm = 1
     endelse
  endif
  
  if keyword_set(lapalma) then begin
     search_dir = "/data/disk?/*/"
     found_dir = file_search(search_dir+date, count = Nfound)
     if Nfound eq 0 then begin
        search_dir = "/data/camera?/*/"
        found_dir = file_search(search_dir+date, count = Nfound)
     endif
  endif else begin
     if keyword_set(stockholm) then begin
        search_dir = "/mnt/sand??/"
     endif else begin
         if ~strmatch(root_dir,'*/') then root_dir += '/'
        search_dir = root_dir
     endelse
        ;;; is search_dir already the one we look for?
     IF file_basename(search_dir) EQ date THEN BEGIN
         found_dir = search_dir
         Nfound = 1
     ENDIF ELSE $
       found_dir = file_search(search_dir+date, count = Nfound)
  endelse
  
  if Nfound eq 1 then begin
     root_dir = found_dir
     if ~strmatch(root_dir,'*/') then root_dir += '/'
  endif else if Nfound eq 0 then begin
     print, 'Cannot fine data from '+date+' in '+search_dir
     return
  endif else begin
     print, 'The root directory is not unique.'
     print, '"'+search_dir+'" matches: '
     print, found_dir
     stop
     ;; And here we need to figure out what to do...
     ;; Could happen in La Palma
  endelse

  Ncpu = !cpu.hw_ncpu
  If Ncpu le 2 then Nthreads = 2 else Nthreads = round(Ncpu*.75) <20

;  ;; Download position log
;  pfile = 'positionLog_'+date
;  tmp = file_search(pfile, count = Nlog)
;  if Nlog eq 0 then spawn, "scp obs@royac27.royac.iac.es:/usr/turret/logs/position/"+pfile+" ./"

  ;; Open two files for writing. Use logical unit Clun for a Config
  ;; file and Slun for a Script file.
  openw, Clun, cfgfile, /get_lun
  openw, Slun, scriptfile, /get_lun

  ;; Specify the date in the config file, ISO format.
  print, 'Date'
  printf, Clun, '#'
  printf, Clun, '# --- Date'
  printf, Clun, '#'
  printf, Clun,'isodate = '+isodate

  ;; printf, Slun, '.r crispred'
  printf, Slun, 'a = crispred("config.txt")' 
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
  printf, Clun, 'cam_t = Crisp-T'
  printf, Clun, 'cam_r = Crisp-R'
  printf, Clun, 'cam_wb = Crisp-W'
  printf, Clun, '#'
  printf, Clun, 'root_dir = ' + root_dir
  printf, Clun, '#'

  print, 'Output'
  printf, Clun, '#'
  printf, Clun, '# --- Output'
  printf, Clun, '#'
  printf, Clun, 'out_dir = ' + out_dir

  print, 'Darks'
  printf, Clun, '#'
  printf, Clun, '# --- Darks'
  printf, Clun, '#'
  darkdirs = file_search(root_dir+'/dark*/*', count = Ndirs, /fold)
  for i = 0, Ndirs-1 do begin
     darksubdirs = file_search(darkdirs[i]+'/crisp*', count = Nsubdirs, /fold)
     if Nsubdirs gt 0 then begin
        printf, Clun, 'dark_dir = '+red_strreplace(darkdirs[i], root_dir, '')
     endif
  endfor
  printf, Slun, 'a -> sumdark, /check' 
  

  print, 'Flats'
  printf, Clun, '#'
  printf, Clun, '# --- Flats'
  printf, Clun, '#'
  flatdirs = file_search(root_dir+'/flat*/*', count = Ndirs, /fold)
  dirarr = strarr(1)
  Nflat = 0
  prefilters = strarr(1)
  Nprefilters = 0
  for i = 0, Ndirs-1 do begin
     flatsubdirs = file_search(flatdirs[i]+'/crisp*', count = Nsubdirs, /fold)
     if Nsubdirs gt 0 then begin
        dirarr = [dirarr, red_strreplace(flatdirs[i], root_dir, '')]
        Nflat += 1
        printf, Clun, 'flat_dir = '+red_strreplace(flatdirs[i], root_dir, '')
        ;; Camera dirs and wavelengths to print to script file
        camdirs = strarr(Nsubdirs)
        wavelengths = strarr(Nsubdirs)
        for idir = 0, Nsubdirs-1 do begin
           camdirs[idir] = (strsplit(flatsubdirs[idir],  '/',/extract,count=nn))[nn-1]
           if camdirs[idir] eq 'Crisp-W' then rf = '3' else rf = '5'
           spawn, "ls "+flatsubdirs[idir]+"|grep cam|wc -l", Nfiles
           if long(Nfiles) gt 0 then begin
              spawn, "ls "+flatsubdirs[idir]+"|grep cam|rev|cut -d. -f"+rf+"|uniq|rev", wls
              wavelengths[idir] = strjoin(wls, ' ')
           endif
        endfor
        ;; Print to script file
        printf, Slun, 'a -> setflatdir, root_dir+"' + red_strreplace(flatdirs[i], root_dir, '') $
                + '"  ; ' + strjoin(camdirs+' ('+wavelengths+')', ' ')
        printf, Slun, 'a -> sumflat, /check'
     endif else begin
        flatsubdirs = file_search(flatdirs[i]+'/*', count = Nsubdirs)
        for j = 0, Nsubdirs-1 do begin
           flatsubsubdirs = file_search(flatsubdirs[j]+'/crisp*', count = Nsubsubdirs, /fold)
           if Nsubsubdirs gt 0 then begin
              dirarr = [dirarr, red_strreplace(flatsubdirs[j], root_dir, '')]
              Nflat += 1
              printf, Clun, 'flat_dir = '+red_strreplace(flatsubdirs[j], root_dir, '')
              ;; Camera dirs and wavelengths to print to script file
              camdirs = strarr(Nsubsubdirs)
              wavelengths = strarr(Nsubsubdirs)
              for idir = 0, Nsubsubdirs-1 do begin
                 camdirs[idir] = (strsplit(flatsubsubdirs[idir],  '/',/extract,count=nn))[nn-1]
                 if camdirs[idir] eq 'Crisp-W' then rf = '3' else rf = '5'
                 spawn, "ls "+flatsubsubdirs[idir]+"|grep cam|wc -l", Nfiles      
                 if long(Nfiles) gt 0 then begin
                    spawn, "ls "+flatsubsubdirs[idir]+"|grep cam|rev|cut -d. -f"+rf+"|uniq|rev", wls
                    wavelengths[idir] = strjoin(wls, ' ')
                 endif
              endfor
              ;; Print to script file
              printf, Slun, 'a -> setflatdir, root_dir+"' + red_strreplace(flatsubdirs[j], root_dir, '') $
                      + '" ; ' + strjoin(camdirs+' ('+wavelengths+')', ' ')
              printf, Slun, 'a -> sumflat, /check'
              ;; The following does not give the right result, 6300
              ;; instead of 6302! 
              prefilters = [prefilters, (strsplit(flatdirs[i],'/',/extract, count = nn))[nn-1]]
              Nprefilters += 1
           endif
        endfor
     endelse
  endfor
  if Nprefilters eq 0 then begin
     ;; This can happen if flats were already summed in La Palma. Look
     ;; for prefilters in the summed flats directory instead.
     spawn, 'ls flats/cam*.flat | cut -d. -f2|sort|uniq', prefilters
     Nprefilters = n_elements(prefilters)
  endif else begin
     prefilters = prefilters[1:*]
  endelse

  printf, Slun, 'a -> makegains' 


  print, 'Pinholes'
  printf, Clun, '#'
  printf, Clun, '# --- Pinholes'
  printf, Clun, '#'
  pinhdirs = file_search(root_dir+'/pinh*/*', count = Ndirs, /fold)
  for i = 0, Ndirs-1 do begin
     pinhsubdirs = file_search(pinhdirs[i]+'/crisp*', count = Nsubdirs, /fold)
     if Nsubdirs gt 0 then begin
        printf, Clun, 'pinh_dir = '+red_strreplace(pinhdirs[i], root_dir, '')
        printf, Slun, 'a -> setpinhdir, root_dir+"'+red_strreplace(pinhdirs[i], root_dir, '')+'"'
;        printf, Slun, 'a -> sumpinh_new'
        printf, Slun, 'a -> sumpinh,/pinhole_align'
     endif else begin
        pinhsubdirs = file_search(pinhdirs[i]+'/*', count = Nsubdirs)
        for j = 0, Nsubdirs-1 do begin
           pinhsubsubdirs = file_search(pinhsubdirs[j]+'/crisp*', count = Nsubsubdirs, /fold)
           if Nsubsubdirs gt 0 then begin
              printf, Clun, 'pinh_dir = '+red_strreplace(pinhsubdirs[j], root_dir, '')
              printf, Slun, 'a -> setpinhdir, root_dir+"'+red_strreplace(pinhsubdirs[j], root_dir, '')+'"'
;              printf, Slun, 'a -> sumpinh_new'
              printf, Slun, 'a -> sumpinh,/pinhole_align'
           endif
        endfor
     endelse
  endfor
  
  print, 'Polcal'
  printf, Clun, '#'
  printf, Clun, '# --- Polcal'
  printf, Clun, '#'
  Npol = 0
  polcaldirs = file_search(root_dir+'/polc*/*', count = Ndirs, /fold)
  for i = 0, Ndirs-1 do begin
     polcalsubdirs = file_search(polcaldirs[i]+'/crisp*', count = Nsubdirs, /fold)
     if Nsubdirs gt 0 then begin
        printf, Clun, 'polcal_dir = '+red_strreplace(polcaldirs[i], root_dir, '')
        Npol += 1
        printf, Slun, 'a -> setpolcaldir, root_dir+"' + red_strreplace(polcaldirs[i], root_dir, '')+'"'
        printf, Slun, 'a -> sumpolcal, /check'
     endif else begin
        polcalsubdirs = file_search(polcaldirs[i]+'/*', count = Nsubdirs)
        for j = 0, Nsubdirs-1 do begin
           polcalsubsubdirs = file_search(polcalsubdirs[j]+'/crisp*', count = Nsubsubdirs, /fold)
           if Nsubsubdirs gt 0 then begin
              printf, Clun, 'polcal_dir = '+red_strreplace(polcalsubdirs[j], root_dir, '')
              Npol += 1
              printf, Slun, 'a -> setpolcaldir, root_dir+"' + red_strreplace(polcalsubdirs[j], root_dir, '')+'"'
              printf, Slun, 'a -> sumpolcal, /check' 
           endif
        endfor
     endelse
  endfor
  ;; Find out the prefilters for which polcal needs to be run
  if Npol gt 0 then begin
     printf, Slun, 'a -> polcalcube' 
     printf, Slun, "spawn, 'ls polcal_cubes/* | grep 3d.f | cut -d. -f 2 | sort| uniq', polprefs"
     printf, Slun, "for i = 0, n_elements(polprefs)-1 do a -> polcal, pref=polprefs[i], nthreads="+strtrim(Nthreads, 2)
  endif


  print, 'PSF scan'
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
           pfscansubsubdirs = file_search(pfscansubdirs[j]+'/crisp*', count = Nsubsubdirs, /fold)
           if Nsubsubdirs gt 0 then begin
              printf, Clun, '# pfscan_dir = '+red_strreplace(pfscansubdirs[j], root_dir, '')
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
  nonsciencedirs = [darkdirs, flatdirs, pinhdirs, polcaldirs, pfscandirs]
  sciencedirs = file_search(root_dir+'/*/*', count = Ndirs)
  dirarr = strarr(1)
  Nsci = 0

  for i = 0, Ndirs-1 do begin

     if total(sciencedirs[i] eq nonsciencedirs) eq 0 then begin

        sciencesubdirs = file_search(sciencedirs[i]+'/crisp*', count = Nsubdirs, /fold)
        if Nsubdirs gt 0 then begin
           dirarr = [dirarr, red_strreplace(sciencedirs[i], root_dir, '')]
           Nsci += 1
        endif else begin
           sciencesubdirs = file_search(sciencedirs[i]+'/*', count = Nsubdirs)
           for j = 0, Nsubdirs-1 do begin
              sciencesubsubdirs = file_search(sciencesubdirs[j]+'/crisp*', count = Nsubsubdirs, /fold)
              if Nsubsubdirs gt 0 then begin
                 dirarr = [dirarr, red_strreplace(sciencesubdirs[j], root_dir, '')]
                 Nsci += 1
              endif
           endfor               ; j
        endelse 
     endif
  endfor
  if Nsci gt 0 then printf, Clun, "data_dir = ['"+strjoin(dirarr[1:*], "','")+"']"

  printf, Slun, 'a -> link_data' 

  if Npol gt 0 then begin
     printf, Slun, 'a -> prepflatcubes          ; For polarimetry data sets' 
  endif
  if Npol ne Nprefilters then begin
     printf, Slun, 'a -> prepflatcubes_lc4      ; For non-polarimetry data sets' 
  endif


  printf, Slun, ''
  printf, Slun, 'a -> getalignclips_new' 
  printf, Slun, 'a -> getoffsets' 
  
  printf, Slun, ''
  printf, Slun, 'a -> pinholecalib'
  
  printf, Slun, ''
  printf, Slun, ';; -----------------------------------------------------'
  printf, Slun, ';; This is how far we should be able to run unsupervised'
  printf, Slun, 'stop'          
  printf, Slun, ''

  print, 'Descatter (not implemented yet)'
  printf, Clun, '#'
  printf, Clun, '# --- 8542 descatter'
  printf, Clun, '#'
  printf, Clun, '#descatter_dir = '
  printf, Clun, '#'

  printf, Slun, 'a -> fitgains, npar = 2, res=res' 
  printf, Slun, 'If you need per-pixel reflectivities for your analysis'
  printf, Slun, '(e.g. for atmospheric inversions) you can set the /fit_reflectivity'
  printf, Slun, 'keyword:'
  printf, Slun, '; a -> fitgains, npar = 3, res=res, /fit_reflectivity  '
  printf, Slun, 'However, running without /fit_reflectivity is safer. In should not'
  printf, Slun, 'be used for chromospheric lines like 6563 and 8542.'

  for iline = 0, Nprefilters-1 do begin
     if long(prefilters[iline]) gt 7700 then maybedescatter = ', /descatter' else maybedescatter = ''
     printf, Slun, "a -> sum_data_intdif, pref = '" + prefilters[iline] $
             + "', cam = 'Crisp-T', /verbose, /show, /overwrite" + maybedescatter
     printf, Slun, "a -> sum_data_intdif, pref = '" + prefilters[iline] $
             + "', cam = 'Crisp-R', /verbose, /show, /overwrite" + maybedescatter
     printf, Slun, "a -> make_intdif_gains3, pref = '" + prefilters[iline] $
             + "', min=0.1, max=4.0, bad=1.0, smooth=3.0, timeaver=1L, /smallscale"
     if strmid(prefilters[iline], 0, 2) eq '63' then begin
        printf, Slun, "a -> fitprefilter, fixcav = 2.0d, pref = '"+prefilters[iline]+"', shift=-0.5"
     endif else begin
        printf, Slun, "a -> fitprefilter, fixcav = 2.0d, pref = '"+prefilters[iline]+"'"
     endelse
     printf, Slun, "a -> prepmomfbd, /newgains, /wb_states, date_obs = " + date_momfbd $
             + ", numpoints = '88', outformat = 'MOMFBD', pref = '"+prefilters[iline]+"'"
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
