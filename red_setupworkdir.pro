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
;    search_dir : in, string
; 
;      The top directory of your saved data, with or wthout a date. Or
;      a regular expression that matches the top directory. Or an
;      array of directories and regular expressions. If not given,
;      setupdir will look for the root_dir based on the site.
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
;    date : in, optional, type=string, default='From out_dir if possible' 
; 
;      The date (in iso format) the data was collected.
; 
;    out_dir : in, optional, type=string, default='Current directory'
; 
;      The output directory to be used by crispred.
; 
;    exclude_chromis : in, optional, type=boolean
;
;       Set this to not setup for CHROMIS data only.
; 
;    exclude_crisp : in, optional, type=boolean
;
;       Set this to not setup for CRISP data. 
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
;    2015-09-01 : MGL. Changed faulty text output when handling prefilter
;                 scan data. 
;
;    2016-02-15 : MGL. Remove (incomplete) definition of descatter_dir
;                 from config file.
;
;    2016-02-16 : MGL. Don't use /descatter keyword. Call sumpinh,
;                 polcalcube, and polcal separately for prefilters so
;                 user can easily add /no_descatter if 7772
;                 backscatter data is not available for the cameras in
;                 use.
;
;    2016-02-17 : MGL. Remove keywords newgain and outformat from
;                 prepmomfbd call. Add quotes around the date in the
;                 same call. Discussion about fitgais in doit.pro is
;                 now proper IDL comments.
;
;    2016-02-18 : MGL. Add a commented-out /all to sum_data_intdif and
;                 make_intdif_gains3 calls. Add a commented-out
;                 /no_descatter to some method calls involving the
;                 7772 prefilter. Remove duplicate wavelengths in the
;                 lists given as comments to the setflatdir calls.
;
;    2016-02-22 : MGL. Use red_extractstates for getting the
;                 prefilters from the flatdirs. Collect prefilters
;                 from flat directories also for the case when the
;                 camera directories are directly below the time-stamp
;                 directories without a prefilter level in between.
;
;    2016-02-24 : MGL. Added another root_dir to search when in
;                 Stockholm. 
;
;    2016-05-03 : THI. Get prefilter from filenames instead of dirname
;                 so we don't rely on a specific directory structure. 
;
;    2016-05-04 : MGL. Removed keywords stockholm and lapalma. Renamed
;                 root_dir to search_dir. Cleaned up the search for a
;                 root_dir, decide where to look based on the
;                 hostname. Also cleaned up the date handling a bit. 
;
;    2016-05-05 : MGL. Separate the CRISP and CHROMIS setups into two
;                 subdirectories under out_dir. Clean up in the
;                 searching for subdirectories for darks and flats.  
;
;    2016-05-17 : MGL. Changed Stockholm search dirs to accommodate
;                 "sand15n" type mounted disk names. Removed some
;                 polarization and descatter things from the CHROMIS
;                 part. New keywords exclude_chromis and
;                 exclude_crisp.
;
;
;-
pro red_setupworkdir, search_dir = search_dir $
                      , out_dir = out_dir $
                      , cfgfile = cfgfile $
                      , scriptfile = scriptfile $
                      , download_all = download_all $
                      , date = date $
                      , setup_chromis = setup_chromis $
                      , setup_crisp = setup_crisp

  if n_elements(out_dir) eq 0 then out_dir = getenv('PWD')  
  if ~strmatch(out_dir,'*/') then out_dir += '/'


  if ~keyword_set(exclude_crisp) then crisp_dir = out_dir + 'CRISP/'
  if ~keyword_set(exclude_chromis) then chromis_dir = out_dir + 'CHROMIS/'
  
  if n_elements(cfgfile) eq 0 then cfgfile = 'config.txt'
  if n_elements(scriptfile) eq 0 then scriptfile = 'doit.pro'

  if n_elements(date) eq 0 then begin
     ;; Date not specified. Does search_dir include the date?
     if n_elements(search_dir) gt 0 then begin
        for i = 0, n_elements(search_dir)-1 do begin
           pos = stregex(search_dir[i],'/[12][0-9][0-9][0-9][.-][01][0-9][.-][0-3][0-9]')
           if pos ne -1 then begin
              ;; Date present in search_dir[i]
              date = strmid(search_dir[i], pos+1, 10)
              break
           endif
        endfor                  ; i
     endif                      ; search_dir given?
  endif                         ; date given?

  if n_elements(date) eq 0 then begin
     ;; Get the date from out_dir?
     pos = stregex(out_dir,'/[12][0-9][0-9][0-9][.-][01][0-9][.-][0-3][0-9]')
     if pos eq -1 then begin
        print, 'sst_makeconfig : No date in either root_dir or out_dir.'
        retall
     endif
     date = strmid(out_dir, pos+1, 10)
     date = red_strreplace(date, '-', '.', n = 2)
  endif

  ;; Just in case date became a (one-element) array.
  date = date[0]

  ;; Date in ISO format and with dots.
  isodate = red_strreplace(date, '.', '-', n = 2)
  date = red_strreplace(isodate, '-', '.', n = 2)

  if n_elements(search_dir) eq 0 then begin
 
     ;; No search_dir specified. Try to find out where we are from the
     ;; hostname and search for a root_dir based on that.
     
     spawn, 'hostname -f', hostname

     case 1 of
        strpos(hostname,'royac.iac.es') ne -1 : begin
           ;; At the SST in La Palma
           search_dir = "/data/disk?/*/"
           ;; We could search the camera directories as well
           ;; but then we'd end up with multiple root_dirs, which I'd
           ;; now like to disallow. - Mats
           ;; search_dir = ["/data/disk?/*/", "/data/camera?/*/"]
        end
        strpos(hostname,'astro.su.se') ne -1 : begin
           ;; At the ISP in Stockholm
           search_dir = ["/mnt/sand*/", "/mnt/sand*/Incoming/", "/mnt/sand*/Incoming/Checked/"]
        end
     endcase

     ;; Make sure search_dir ends with a slash before we append the date
     for i = 0, n_elements(search_dir)-1 do begin
        if ~strmatch(search_dir[i],'*/') then search_dir[i] += '/'
     endfor

  endif                         ; No search_dir given
  
  ;; We now have a search_dir, but it could be a regular expression or
  ;; an array of directories and/or regular expressions. It could
  ;; include the date but then again it might not.

  for i = 0, n_elements(search_dir)-1 do begin
     if ~strmatch(search_dir[i],'*/') then search_dir[i] += '/'
     if file_basename(search_dir[i]) ne date then search_dir[i] += date
  endfor
  
     
  ;; Now search
  found_dir = file_search(search_dir, count = Nfound)

  ;; Maybe some searches returned the same result?
  if Nfound gt 1 then begin
     found_dir = found_dir[uniq(found_dir,sort(found_dir))]
     Nfound = n_elements(found_dir)
  endif

  print, 'Looked for data from '+date+' in:'
  for i = 0, n_elements(search_dir)-1 do print, search_dir[i]
  
  case Nfound of
     0: begin
        print, "Didn't find any data."
        return
     end
     1: begin
        print, 'Found one dir, which we will use:'
        print, found_dir
        root_dir = found_dir
     end
     else: begin
        print, 'Found several possible locations:'
        for i = 0, Nfound-1 do print, found_dir[i]
        print, 'Please call red_setupworkdir again with one of them specified as search_dir.'
        ;; We'd have to do something else here if we want to
        ;; allow the camera directories in La Palma.
        return
     end
  endcase
  
  ;; Make sure root_dir ends with a slash.
  if ~strmatch(root_dir,'*/') then root_dir += '/'


  Ncpu = !cpu.hw_ncpu
  If Ncpu le 2 then Nthreads = 2 else Nthreads = round(Ncpu*.75) <20

;  ;; Download position log
;  pfile = 'positionLog_'+date
;  tmp = file_search(pfile, count = Nlog)
;  if Nlog eq 0 then spawn, "scp obs@royac27.royac.iac.es:/usr/turret/logs/position/"+pfile+" ./"


  
  ;; CHROMIS ---------------------------------------------------------------------------------------
  
  if ~keyword_set(exclude_chromis) then begin
  
     ;; Open two files for writing. Use logical unit Clun for a Config
     ;; file and Slun for a Script file.
     file_mkdir, chromis_dir
     openw, Clun, chromis_dir + cfgfile, /get_lun
     openw, Slun, chromis_dir + scriptfile , /get_lun

     ;; Specify the date in the config file, ISO format.
     print, 'Date'
     printf, Clun, '#'
     printf, Clun, '# --- Date'
     printf, Clun, '#'
     printf, Clun,'isodate = '+isodate

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
     printf, Clun, 'cam_w = Chromis-W'
     printf, Clun, 'cam_p = Chromis-P'
     printf, Clun, 'cam_n = Chromis-N'
     printf, Clun, '#'
     printf, Clun, 'root_dir = ' + root_dir
     printf, Clun, '#'

     print, 'Output'
     printf, Clun, '#'
     printf, Clun, '# --- Output'
     printf, Clun, '#'
     printf, Clun, 'out_dir = ' + chromis_dir

     print, 'Darks'
     printf, Clun, '#'
     printf, Clun, '# --- Darks'
     printf, Clun, '#'
     
     darksubdirs = red_find_instrumentdirs(root_dir, 'chromis', 'dark', count = Nsubdirs)
     if Nsubdirs gt 0 then begin
        darkdirs = file_dirname(darksubdirs)
        darkdirs = darkdirs[uniq(darkdirs, sort(darkdirs))]
        for idir = 0, n_elements(darkdirs)-1 do begin
           printf, Clun, 'dark_dir = '+red_strreplace(darkdirs[idir], root_dir, '')
           printf, Slun, 'a -> setdarkdir, root_dir+"' + red_strreplace(darkdirs[idir], root_dir, '')
           printf, Slun, 'a -> sumdark, /check'
        endfor                  ; idir
     endif                      ; Nsubdirs
     
     print, 'Flats'
     printf, Clun, '#'
     printf, Clun, '# --- Flats'
     printf, Clun, '#'

     flatsubdirs = red_find_instrumentdirs(root_dir, 'chromis', 'flat', count = Nsubdirs)
     if Nsubdirs gt 0 then begin
        ;; There are CHROMIS flats!

        ;; Directories with camera dirs below:
        flatdirs = file_dirname(flatsubdirs)
        flatdirs = flatdirs[uniq(flatdirs, sort(flatdirs))]
        Nflatdirs = n_elements(flatdirs)

        ;; Loop over the flatdirs, write each to the config file and
        ;; to the script file
        for idir = 0, Nflatdirs-1 do begin

           ;; Config file

           printf, Clun, 'flat_dir = '+red_strreplace(flatdirs[idir], root_dir, '')

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
              print, 'a -> setflatdir, root_dir+"' + red_strreplace(flatdirs[idir], root_dir, '')$
                      + '"  ; ' + camdirs+' ('+wavelengths+')'
              print, 'a -> sumflat, /check'

              printf, Slun, 'a -> setflatdir, root_dir+"' + red_strreplace(flatdirs[idir], root_dir, '')$
                      + '"  ; ' + camdirs+' ('+wavelengths+')'
              printf, Slun, 'a -> sumflat, /check'

              red_append, prefilters, wls

           endif                ; Nfiles
        endfor                  ; idir
     endif                      ; Nsubdirs

     if n_elements(prefilters) gt 0 then begin
        prefilters = prefilters[uniq(prefilters, sort(prefilters))]
        Nprefilters = n_elements(prefilters)
     endif

     if Nprefilters eq 0 then begin
        ;; This can happen if flats were already summed in La Palma. Look
        ;; for prefilters in the summed flats directory instead.
        spawn, 'ls flats/cam*.flat | cut -d. -f2|sort|uniq', prefilters
        Nprefilters = n_elements(prefilters)
     endif
     
     for ipref = 0, Nprefilters-1 do begin
        printf, Slun, "a -> makegains, pref='" + prefilters[ipref] $
                + "' "
     endfor

     print, 'Pinholes'
     printf, Clun, '#'
     printf, Clun, '# --- Pinholes'
     printf, Clun, '#'
     pinhdirs = file_search(root_dir+'/pinh*/*', count = Ndirs, /fold)
     for i = 0, Ndirs-1 do begin
        pinhsubdirs = file_search(pinhdirs[i]+'/chromis*', count = Nsubdirs, /fold)
        if Nsubdirs gt 0 then begin
           printf, Clun, 'pinh_dir = '+red_strreplace(pinhdirs[i], root_dir, '')
           printf, Slun, 'a -> setpinhdir, root_dir+"'+red_strreplace(pinhdirs[i], root_dir, '')+'"'
;        printf, Slun, 'a -> sumpinh_new'
           for ipref = 0, Nprefilters-1 do begin
              printf, Slun, "a -> sumpinh, /pinhole_align, pref='"+prefilters[ipref]+"'"
           endfor
        endif else begin
           pinhsubdirs = file_search(pinhdirs[i]+'/*', count = Nsubdirs)
           for j = 0, Nsubdirs-1 do begin
              pinhsubsubdirs = file_search(pinhsubdirs[j]+'/chromis*', count = Nsubsubdirs, /fold)
              if Nsubsubdirs gt 0 then begin
                 printf, Clun, 'pinh_dir = '+red_strreplace(pinhsubdirs[j], root_dir, '')
                 printf, Slun, 'a -> setpinhdir, root_dir+"'+red_strreplace(pinhsubdirs[j], root_dir, '')+'"'
;              printf, Slun, 'a -> sumpinh_new'
                 for ipref = 0, Nprefilters-1 do begin
                    printf, Slun, "a -> sumpinh, /pinhole_align, pref='"+prefilters[ipref]+"'" 
                 endfor
              endif
           endfor
        endelse
     endfor

     
     print, 'Prefilter scan'
     printf, Clun, '#'
     printf, Clun, '# --- Prefilter scan'
     printf, Clun, '#'
     Npfs = 0
     pfscandirs = file_search(root_dir+'/pfscan*/*', count = Ndirs, /fold)
     for i = 0, Ndirs-1 do begin
        pfscansubdirs = file_search(pfscandirs[i]+'/chromis*', count = Nsubdirs, /fold)
        if Nsubdirs gt 0 then begin
           printf, Clun, '# pfscan_dir = '+red_strreplace(pfscandirs[i], root_dir, '')
           Npfs += 1
        endif else begin
           pfscansubdirs = file_search(pfscandirs[i]+'/*', count = Nsubdirs)
           for j = 0, Nsubdirs-1 do begin
              pfscansubsubdirs = file_search(pfscansubdirs[j]+'/chromis*', count = Nsubsubdirs, /fold)
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
     nonsciencedirs = [darkdirs, flatdirs, pinhdirs, pfscandirs]
     sciencedirs = file_search(root_dir+'/*/*', count = Ndirs)

     for i = 0, Ndirs-1 do begin

        if total(sciencedirs[i] eq nonsciencedirs) eq 0 then begin
           sciencesubdirs = file_search(sciencedirs[i]+'/chromis*', count = Nsubdirs, /fold)
           if Nsubdirs gt 0 then begin
              red_append, dirarr, red_strreplace(sciencedirs[i], root_dir, '')
           endif else begin
              sciencesubdirs = file_search(sciencedirs[i]+'/*', count = Nsubdirs)
              for j = 0, Nsubdirs-1 do begin
                 sciencesubsubdirs = file_search(sciencesubdirs[j]+'/chromis*', count = Nsubsubdirs, /fold)
                 if Nsubsubdirs gt 0 then begin
                    red_append, dirarr, red_strreplace(sciencesubdirs[j], root_dir, '')
                 endif
              endfor            ; j
           endelse 
        endif
     endfor
     if n_elements(dirarr) gt 0 then printf, Clun, "data_dir = ['"+strjoin(dirarr, "','")+"']"

     printf, Slun, 'a -> link_data' 
     
     for ipref = 0, Nprefilters-1 do begin
        printf, Slun, "a -> prepflatcubes_lc4, pref='"+prefilters[ipref]+"'"
     endfor                     ; ipref
     
     

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

     
     free_lun, Clun
     free_lun, Slun

  endif                         ; setup_chromis

  
  ;; CRISP -----------------------------------------------------------------------------------------

  if ~keyword_set(exclude_crisp) then begin
     
     ;; Open two files for writing. Use logical unit Clun for a Config
     ;; file and Slun for a Script file.
     file_mkdir, crisp_dir
     openw, Clun, crisp_dir + cfgfile, /get_lun
     openw, Slun, crisp_dir + scriptfile , /get_lun

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
     printf, Clun, 'out_dir = ' + crisp_dir

     print, 'Darks'
     printf, Clun, '#'
     printf, Clun, '# --- Darks'
     printf, Clun, '#'
     
     darksubdirs = red_find_instrumentdirs(root_dir, 'crisp', 'dark', count = Nsubdirs)
     if Nsubdirs gt 0 then begin
        darkdirs = file_dirname(darksubdirs)
        darkdirs = darkdirs[uniq(darkdirs, sort(darkdirs))]
        for idir = 0, n_elements(darkdirs)-1 do begin
           printf, Clun, 'dark_dir = '+red_strreplace(darkdirs[idir], root_dir, '')
           printf, Slun, 'a -> setdarkdir, root_dir+"' + red_strreplace(darkdirs[idir], root_dir, '')
           printf, Slun, 'a -> sumdark, /check'
        endfor                  ; idir
     endif                      ; Nsubdirs
     
     print, 'Flats'
     printf, Clun, '#'
     printf, Clun, '# --- Flats'
     printf, Clun, '#'

     flatsubdirs = red_find_instrumentdirs(root_dir, 'crisp', 'flat', count = Nsubdirs)
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

           printf, Clun, 'flat_dir = '+red_strreplace(flatdirs[idir], root_dir, '')

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
              print, 'a -> setflatdir, root_dir+"' + red_strreplace(flatdirs[idir], root_dir, '')$
                      + '"  ; ' + camdirs+' ('+wavelengths+')'
              print, 'a -> sumflat, /check'

              printf, Slun, 'a -> setflatdir, root_dir+"' + red_strreplace(flatdirs[idir], root_dir, '')$
                      + '"  ; ' + camdirs+' ('+wavelengths+')'
              printf, Slun, 'a -> sumflat, /check'

              red_append, prefilters, wls

           endif                ; Nfiles
        endfor                  ; idir
     endif                      ; Nsubdirs

     
     if n_elements(prefilters) gt 0 then begin
        prefilters = prefilters[uniq(prefilters, sort(prefilters))]
        Nprefilters = n_elements(prefilters)
     endif

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
           printf, Slun, 'a -> setpinhdir, root_dir+"'+red_strreplace(pinhdirs[i], root_dir, '')+'"'
;        printf, Slun, 'a -> sumpinh_new'
           for ipref = 0, Nprefilters-1 do begin
              printf, Slun, "a -> sumpinh, /pinhole_align, pref='"+prefilters[ipref]+"'" $
                      + maybe_nodescatter[ipref]
           endfor
        endif else begin
           pinhsubdirs = file_search(pinhdirs[i]+'/*', count = Nsubdirs)
           for j = 0, Nsubdirs-1 do begin
              pinhsubsubdirs = file_search(pinhsubdirs[j]+'/crisp*', count = Nsubsubdirs, /fold)
              if Nsubsubdirs gt 0 then begin
                 printf, Clun, 'pinh_dir = '+red_strreplace(pinhsubdirs[j], root_dir, '')
                 printf, Slun, 'a -> setpinhdir, root_dir+"'+red_strreplace(pinhsubdirs[j], root_dir, '')+'"'
;              printf, Slun, 'a -> sumpinh_new'
                 for ipref = 0, Nprefilters-1 do begin
                    printf, Slun, "a -> sumpinh, /pinhole_align, pref='"+prefilters[ipref]+"'" $
                            + maybe_nodescatter[ipref]
                 endfor
              endif
           endfor
        endelse
     endfor
     
     print, 'Polcal'
     printf, Clun, '#'
     printf, Clun, '# --- Polcal'
     printf, Clun, '#'
;  Npol = 0
     polcaldirs = file_search(root_dir+'/polc*/*', count = Npol, /fold)
     if Npol gt 0 then begin
        polprefs = file_basename(polcaldirs)
        for i = 0, Npol-1 do begin
           polcalsubdirs = file_search(polcaldirs[i]+'/crisp*', count = Nsubdirs, /fold)
           if Nsubdirs gt 0 then begin
              printf, Clun, 'polcal_dir = '+red_strreplace(polcaldirs[i], root_dir, '')
;        Npol += 1
              printf, Slun, 'a -> setpolcaldir, root_dir+"' + red_strreplace(polcaldirs[i], root_dir, '')+'"'
              printf, Slun, 'a -> sumpolcal, /check'
           endif else begin
              polcalsubdirs = file_search(polcaldirs[i]+'/*', count = Nsubdirs)
              for j = 0, Nsubdirs-1 do begin
                 polcalsubsubdirs = file_search(polcalsubdirs[j]+'/crisp*', count = Nsubsubdirs, /fold)
                 if Nsubsubdirs gt 0 then begin
                    printf, Clun, 'polcal_dir = '+red_strreplace(polcalsubdirs[j], root_dir, '')
;              Npol += 1
                    printf, Slun, 'a -> setpolcaldir, root_dir+"' + red_strreplace(polcalsubdirs[j], root_dir, '')+'"'
                    printf, Slun, 'a -> sumpolcal, /check' 
                 endif
              endfor            ; j
           endelse
        endfor                  ; i

        for ipref = 0, Npol-1 do begin
           printf, Slun, "a -> polcalcube, pref='"+polprefs[ipref]+"' " $
                   + maybe_nodescatter[ipref]
           printf, Slun, "a -> polcal, pref='"+polprefs[ipref]+"', nthreads=" $
                   + strtrim(Nthreads, 2)
        endfor                  ; ipref
        
     endif else begin
        polprefs = ''
     endelse                    ; Npol

     
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

     for i = 0, Ndirs-1 do begin

        if total(sciencedirs[i] eq nonsciencedirs) eq 0 then begin
           sciencesubdirs = file_search(sciencedirs[i]+'/crisp*', count = Nsubdirs, /fold)
           if Nsubdirs gt 0 then begin
              red_append, dirarr, red_strreplace(sciencedirs[i], root_dir, '')
           endif else begin
              sciencesubdirs = file_search(sciencedirs[i]+'/*', count = Nsubdirs)
              for j = 0, Nsubdirs-1 do begin
                 sciencesubsubdirs = file_search(sciencesubdirs[j]+'/crisp*', count = Nsubsubdirs, /fold)
                 if Nsubsubdirs gt 0 then begin
                    red_append, dirarr, red_strreplace(sciencesubdirs[j], root_dir, '')
                 endif
              endfor            ; j
           endelse 
        endif
     endfor
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
     endfor                     ; ipref
     
     

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
  endif                         ; setup_crisp

  

  ;; Write message and then we are done.
  
  print
  
  if ~keyword_set(exclude_crisp) then begin
     print, 'CRISP setup in ' + crisp_dir
  endif else begin
     print, 'No CRISP data'
  endelse
  
  if ~keyword_set(exclude_chromis) then begin
     print, 'CHROMIS setup in ' + chromis_dir
  endif else begin
     print, 'No CHROMIS data'
  endelse
  
end
