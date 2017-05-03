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
;      The output directory, under which instrument-specific
;      directories are created.
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
;                 exclude_crisp. Use the new camera in the
;                 config file. Make an instance of chromisred rather
;                 than crispred in the CHROMIS script.  
;
;    2016-05-18 : MGL. Use sumdark's dirs keyword instead of
;                 setdarkdir.  
;
;    2016-05-30 : MGL. Don't zero Sarnoff tap borders for
;                 CHROMIS cameras.   
;
;    2016-06-04 : MGL. CHROMIS image scale.   
;
;    2016-08-12 : MGL. Detect the presence of CRISP and CHROMIS data.
;                 Remove keywords exclude_crisp and exclude_chromis.
;
;    2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                 so the names match those of the corresponding SolarNet
;                 keywords.
; 
;    2016-09-08 : MGL. Add analyze_directories and red_plot_r0 to the
;                 script file, but commented out.
; 
;    2016-09-19 : MGL. Add hrz_zeropoint.
;
;    2017-01-25 : MGL. Added (nominal) diversity.
;
;    2017-03-03 : MGL. Search directories based on ip address rather
;                 than hostname. Check for existing workdirs.
;
;    2017-03-06 : MGL. Add storing of metadata in the work
;                 directories.
;
;    2017-03-07 : MGL. Remove calls to getalignclips and getoffsets
;                 methods, not needed with Tomas' getalignclips.
;
;   2017-03-09 : MGL. Make Sun the OBJECT.
;
;   2017-03-13 : MGL. Moved CHROMIS and CRISP setups to separate
;                subprograms. Make ready for setting up also for
;                TRIPPEL and SLITJAW data.
;
;   2017-05-03 : MGL. Use hostname -I rather than hostname -i to
;                figure out where we are.
;
;-
pro red_setupworkdir, search_dir = search_dir $
                      , out_dir = out_dir $
                      , cfgfile = cfgfile $
                      , instruments = instruments $
                      , scriptfile = scriptfile $
                      , download_all = download_all $
                      , date = date

  if n_elements(instruments) eq 0 then instruments = ['CHROMIS']
  
  if n_elements(out_dir) eq 0 then out_dir = getenv('PWD')  
  if ~strmatch(out_dir,'*/') then out_dir += '/'

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
      endfor                    ; i
    endif                       ; search_dir given?
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
    ;; ip address and search for a root_dir based on that.
    
    spawn, 'hostname -I', ipaddress

    case 1 of
      strmatch(ipaddress,'*161.72.15.*') : begin
        ;; At the SST in La Palma
        search_dir = "/data/disk?/*/"
        ;; We could search the camera directories as well but then
        ;; we'd end up with multiple root_dirs, which I'd now like to
        ;; disallow. - Mats
      end
      strmatch(ipaddress,'*130.237.166.*') : begin
        ;; At the ISP in Stockholm
        search_dir = '/storage/sand*/' + ['', 'Incoming/', 'Incoming/Checked/']
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

  ;; Used for detecting instruments
  testdirs = file_search(root_dir+'/*/*/*', count = Ndirs)

  ;; Work directories

  if total(strmatch(testdirs,'*spec*',/fold)) and strmatch(instruments, 'TRIPPEL') then begin
    setup_trippel = 1
    odirs = file_search(out_dir + 'TRIPPEL*', count = Nodirs)
    if Nodirs eq 0 then begin
      trippel_dir = out_dir + 'TRIPPEL/'  
    endif else begin
      print
      print, 'Existing TRIPPEL work dirs in '+out_dir+' :'
      print, file_basename(odirs), format = '(a0)'
      trippel_dir = ''
      print, 'Use an existing directory or create a new one.'
      read, 'Specify TRIPPEL workdir name: ', trippel_dir
      trippel_dir = out_dir + trippel_dir + '/'
    endelse
    red_append, workdirs, trippel_dir
  endif else setup_trippel = 0

  if total(strmatch(testdirs,'*slit*',/fold)) and strmatch(instruments, 'TRIPPEL') then begin
    setup_slitjaw = 1
    odirs = file_search(out_dir + 'SLITJAW*', count = Nodirs)
    if Nodirs eq 0 then begin
      slitjaw_dir = out_dir + 'SLITJAW/'  
    endif else begin
      print
      print, 'Existing SLITJAW work dirs in '+out_dir+' :'
      print, file_basename(odirs), format = '(a0)'
      slitjaw_dir = ''
      print, 'Use an existing directory or create a new one.'
      read, 'Specify SLITJAW workdir name: ', slitjaw_dir
      slitjaw_dir = out_dir + slitjaw_dir + '/'
    endelse
    red_append, workdirs, slitjaw_dir
  endif else setup_slitjaw = 0
  
  if total(strmatch(testdirs,'*chromis*',/fold)) and strmatch(instruments, 'CHROMIS') then begin
    setup_chromis = 1
    odirs = file_search(out_dir + 'CHROMIS*', count = Nodirs)
    if Nodirs eq 0 then begin
      chromis_dir = out_dir + 'CHROMIS/'  
    endif else begin
      print
      print, 'Existing CHROMIS work dirs in '+out_dir+' :'
      print, file_basename(odirs), format = '(a0)'
      chromis_dir = ''
      print, 'Use an existing directory or create a new one.'
      read, 'Specify CHROMIS workdir name: ', chromis_dir
      chromis_dir = out_dir + chromis_dir + '/'
    endelse
    red_append, workdirs, chromis_dir
  endif else setup_chromis = 0

  if total(strmatch(testdirs,'*crisp*',/fold)) and strmatch(instruments, 'CRISP') then begin
    setup_crisp = 1
    odirs = file_search(out_dir + 'CRISP*', count = Nodirs)
    if Nodirs eq 0 then begin
      crisp_dir = out_dir + 'CRISP/'  
    endif else begin
      print, 'Existing CRISP work dirs in '+out_dir+' :'
      print, file_basename(odirs), format = '(a0)'
      crisp_dir = ''
      print, 'Use an existing directory or create a new one.'
      read, 'Specify CRISP workdir name: ', crisp_dir
      crisp_dir = out_dir + crisp_dir + '/'
    endelse
    red_append, workdirs, crisp_dir
  endif else setup_crisp = 0

  if n_elements(workdirs) eq 0 then return
  
  ;; Common setup tasks

  ;; Telescope location:
  ;; wikipedia, geo:28.759733,-17.880736
  ;; Mats C:        28.759693, -17.880757
  ;; wikipedia says altitude is 2360 m. Should add 20 m for height of tower?
  obsgeo_xyz = round(red_obsgeo(28.759733d,-17.880736d, 2360d))

  for idir = 0, n_elements(workdirs)-1 do begin 

    file_mkdir, workdirs[idir]

    ;; Write string metadata
    red_metadata_store, fname = workdirs[idir] + '/info/metadata.fits' $
                        , [{keyword:'OBSRVTRY', value:'Observatorio del Roque de los Muchachos (ORM)' $
                            , comment:'Name of observatory'} $
                           , {keyword:'TELESCOP', value:'Swedish 1-meter Solar Telescope (SST)' $
                              , comment:'Name of telescope'} $
                           , {keyword:'OBJECT', value:'Sun', comment:''} $
                          ]
    
    ;; Write numerical metadata
    red_metadata_store, fname = workdirs[idir] + '/info/metadata.fits' $
                        , [{keyword:'OBSGEO-Z', value:obsgeo_xyz[2] $
                            , comment:'[m] SST location'}, $
                           {keyword:'OBSGEO-Y', value:obsgeo_xyz[1] $
                            , comment:'[m] SST location'}, $
                           {keyword:'OBSGEO-X', value:obsgeo_xyz[0] $
                            , comment:'[m] SST location'}]
    
  endfor                        ; idir
  
  ;; Setup the different instruments
  ;; -----------------------------------------------------------------------------------------
  if setup_slitjaw then red_setupworkdir_slitjaw, slitjaw_dir, root_dir, cfgfile, scriptfile, isodate
  if setup_trippel then red_setupworkdir_trippel, trippel_dir, root_dir, cfgfile, scriptfile, isodate
  if setup_chromis then red_setupworkdir_chromis, chromis_dir, root_dir, cfgfile, scriptfile, isodate
  if setup_crisp then red_setupworkdir_crisp, crisp_dir, root_dir, cfgfile, scriptfile, isodate

  ;; Write message and then we are done.
  print
  if strmatch(instruments, 'TRIPPEL') then begin
    print, 'TRIPPEL setup in ' + trippel_dir
  endif else begin
    print, 'No TRIPPEL setup.'
  endelse
  if strmatch(instruments, 'SLITJAW') then begin
    print, 'SLITJAW setup in ' + slitjaw_dir
  endif else begin
    print, 'No slitjaw setup.'
  endelse
  if strmatch(instruments, 'CHROMIS') then begin
    print, 'CHROMIS setup in ' + chromis_dir
  endif else begin
    print, 'No CHROMIS setup.'
  endelse
  if strmatch(instruments, 'CRISP') then begin
    print, 'CRISP setup in ' + crisp_dir
  endif else begin
    print, 'No CRISP setup.'
  endelse
  
end
