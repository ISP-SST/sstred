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
;
; :Keywords:
;
;    old_dir : in, optional, type = string
;
;      Copy files from this directory, in particular summed
;      calibration data. Useful if you summed the calibration data in
;      La Palma.
;
;    search_dirs : in, optional, type = array of strings
;
;      The top directory of your saved data, with or wthout a date. Or
;      a regular expression that matches the top directory. Or an
;      array of directories and regular expressions. If not given,
;      setupdir will look for the root_dir based on the site.
;
;    calibrations_only : in, optional, type=boolean
;
;      Set up to process calibration data only.
;
;    cfgfile : in, optional, type=string, default='config.txt'
;
;      The name of the generated config file.
;
;    download_all : in, optional, type=boolean
;
;      Set this to download auxiliary data, like SDO/HMI images and AR
;      maps. Otherwise download SST log file only.
;      AVS: obsolete keyword, not in use anymore.
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
;   2017-07-05 : MGL. New keyword calibrations_only.
;
;   2017-08-11 : MGL. Write a 0README file in the created work
;                directories when run in /calibrations_only mode.
;
;   2017-08-14 : MGL. Use red_currentsite.
;
;   2017-08-14 : MGL. Use loops instead of duplicating code for
;                different instruments.
;
;   2017-08-21 : AVS. Call red_setupworkdir_* using call_procedure.
;
;   2017-08-23 : AVS. Loop over all root_dirs from found_dirs.
;
;   2017-08-23 : AVS. Function add_suffix is added.
;
;   2017-08-24 : MGL. Rename add_suffix to red_add_suffix and move to
;                a file of its own.
;
;   2017-08-24 : AVS. Make search_dir plural. If it is a single
;                string, IDL allows to use the array indexing such as
;                seach_dirs[ 0 ].
;
;   2017-11-27 : MGL. Set altitude to 2380 m.
;
;   2018-06-15 : MGL. Add fits keyword AO_NMODE.
;
;   2018-12-06 : MGL. Add CRISP to default instruments.
;
;   2019-07-01 : MGL. New keyword old_dir.
;
;-
pro red_setupworkdir_alt, calibrations_only = calibrations_only $
                          , cfgfile = cfgfile $
                          , old_dir = old_dir $
                          , link_old = link_old $
                          , date = date $
                          , download_all = download_all $
                          , instruments = instruments $
                          , out_dir = out_dir $
                          , scriptfile = scriptfile $
                          , search_dirs = search_dirs                    

  if n_elements(instruments) eq 0 then instruments = ['CHROMIS', 'CRISP']

  if n_elements(out_dir) eq 0 then out_dir = getenv('PWD')
  if ~strmatch(out_dir,'*/') then out_dir += '/'

  if n_elements(cfgfile) eq 0 then cfgfile = 'config.txt'
  if n_elements(scriptfile) eq 0 then scriptfile = 'doit.pro'

  if n_elements(date) eq 0 then begin
    ;; Date not specified.  Do search_dirs include the date?
    if n_elements( search_dirs ) gt 0 then begin
      for i = 0, n_elements( search_dirs ) - 1 do begin
        pos = stregex( search_dirs[ i ], '/[12][0-9][0-9][0-9][.-][01][0-9][.-][0-3][0-9]')
        if pos ne -1 then begin
          ;; Date present in search_dirs[ i ]
          date = strmid( search_dirs[ i ], pos + 1, 10 )
          break
        endif
      endfor                    ; i
    endif                       ; search_dir given?
  endif                         ; date given?

  if n_elements(date) eq 0 then begin
    ;; Get the date from out_dir?
    pos = stregex(out_dir,'/[12][0-9][0-9][0-9][.-][01][0-9][.-][0-3][0-9]')
    if pos eq -1 then begin
      print, 'red_setupworkdir : No date in either root_dir or out_dir.'
      retall
    endif
    date = strmid(out_dir, pos+1, 10)
    date = red_strreplace(date, '-', '.', n = 2)
  endif

  ;; Just in case date became a (one-element) array.
  date = date[0]

  ;; Date in ISO format and with dots.
  isodate = red_strreplace( date,    '.', '-', n = 2 )
  date    = red_strreplace( isodate, '-', '.', n = 2 )

  year = long((strsplit(isodate, '-', /extract))[0])
  if year lt 2013 then begin
    ao_nmode = 37L
    modename = 'Karhunen-Loeve'
  endif else begin
    ao_nmode = 84L
    modename = 'Mirror'
  endelse

  ;; Existing instruments and how to recognice them
  all_instruments =       [ 'CHROMIS', 'CRISP', 'TRIPPEL', 'SLITJAW' ]
  all_regexps     = '*' + [ 'chromis', 'crisp', 'spec',    'slit'    ] + '*'
  Ninstruments    = n_elements( all_instruments )

  ;; Telescope location:
  ;; wikipedia, geo:28.759733,-17.880736
  ;; Mats C:        28.759693, -17.880757
  ;; wikipedia says altitude is 2360 m. Should add 20 m for height
  ;; of tower?
  ;; Altitude 2360 plus a few meters according to Google's
  ;; elevation service. Allowing for the tower we get 2380 m.
  obsgeo_xyz = round( red_obsgeo( 28.759733d, -17.880736d, 2380d ) )

  ;; If search_dirs are not given, pick them depending on the current location.
  if ( n_elements( search_dirs ) eq 0 ) then begin

    message, 'search_dirs are not given and are set depending on the ' + $
      'current computer location.', /informational

    red_currentsite, site = site, search_dirs = search_dirs, date = date

  endif ; no search_dirs are given

  ;; search_dirs might be a single path or an array of paths that, in turn,
  ;; could be either a regular directory name or a regular expression.
  for i = 0, n_elements( search_dirs ) - 1 do begin

    ;; Each search directory must end with a slash.
    if ~strmatch( search_dirs[ i ], '*/' ) then search_dirs[ i ] += '/'

    ;; Add the date directory at the end if it is not added yet.
    if ( file_basename( search_dirs[ i ] ) ne date ) then begin
      search_dirs[ i ] += date
    endif

  endfor ; all search_dirs

  ;; Search all search_dirs for the specified date.  This might be slow on La
  ;; Palma for non-mounted /data/camera? or /data/disk? directories.
;  found_dir = file_search( search_dirs, count = nfound_dirs )

  for iinstr = 0, Ninstruments - 1 do begin

    ;; The current instrument and the corresponding regexp.
    instrument = all_instruments[ iinstr ]
    regexp     = all_regexps[     iinstr ]

    if max( strmatch( instruments, instrument ) ) eq 0 then continue
    
    ;; Camera directories are named Chromis-D, Chromis-N, or Chromis-W
    ;; for CHROMIS and Crisp-R, Crisp-T, or Crisp-W for CRISP and so
    ;; on. They are (mostly) nested at the third level from root_dir,
    ;; for example: 
    ;;   .../2017.05.02/CHROMIS-flats/8545/Chromis-W
    ;;   .../2017.05.02/Darks/12:09:18/Crisp-T
    ;; where .../2017.05.02 is the current root_dir.

    ;; Search 
    found_dir = file_search( search_dirs + '/*/??:??:??/'+regexp $
                             , /fold, count = nfound_dirs )
    
    for i = 0, n_elements(found_dir)-1 do $
       found_dir[i] = strjoin((strsplit(found_dir[i], '/', /extract))[0:-4], '/')
    found_dir = '/'+found_dir[uniq(found_dir, sort(found_dir))]

    nfound_dirs = n_elements(found_dir)
    if nfound_dirs eq 0 then begin
      print, "Didn't find any "+instrument+" data from " + date + '.'
      continue
    end

    ;; We asked for this instrument and there seems to be data
    odirs = file_search( out_dir + instrument + '*', count = Nodirs )
    if Nodirs eq 0 then begin
      workdir = out_dir + instrument + '/'
    endif else begin
      print
      print, 'Existing ' + instrument + ' work dirs in ' + out_dir + ' :'
      print, file_basename( odirs ), format = '(a0)'
      workdir = ''
      print, 'Use an existing directory or create a new one.'
      read, 'Specify ' + instrument + ' workdir name: ', workdir
      workdir = out_dir + workdir + '/'
    endelse
    file_mkdir, workdir

    ;; Write string metadata.
    red_metadata_store, fname = workdir + '/info/metadata.fits',                        $
                        [ { keyword : 'OBSRVTRY',                                       $
                            value   : 'Observatorio del Roque de los Muchachos (ORM)',  $
                            comment : 'Name of observatory' },                          $
                          { keyword : 'TELESCOP',                                       $
                            value   : 'Swedish 1-meter Solar Telescope (SST)',          $
                            comment : 'Name of telescope' },                            $
                          { keyword : 'OBJECT',                                         $
                            value   : 'Sun',                                            $
                            comment : '' } ] ;
    
    ;; Write numerical (LONG) metadata.
    red_metadata_store, fname = workdir + '/info/metadata.fits',                        $
                        [ { keyword : 'OBSGEO-Z',                                       $
                            value   : obsgeo_xyz[ 2 ],                                  $
                            comment : '[m] SST location' },                             $
                          { keyword : 'OBSGEO-Y',                                       $
                            value   : obsgeo_xyz[ 1 ],                                  $
                            comment : '[m] SST location' },                             $
                          { keyword : 'OBSGEO-X',                                       $
                            value   : obsgeo_xyz[ 0 ],                                  $
                            comment : '[m] SST location' },                             $
                          { keyword : 'AO_NMODE',                                       $
                            value   : ao_nmode,                                         $
                            comment : 'Number of AO corrected '+modename+' modes' } ] ;

    if nfound_dirs gt 1 then begin
      print, 'Found several possible locations:'
      for i = 0, nfound_dirs - 1 do print, found_dir[ i ]
      print, 'Several doit.pro/config.txt will be generated. ' + $
             'Please run them separately'
    endif

    ;; Loop over all found directories.
    for irootdir = 0, nfound_dirs - 1 do begin

      root_dir = found_dir[ irootdir ]
      ;; Make sure root_dir ends with a slash.
      if ~strmatch( root_dir, '*/' ) then root_dir += '/'

      print, 'Setting up for ' + instrument + ' data from '+root_dir+ ' in ' + workdir


      ;; If there are several root_dirs, a corresponding suffix must be added
      ;; to the doit.pro script and the config.txt file to separate them.
      ;;suffix = ( nfound_dirs gt 1 ) ? string( irootdir + 1, format = '( i01 )' ) : ''
      if ( nfound_dirs gt 1 ) then begin
        suffix = string( irootdir + 1, format = '( i01 )' )
      endif else begin
        suffix = ''
      endelse
      config_file = red_add_suffix( cfgfile,    suffix = suffix )
      script_file = red_add_suffix( scriptfile, suffix = suffix )

      
      ;; Setup the different instruments.
      call_procedure, 'red_setupworkdir_' + instrument $
                      , workdir, root_dir, config_file, script_file, isodate $
                      , calibrations_only = calibrations_only  $
                      , old_dir = old_dir

      
    endfor                      ; irootdir


  endfor                        ; iinstr
  

end
