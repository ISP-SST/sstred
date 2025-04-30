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
;    ampm_cutoff : in, optional, type=string, default='13:00:00'
;
;      Data collected before this time belongs to AM observations,
;      everything else to PM.
;
;    calibrations_only : in, optional, type=boolean
;
;      Set up to process calibration data only.
;
;    cfgfile : in, optional, type=string, default='config.txt'
;
;      The name of the generated config file.
;
;    date : in, optional, type=string, default='From out_dir if possible'
;
;      The date (in iso format) the data was collected.
;
;    no_observer_metadata : in, optional, type=boolean
;
;      Don't do anything to get OBSERVER metadata. 
;
;    no_lapalma : in, optional, type=boolean
;
;      By default, if you are in La Palma, the workdirs will be set up
;      to sum all calibration data into timestamped directories, like
;      flats/hh:mm:ss plus make quicklooks. Set this keyword to
;      generate regular workdirs.
;
;    old_dir : in, optional, type = string
;
;      Copy files from this directory, in particular summed
;      calibration data. Useful if you summed the calibration data in
;      La Palma.
;
;    out_dir : in, optional, type=string, default='Current directory'
;
;      The output directory, under which instrument-specific
;      directories are created.
;
;    scriptfile : in, optional, type=string, default='doit.pro'
;
;      The name of the generated script file. The script file can be
;      run in an idl session with "@doit.pro" (assuming the default
;      name). It will perform the basic things, like co-adding of
;      darks, flats, etc. Later commands, that involve human
;      interaction are present in the file but commented out.
;
;    search_dirs : in, optional, type = array of strings
;
;      The top directory of your saved data, with or wthout a date. Or
;      a regular expression that matches the top directory. Or an
;      array of directories and regular expressions. If not given,
;      setupdir will look for the root_dir based on the site.
;
;
; :History:
;
;   2013-07-10 : MGL. Will now get a date from out_dir if none is
;                present in root_dir.
;
;   2013-08-29 : MGL. Any subdirectory that is not a known calibration
;                data directory (darks, flats, etc.) is a potential
;                science data directory.
;
;   2013-08-30 : MGL. Take care of prefilter scans.
;
;   2013-11-25 : MGL. Renamed from original name sst_makeconfig and
;                moved to Crispred repository. You can now use regular
;                expressions for root_dir. New keywords "date",
;                "lapalma" and "stockholm". Made root_dir a keyword.
;
;   2013-11-26 : MGL. Various improvements. The "new" flat field
;                procedure. Find out which cameras and wavelengths are
;                present in the raw flats directories.
;
;   2013-11-29 : MGL. Changed the order of some commands in the
;                doit.pro file. Add "/descatter" keyword to
;                sum_data_intdif call only for wavelengths > 7700.
;
;   2013-12-02 : MGL. Bugfixes. Add slash at end of root_dir. Can now
;                deal with empty raw flats directories. Gets the
;                prefilters from the summed flats instead.
;
;   2013-12-02 : MGL. Move prepflatcubes[_lc4] to the unsupervised
;                part. Default parameters for fitgains[_ng]. Deal with
;                data sets with more than a single polcal directory
;                and prefilter. Set number of threads to use based on
;                the number of CPUs.
;
;   2013-12-09 : Pit. Allow root_dir to be the actual date directory.
;
;   2013-12-19 : MGL. Download SST logfiles and some other data from
;                the web.
;
;   2013-12-20 : MGL. Change calls from fitgains or fitgains_ng to
;                fitgains or fitgains,/fit_reflectivity and with
;                npar=2.
;
;   2013-12-22 : MGL. New keyword: download_all. Make downloading log
;                files only the default.
;
;   2014-01-08 : MGL. Don't do downloading of log files directly, put
;                the command to do it in the script file. Add isodate
;                to the config file.
;
;   2014-01-09 : MGL. Bugfix isodate in config file.
;
;   2014-01-10 : MGL. The download command is now a method, write it
;                in that form in the script.;
;
;   2014-01-22 : MGL. Adapt to string functions moved to the str_
;                namespace.
;
;   2014-01-23 : MGL. No need to give /lapalma keyword, we'll know by
;                examining the host name.
;
;   2014-09-08 : MGL. Changed the wording of comments on the fitgains
;                method written to doit.pro.
;
;   2015-04-07 : MGL. Changed the default path for data in Stockholm
;                to "/mnt/sand??/" (and not its subdirectory
;                "Incoming/".
;
;   2015-08-12 : THI. Use demodulate rather than (recently renamed)
;                demodulate2.
;
;   2015-09-01 : MGL. Changed faulty text output when handling
;                prefilter scan data.
;
;   2016-02-15 : MGL. Remove (incomplete) definition of descatter_dir
;                from config file.
;
;   2016-02-16 : MGL. Don't use /descatter keyword. Call sumpinh,
;                polcalcube, and polcal separately for prefilters so
;                user can easily add /no_descatter if 7772 backscatter
;                data is not available for the cameras in use.
;
;   2016-02-17 : MGL. Remove keywords newgain and outformat from
;                prepmomfbd call. Add quotes around the date in the
;                same call. Discussion about fitgais in doit.pro is
;                now proper IDL comments.
;
;   2016-02-18 : MGL. Add a commented-out /all to sum_data_intdif and
;                make_intdif_gains3 calls. Add a commented-out
;                /no_descatter to some method calls involving the 7772
;                prefilter. Remove duplicate wavelengths in the lists
;                given as comments to the setflatdir calls.
;
;   2016-02-22 : MGL. Use red_extractstates for getting the prefilters
;                from the flatdirs. Collect prefilters from flat
;                directories also for the case when the camera
;                directories are directly below the time-stamp
;                directories without a prefilter level in between.
;
;   2016-02-24 : MGL. Added another root_dir to search when in
;                Stockholm.
;
;   2016-05-03 : THI. Get prefilter from filenames instead of dirname
;                so we don't rely on a specific directory structure.
;
;   2016-05-04 : MGL. Removed keywords stockholm and lapalma. Renamed
;                root_dir to search_dir. Cleaned up the search for a
;                root_dir, decide where to look based on the hostname.
;                Also cleaned up the date handling a bit.
;
;   2016-05-05 : MGL. Separate the CRISP and CHROMIS setups into two
;                subdirectories under out_dir. Clean up in the
;                searching for subdirectories for darks and flats.
;
;   2016-05-17 : MGL. Changed Stockholm search dirs to accommodate
;                "sand15n" type mounted disk names. Removed some
;                polarization and descatter things from the CHROMIS
;                part. New keywords exclude_chromis and exclude_crisp.
;                Use the new camera in the config file. Make an
;                instance of chromisred rather than crispred in the
;                CHROMIS script.
;
;   2016-05-18 : MGL. Use sumdark's dirs keyword instead of
;                setdarkdir.
;
;   2016-05-30 : MGL. Don't zero Sarnoff tap borders for CHROMIS
;                cameras.
;
;   2016-06-04 : MGL. CHROMIS image scale.
;
;   2016-08-12 : MGL. Detect the presence of CRISP and CHROMIS data.
;                Remove keywords exclude_crisp and exclude_chromis.
;
;   2016-08-23 : THI. Rename camtag to detector and channel to camera,
;                so the names match those of the corresponding
;                SolarNet keywords.
;
;   2016-09-08 : MGL. Add analyze_directories and red_plot_r0 to the
;                script file, but commented out.
;
;   2016-09-19 : MGL. Add hrz_zeropoint.
;
;   2017-01-25 : MGL. Added (nominal) diversity.
;
;   2017-03-03 : MGL. Search directories based on ip address rather
;                than hostname. Check for existing workdirs.
;
;   2017-03-06 : MGL. Add storing of metadata in the work directories.
;
;   2017-03-07 : MGL. Remove calls to getalignclips and getoffsets
;                methods, not needed with Tomas' getalignclips.
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
;   2020-09-10 : MGL. New keyword old_dir.
;
;   2020-10-25 : JdlCR. old_dir was not defined as a keyword.
; 
;   2021-01-25 : MGL. New keyword no_observer_metadata.
; 
;   2022-08-02 : MGL. Support CRISP's new cameras as well as future
;                CRISP2 instrument.
; 
;   2024-07-04 : MGL. Automatically copy summed calibration data from
;                a reduc/ subdirectory in the raw data root directory.
;
;   2025-03-27 : MGL. New keyword no_lapalma. 
;
;   2025-04-03 : MGL. New keyword ampm_cutoff 
;   
;-
pro red_setupworkdir, ampm_cutoff = ampm_cutoff $
                      , cfgfile = cfgfile $
                      , calibrations_only = calibrations_only $
                      , date = date $
                      , instruments = instruments $
                      , no_observer_metadata = no_observer_metadata $
                      , no_lapalma = no_lapalma $
                      , lapalma_setup = lapalma_setup $
                      , old_dir = old_dir $
                      , out_dir = out_dir $
                      , scriptfile = scriptfile $
                      , search_dirs = search_dirs

  ;; Name of this program
  inam = red_subprogram(/low, calling = inam1)

  if n_elements(ampm_cutoff) eq 0 then ampm_cutoff = '13:00:00'

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

  ;; Make sure CHROMIS is done before CRISP, this is so the CRISP
  ;; setup can peek into the CHROMIS directory and perhaps find
  ;; metadata. 
  instruments = instruments[sort(instruments)] 
  
  ;; If search_dirs are not given, pick them depending on the current location.
  if ( n_elements( search_dirs ) eq 0 ) then begin

    message, 'search_dirs are not given and are set depending on the ' + $
             'current computer location.', /informational

    red_currentsite, site = site, search_dirs = search_dirs, date = date

  endif else red_currentsite, site = site

  ;; Should we generate workdirs for standard La Palma use, i.e., with
  ;; calibrations data summing and quicklook only?
  if ~keyword_set(lapalma_setup) then lapalma_setup = ~keyword_set(no_lapalma) && n_elements(site) gt 0 && site eq 'La Palma'
  
  ;; search_dirs might be a single path or an array of paths that, in
  ;; turn, could be either a regular directory name or a regular
  ;; expression.
  for i = 0, n_elements( search_dirs ) - 1 do begin

    ;; Each search directory must end with a slash.
    if ~strmatch( search_dirs[ i ], '*/' ) then search_dirs[ i ] += '/'

    ;; Add the date directory at the end if it is not added yet.
    if ( file_basename( search_dirs[ i ] ) ne date and file_basename( search_dirs[ i ] ) ne isodate ) then begin
      search_dirs[ i ] += date
    endif

  endfor                        ; all search_dirs

  ;; Search all search_dirs for the specified date.  This might be slow on La
  ;; Palma for non-mounted /data/camera? or /data/disk? directories.
  found_dir = file_search( search_dirs, count = nfound_dirs )

  ;; Maybe some searches returned the same result?
  if ( nfound_dirs gt 1 ) then begin
    found_dir   = found_dir[ uniq( found_dir, sort( found_dir ) ) ]
    nfound_dirs = n_elements( found_dir )
  endif

  print, 'Looked for data from ' + date + ' in:'
  for i = 0, n_elements( search_dirs ) - 1 do print, search_dirs[ i ]

  case nfound_dirs of
    0: begin
      print, "Didn't find any data from " + date + '.'
      return
    end
    1: begin
      print, 'Found one dir, which we will use:'
      print, found_dir
    end
    else: begin
      print, 'Found several possible locations:'
      for i = 0, nfound_dirs - 1 do print, found_dir[ i ]
      print, 'Several doit.pro/config.txt will be generated. ' + $
             'Please run them separately'
    end
  endcase

  ;; Telescope location:
  ;; wikipedia, geo:28.759733,-17.880736
  ;; Mats C:        28.759693, -17.880757
  ;; wikipedia says altitude is 2360 m. Should add 20 m for height
  ;; of tower?
  ;; Altitude 2360 plus a few meters according to Google's
  ;; elevation service. Allowing for the tower we get 2380 m.
  obsgeo_xyz = round( red_obsgeo( 28.759733d, -17.880736d, 2380d ) )

  all_instruments =       [ 'CHROMIS', 'CRISP2', 'CRISP',  'TRIPPEL', 'SLITJAW', 'HeSP']
  all_regexps     = '*' + [ 'chromis', 'crisp2', 'crisp',  'spec',    'slit',    'HeSP'] + '*'
  all_classes     =       [ 'CHROMIS', 'CRISP2', 'CRISP2', 'TRIPPEL', 'SLITJAW', 'HESP']
  Ninstruments    = n_elements( all_instruments )

  dirs = file_search(found_dir+'/*', count = Ndirs, /test_directory)

  ;; Initialize list of existing directories for instruments
  existing_instrument_dirs = list()
  
  ;; Identify all instruments with existing data except CRISP with old cameras
  for iinstrument = 0, Ninstruments-1 do begin

    if Ndirs eq 0 then break
    
;    indx = where(strmatch(file_basename(dirs), all_regexps[iinstrument], /fold), Nmatch $
;                 , COMPLEMENT=indx_complement, NCOMPLEMENT=Nnomatch)
    indx = where(strmatch(dirs, all_regexps[iinstrument], /fold), Nmatch $
                 , COMPLEMENT=indx_complement, NCOMPLEMENT=Nnomatch)

    if Nmatch gt 0 then begin
      existing_instrument_dirs.add, dirs[indx]
      red_append, existing_instruments, all_instruments[iinstrument]
      red_append, existing_classes, all_classes[iinstrument]

      if Nnomatch eq 0 then begin
        Ndirs = 0
      endif else begin
        Ndirs = Nnomatch
        dirs = dirs[indx_complement]
      endelse
      
    endif
    
  endfor                        ; iinstrument
  
  ;; Identify CRISP with old cameras as well as non-instrument
  ;; directories. 
  if Ndirs gt 0 then begin

    ;; Clean from non-instrument directories
    non_instrument_dirs = ['reduc']
;    indx = where(strmatch(file_basename(dirs), non_instrument_dirs+'*', /fold), Nmatch $
;                 , COMPLEMENT=indx_complement, NCOMPLEMENT=Nnomatch)
    indx = where(strmatch(file_basename(dirs), non_instrument_dirs+'*', /fold), Nmatch $
                 , COMPLEMENT=indx_complement, NCOMPLEMENT=Nnomatch)

    print, 'Skipping non-instrument direcrories:'
    print, dirs[indx]
    
    ;; Anything not recognized by now should be CRISP with old cameras
    if Nnomatch gt 0 then begin
      existing_instrument_dirs.add, dirs[indx_complement]
      red_append, existing_instruments, 'CRISP'
      red_append, existing_classes, 'CRISP'
    endif
    
  endif
    
  if nfound_dirs gt 1 then stop ; Not supported yet - are multiple root directories ever needed?
  irootdir = 0
  root_dir = found_dir[ irootdir ]
  if strmid(root_dir,0,1,/rev) ne '/' then root_dir += '/'

  print, 'Will try to setup for the following instruments: ', instruments
  
  ;; Loop over wanted instruments
  for iinstrument = 0, n_elements(instruments)-1 do begin

    instrument = instruments[iinstrument]
    
    iexisting = (where(instrument eq existing_instruments, Nwhere))[0]

    if Nwhere eq 0 then begin
      print, 'No data for ', instrument
      continue
    endif

    class = existing_classes[iexisting]
    instrument_dirs = reform(existing_instrument_dirs[iexisting])

    case 1 of

      keyword_set(lapalma_setup) : begin
        ;; Setup for processing on La Palma
        print, 'La Palma setup up for ' + instrument + ' with class ' + class + '.'
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
      end
      
      keyword_set(calibrations_only) : begin
        ;; Probably not used anymore, keeping for backwards
        ;; compatibility.
        if ( instrument eq 'CRISP' || instrument eq 'CHROMIS' ) then begin
          print, 'Setting up for ' + instrument + $
                 ', calibration data processing only!'
          workdir = out_dir + instrument + '-calibrations/'
        endif else begin
          workdir = ''
        endelse
      end

      else : begin
        ;; Regular setup
        print, 'Setting up for ' + instrument + ' with class ' + class + '.'
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
      end
      
    endcase
    
    
;    if keyword_set( calibrations_only ) then begin
;      
;      if ( instrument eq 'CRISP' || instrument eq 'CHROMIS' ) then begin
;        
;        print, 'Setting up for ' + instrument + $
;               ', calibration data processing only!'
;        
;        workdir = out_dir + instrument + '-calibrations/'
;        
;      endif else begin
;        workdir = ''
;      endelse
;      
;    endif else begin ;; full processing mode (not calibration only).
;      
;      print, 'Setting up for ' + instrument + ' with class ' + class + '.'
;      
;      ;; We asked for this instrument and there seems to be data
;      odirs = file_search( out_dir + instrument + '*', count = Nodirs )
;      if Nodirs eq 0 then begin
;        workdir = out_dir + instrument + '/'
;      endif else begin
;        print
;        print, 'Existing ' + instrument + ' work dirs in ' + out_dir + ' :'
;        print, file_basename( odirs ), format = '(a0)'
;        workdir = ''
;        print, 'Use an existing directory or create a new one.'
;        read, 'Specify ' + instrument + ' workdir name: ', workdir
;        workdir = out_dir + workdir + '/'
;      endelse
;      
;    endelse
    
    file_mkdir, workdir

    if keyword_set(lapalma_setup) then begin
      ;; Write a warning message in 0README
      openw, rlun, /get_lun, workdir+'/0README'
      red_strflow, lun = rlun, width = 50 $
                   , ['This work directory was created for default La Palma processing,' $
                      , 'including summation of calibration data for transport to the' $
                      , 'home institute, and for making quicklooks. The summed calibrations' $
                      , 'data is stored in time-stamped subdirectories to the usual' $
                      , 'directories, i.e., flats/hh:mm:ss/ rather than just flats/.']

      red_strflow, lun = rlun, width = 50, indent = '   ' $
                   , ['If you are in La Palma and this is not what you wanted,' $
                      , 'please run red_setupworkdir with /no_lapalma.']

      red_strflow, lun = rlun, width = 50, indent = '   ' $
                   , ['The last summed darks will be linked to the darks/ directory' $ 
                      , 'and used when summing the other kinds of calibrations.' $
                      , 'Similarly, the last flats will be linked and used for polcal and pinholes.']

      red_strflow, lun = rlun, width = 50, indent = '   ' $
                   , ['If you have several sets of darks or several sets of flats of the same kind,' $
                      , 'you may be better off using regular workdirs and making sure' $
                      , 'to match the calibrations properly. Note that pinholes should be' $
                      , 'rather insensitive to what flats are used.']
      free_lun, rlun
    endif
    
    if keyword_set(calibrations_only) then begin
      ;; Write a warning message in 0README about the limitations of
      ;; these data.
      openw, rlun, /get_lun, workdir+'/0README'
      printf, rlun, 'This work directory was created to be used for the automatic collection and'
      printf, rlun, 'co-adding of calibration data, like darks, flats, etc., using the initial'
      printf, rlun, 'few commands in the data ordinary processing pipeline.'
      printf, rlun, ''
      printf, rlun, 'The summed data could be used, with certain limitations, when processing'
      printf, rlun, 'science data with the same pipeline. In fact, we hope that after some'
      printf, rlun, 'testing and evaluation, we can stop transferring and storing some of the'
      printf, rlun, 'raw calibration data, particularly the flats. But, because it is automatic,'
      printf, rlun, 'the summed data are not necessarily exactly the same as they would be after'
      printf, rlun, 'running the pipeline manually.'
      printf, rlun, ''
      printf, rlun, 'The main limitation comes from the fact that observers often collect'
      printf, rlun, 'several sets of the same kind of calibration data. Sometimes because'
      printf, rlun, 'science data are collected both in the morning and in the late afternoon,'
      printf, rlun, 'each set requiring their own calibrations. But also because of a failed'
      printf, rlun, 'attempt to collect calibration data, leaving a faulty and/or incomplete'
      printf, rlun, 'calibration data set on disk.'
      printf, rlun, ''
      printf, rlun, 'When running manually, operators can consult the observer logs and the'
      printf, rlun, 'correct calibration data can be selected for summing. In the automatic'
      printf, rlun, 'mode we instead sum data from all sets separately, storing the results'
      printf, rlun, 'in timestamp subdirectories below the ordinary directory for the'
      printf, rlun, 'particular kind of calibration data. This means a selection often has'
      printf, rlun, 'to be made, between different versions of the summed data, just like'
      printf, rlun, 'you would otherwise have to do for the raw calibration data before'
      printf, rlun, 'summing.'
      printf, rlun, ''
      printf, rlun, 'Because the summed and averaged flats are stored after dark subtraction,'
      printf, rlun, 'and there could be several versions of the darks, the summed and NOT'
      printf, rlun, 'averaged versions of the flats should be copied over to the ordinary work'
      printf, rlun, 'directory, and dark corrected with the selected dark version.'
      printf, rlun, ''
      printf, rlun, 'The summed pinholes are corrected for both dark and flat and no'
      printf, rlun, 'un-corrected version is stored here. This means you may end up with'
      printf, rlun, 'pinhole data that are dark and flat corrected with non-optimal or even'
      printf, rlun, 'faulty darks and/or flats. The raw pinhole data should therefore always be'
      printf, rlun, 'copied to the home institute, so the summing could be done again using'
      printf, rlun, 'the correct darks and flats. However, it is likely that merely using'
      printf, rlun, 'afternoon darks and flats for morning pinholes will work just fine.'
      printf, rlun, ''
      printf, rlun, 'The darks used here for the flats, as well as the darks and flats used'
      printf, rlun, 'for the pinholes, are always the version that is collected last. This way'
      printf, rlun, 'we at least avoid the common case with faulty data followed by correct data.'
      printf, rlun, ''
      printf, rlun, ''
      printf, rlun, ''
      printf, rlun, ''
      free_lun, rlun
    endif
    
    ;; Write string metadata.
    red_metadata_store, fname = workdir + '/info/metadata.fits',                        $
                        [ { keyword : 'OBSRVTRY',                                      $
                            value   : 'Observatorio del Roque de los Muchachos (ORM)',  $
                            comment : 'Name of observatory' },                          $
                          { keyword : 'TELESCOP',                                      $
                            value   : 'Swedish 1-meter Solar Telescope (SST)',          $
                            comment : 'Name of telescope' },                            $
                          { keyword : 'OBJECT',                                        $
                            value   : 'Sun',                                            $
                            comment : '' } ] ;
    
    year = long((strsplit(isodate, '-', /extract))[0])
    if isodate lt red_dates(tag = 'AO KL') then begin
      ;;  if year lt 2013 then begin
      ao_nmode = 37L
      modename = 'Karhunen-Loeve'
    endif else begin
      ao_nmode = 84L
      modename = 'Mirror'
    endelse
    
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

    if n_elements(old_dir) gt 0 then sum_dir = old_dir else begin
      ;; Look for already summed data in a reduc/ subdirectory.
      sum_dir = root_dir + 'reduc/' + instrument
      if ~file_test(sum_dir, /directory) then undefine, sum_dir
      case 1 of

        file_test(root_dir + 'reduc/AM/' + instrument, /directory) $
           && file_test(root_dir + 'reduc/PM' + instrument, /directory) : begin
          ;; Both reduc/AM/ and reduc/PM/ subdirectories. Ask user
          ;; which one they want to set up for.
          red_message, 'Found both reduc/AM/ and reduc/PM/ subdirectories with summed data.' $
                       + ' Which one do you want to set up for?'
          sel = red_select_subset(['AM', 'PM'], count = count, maxcount = 1, default = '')
          if count eq 0 then undefine, sum_dir else begin
            sum_dir = root_dir + 'reduc/'+sel+'/' + instrument
          endelse
        end

        file_test(root_dir + 'reduc/AM/' + instrument, /directory) : begin
          ;; Just the reduc/AM/ subdirectory. Use it!
          sum_dir = root_dir + 'reduc/AM/' + instrument
        end

        file_test(root_dir + 'reduc/PM/' + instrument, /directory) : begin
          ;; Just the reduc/PM/ subdirectory. Use it!
          sum_dir = root_dir + 'reduc/PM/' + instrument
        end

        file_test(root_dir + 'reduc/' + instrument, /directory) : begin
          ;; Just the reduc/ subdirectory. Use it!
          sum_dir = root_dir + 'reduc/' + instrument
        end

        else : undefine, sum_dir
        
      endcase
      
    endelse
    
    ;; Setup the different instruments.
    call_procedure, 'red_setupworkdir_' + class       $
                    , workdir, root_dir, config_file, script_file, isodate $                    
                    ;;, workdir, root_dir, config_file, script_file, isodate $     
                    , ampm_cutoff = ampm_cutoff $               
                    , calibrations_only = calibrations_only $
                    , lapalma_setup = lapalma_setup $
                    , no_observer_metadata = no_observer_metadata $
                    , old_dir = sum_dir
    
  endfor                        ; iinstrument    

end

red_setupworkdir, /no_obs, search_dirs = '/storage/Incoming/HeSP/2023.10.26/', instr = 'HeSP'



;red_setupworkdir, /no_obs, search_dirs = '/scratch/mats/crisp2-simulated-data/2016-09-19/2016-09-19/'
;red_setupworkdir, /no_obs, search_dirs = '/storage/Incoming/2022.08.24/', instr = 'CRISP'
;red_setupworkdir, /no_obs, search_dirs = '/storage/Incoming/2022.08.26/', instr = 'CRISP'

end
