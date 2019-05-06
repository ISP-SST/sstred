; docformat = 'rst'

;+
; Downloads logfiles and other useful auxiliary data when and if
; needed. Puts them in subdirectory downloads/.
; 
; :Categories:
;
;    SST observations
; 
; 
; :Author:
; 
;    Mats Löfdahl, 2013-12-16
; 
; 
; :Params:
; 
; :Keywords:
; 
;    date : in, out, optional, type=string
; 
;      A string representing the day for which the logged data is
;      wanted. If date is not set, then try to get the date from the
;      current directory name. The format has to be either yyyy.mm.dd
;      or yyyy-mm-dd.
;
;    all : in, optional, type=boolean
;
;      Try to download all types of data.
;
;    logs : in, optional, type=boolean
;
;      Try to download SST log files. This is the default operation in
;      case no particular target is specified and all is also not
;      specified. 
;
;    pig : in, optional, type=boolean
;
;      Try to download the SST/PIG log file. Returns with the
;      downloaded file name (or empty string in case download failed).
;
;    pathpig : out, optional, type=string
;
;      The path to where the pig log file was saved (or the empty
;      string).
;
;    r0 : in, optional, type=boolean
;
;      Try to download the SST/AO r0 log file. Returns with the
;      downloaded file name (or empty string in case download failed).
;
;    pathr0 : out, optional, type=string
;
;      The path to where the r0 log file was saved (or the empty
;      string).
;
;    turret : in, optional, type=boolean
;
;      Try to download the SST/Turret log file. Returns with the
;      downloaded file name (or empty string in case download failed).
;
;    pathturret : out, optional, type=string
;
;      The path to where the turret log file was saved (or the empty
;      string).
;
;    armap : in, optional, type=boolean
;
;      Try to download the Active Regions map. 
;
;    hmi : in, optional, type=boolean
;
;      Try to download HMI images and movies. 
;
;    overwrite : in, optional, type=boolean
;
;      Set this to download without checking if the file already exists
;
;    backscatter  : in, optional, type="integer array"
;
;      Set this to ["8542", "7772"] (or a subset thereof) to download
;      backscatter gain and psf for the corresponding prefilter(s).
;
; :History:
; 
;    2013-12-19 : MGL. Let red_geturl do more of the testing. And also
;                 make softlinks to the log files. Make the overwrite
;                 keyword work.
; 
;    2013-12-20 : MGL. Optionally return file names in r0file,
;                 turretfile, and pigfile keywords. Do not make soft
;                 links. 
;
;    2014-01-04 : MGL. For turret, now downloads the file for the
;                 specified date as well as the first existing file
;                 earlier in time. Do this without writing to disk.
;                 Concatenate the downloaded turret data, select only
;                 the relevant lines, and save to the turret file to
;                 be used by the pipeline.
;
;    2014-01-08 : MGL. Make downloading logs the default.
;
;    2014-01-22 : MGL. Adapt to string functions moved to the str_
;                 namespace.
;
;    2014-04-28 : MGL. An earlier name change of keyword turretfile to
;                 pathturret was apparently not complete.
;
;    2014-10-13 : MGL. Fixed bug in turret logfile dates.
;
;    2016-05-19 : THI.
;
;    2016-05-19 : MGL. Return correct r0path.
;
;    2016-07-01 : MGL. Bugfixes in the downloading and uncompressing
;                 of r0 data.
;
;    2016-08-10 : JLF. Bugfixes in turret log downloading. Now it looks
;		  backward through time properly (it was having trouble
;		  crossing year boundaries). It will not try to write
;		  output when no data exists for a date and will return 
;		  string null in pathturret when that happens.
;
;    2017-04-12 : MGL. Pig log can have isodate or dotdate in the
;                 path. Use new convertlog as DLM and not separate
;                 executable. 
;
;    2017-04-21 : MGL. Implement downloading of log files for the
;                 SHABAR, for the weather station, and for the
;                 temperature sensors.
;
;    2017-05-08 : MGL. Extended backscatter to 2017.
;
;
;-
pro red_download, date = date $
                  , overwrite = overwrite $
                  , all = all $
                  , logs = logs $
                  , shabar = shabar $
                  , pathshabar  = pathshabar  $
                  , weather = weather $
                  , pathweather  = pathweather  $
                  , temp = temp $
                  , pathtemp  = pathtemp  $
                  , pig = pig $
                  , pathpig  = pathpig  $
                  , r0 = r0 $
                  , pathr0  = pathr0  $
                  , turret = turret $
                  , pathturret = pathturret  $
                  , armap = armap $
                  , hmi = hmi $
                  , backscatter = backscatter

  any = n_elements(backscatter) gt 0 $
        or keyword_set(pig) $
        or keyword_set(armap)  $
        or keyword_set(hmi)  $
        or keyword_set(logs) $
        or keyword_set(r0)  $
        or keyword_set(shabar)  $
        or keyword_set(temp)  $
        or keyword_set(turret)  $
        or keyword_set(weather)  $
        or keyword_set(all)

  if ~any then logs = 1

  if keyword_set(all) then begin
    armap = 1
    backscatter = ['8542', '7772']
    hmi = 1
    pig = 1
    r0 = 1
    shabar = 1
    temp = 1
    turret = 1
    weather = 1
  endif

  if keyword_set(logs) then begin
    pig = 1
    r0 = 1
    shabar = 1
    temp = 1
    turret = 1
    weather = 1
  endif


  dir = 'downloads/'            ; Make this part of the crispred class structure?
  logdir = dir+'sstlogs/'

  file_mkdir, dir
  if any then file_mkdir, logdir

  if n_elements(date) gt 0 then begin
    isodate = red_strreplace(date, '.', '-', n = 2)
  endif else begin
    date = stregex(getenv('PWD'),'[12][0-9][0-9][0-9][-.][0-1][0-9][-.][0-3][0-9]',/extr)
    if date eq '' then begin
      print, 'red_download : No date given and PWD does not contain a date.'
      return
    endif
    isodate = red_strreplace(date, '.', '-', n = 2)
  endelse

  datearr = strsplit(isodate, '-', /extract)

  
  ;; Backscatter gain and psf
  if n_elements(backscatter) gt 0 then begin
    for iback = 0, n_elements(backscatter)-1 do begin
      backscatter_dir = dir+'backscatter/'
      file_mkdir, backscatter_dir
      
      downloadOK = red_geturl('http://www.isf.astro.su.se/data1/backscatter_' $
                              + backscatter[iback] + '/sst_backscatter_' $
                              + backscatter[iback] + '.tgz' $
                              , dir = backscatter_dir $
                              , overwrite = overwrite $
                              , path = backscatter_tarfile)

      if downloadOK then begin
        spawn, 'cd '+backscatter_dir+'; tar xzf '+file_basename(backscatter_tarfile)
        file_delete, backscatter_tarfile, /allow_nonexistent
        gfiles = file_search(backscatter_dir + 'cam*.backgain.' + backscatter[iback] $
                             + '_2012.f0', count = Nfiles)
        ;; Due to changes in the way the Sarnoff cameras are read out
        ;; in different years, we have to make versions of the
        ;; backscatter gain for the particular year. The files
        ;; downloaded are for 2012 (and earlier), the
        ;; backscatter_orientations matrix defined below has the value
        ;; of the parameter needed to make the rotate() command do the
        ;; needed transformation for 2013 and later.
        if datearr[0] ne '2012' then begin
          backscatter_cameras = 'cam'+['XVIII', 'XIX', 'XX', 'XXV']
          backscatter_years = '20'+['08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19']
          backscatter_orientations = bytarr(n_elements(backscatter_years) $
                                            , n_elements(backscatter_cameras))
          ;; Change orientations here if needed. Let's hope the
          ;; orientation never changed *during* an observation
          ;; season...
          backscatter_orientations[0, *] = [0, 0, 0, 0]     ; 2008 - same as 2012?
          backscatter_orientations[1, *] = [0, 0, 0, 0]     ; 2009 - same as 2012?
          backscatter_orientations[2, *] = [0, 0, 0, 0]     ; 2010 - same as 2012?
          backscatter_orientations[3, *] = [0, 0, 0, 0]     ; 2011 - same as 2012?
          backscatter_orientations[4, *] = [0, 0, 0, 0]     ; 2012
          backscatter_orientations[5, *] = [0, 0, 7, 0]     ; 2013
          backscatter_orientations[6, *] = [0, 0, 7, 0]     ; 2014 - same as 2013
          backscatter_orientations[7, *] = [0, 7, 7, 0]     ; 2015
          backscatter_orientations[8, *] = [0, 7, 7, 0]     ; 2016 - same as 2015?
          backscatter_orientations[9, *] = [0, 7, 7, 0]     ; 2017 - same as 2015?
          backscatter_orientations[10, *] = [0, 7, 7, 0]     ; 2018 - same as 2015
          backscatter_orientations[11, *] = [0, 7, 7, 0]     ; 2019 - same as 2015
          for ifile = 0, Nfiles-1 do begin
            icam = where(backscatter_cameras $
                         eq (strsplit(file_basename(gfiles[ifile]),'.', /extract))[0], Ncam)
            iyear = where(backscatter_years eq datearr[0], Nyear)
            if Ncam eq 0 or Nyear eq 0 then begin
              print, 'red::download : Backgain orientations unknown for ' $
                     + backscatter_cameras[icam] + ' in year'+datearr[0]+'.'
              print, '       Please do "git pull" in your crispred directory and try again.'
              print, '       Contact crispred maintainers if this does not help.'
              stop
            endif else begin
              ;; Read the downloaded 2012 files
              fzread, bgain, gfiles[ifile], gheader
              pfile = red_strreplace(gfiles[ifile], 'backgain', 'psf')
              fzread, bpsf, pfile, pheader
              ;; Write them with the transformation given by backscatter_orientations.
              red_loadbackscatter, backscatter_cameras[icam], isodate $
                                   , backscatter_dir, backscatter[iback] $
                                   , rotate(bgain, backscatter_orientations[iyear, icam]) $
                                   , rotate(bpsf,  backscatter_orientations[iyear, icam]) $
                                   , /write
            endelse             ; Known year and camera?
          endfor                ; ifile
        endif                   ; 2012?
      endif else begin
        print, "red_download : Couldn't download the backscatter data for " $
               + backscatter[iback]
      endelse        
    endfor                      ; iback
  endif

  
  ;; Temperature log file
  if keyword_set(temp) then begin

    tempfile = 'templog-'+strjoin(datearr, '')

    if ~file_test(logdir+tempfile) or keyword_set(overwrite) then begin
      downloadOK = red_geturl('http://www.royac.iac.es/Logfiles/temperature/' + tempfile $
                              , file = tempfile $
                              , dir = logdir $
                              , overwrite = overwrite $
                              , path = pathtemp)
      
      if ~downloadOK then begin ; also try in subfolder /{year}/
        downloadOK = red_geturl('http://www.royac.iac.es/Logfiles/temperature/' + datearr[0] $
                                + '/' + tempfile $
                                , file = tempfile $
                                , dir = logdir $
                                , overwrite = overwrite $
                                , path = pathtemp) 
      endif
      
      ;; Where did the uncompressed file end up?
      if downloadOK then pathtemp = logdir + file_basename(pathtemp)

    endif else pathtemp = logdir + tempfile
    
  endif

  
  ;; Weather log file
  if keyword_set(weather) then begin

    weatherfile = 'weather.log-' + strjoin(datearr, '')

    if ~file_test(logdir+weatherfile) or keyword_set(overwrite) then begin
      downloadOK = red_geturl('http://www.royac.iac.es/Logfiles/weather/' + weatherfile + '.xz' $
                              , file = weatherfile + '.xz' $
                              , dir = logdir $
                              , overwrite = overwrite $
                              , path = pathweather)
      
      if ~downloadOK then begin ; also try in subfolder /{year}/
        downloadOK = red_geturl('http://www.royac.iac.es/Logfiles/weather/' + datearr[0] $
                                + '/' + weatherfile + '.xz' $
                                , file = weatherfile + '.xz' $
                                , dir = logdir $
                                , overwrite = overwrite $
                                , path = pathweather) 
      endif
      
      if downloadOK then begin
        spawn, 'cd '+logdir+'; xz -df '+file_basename(pathweather)
        file_delete, pathweather, /allow_nonexistent

        ;; Where did the uncompressed file end up?
        pathweather = logdir + file_basename(pathweather,'.xz')
      endif

    endif else pathweather = logdir + weatherfile

  endif

  
  ;; SHABAR log file
  if keyword_set(shabar) then begin

    shabarfile = 'shabar.data-'+strjoin(datearr, '')

    if ~file_test(logdir+shabarfile) or keyword_set(overwrite) then begin
      downloadOK = red_geturl('http://www.royac.iac.es/Logfiles/R0/' + shabarfile + '.xz' $
                              , file = shabarfile + '.xz' $
                              , dir = logdir $
                              , overwrite = overwrite $
                              , path = pathshabar)
      
      if ~downloadOK then begin ; also try in subfolder /{year}/
        downloadOK = red_geturl('http://www.royac.iac.es/Logfiles/R0/' + datearr[0] $
                                + '/' + shabarfile + '.xz' $
                                , file = shabarfile + '.xz' $
                                , dir = logdir $
                                , overwrite = overwrite $
                                , path = pathshabar) 
      endif
      
      if downloadOK then begin
        spawn, 'cd '+logdir+'; xz -df '+file_basename(pathshabar)
        file_delete, pathshabar, /allow_nonexistent

        ;; Where did the uncompressed file end up?
        pathshabar = logdir + file_basename(pathshabar,'.xz')
      endif

    endif else pathshabar = logdir + shabarfile

  endif
  
  
  ;; R0 log file
  if keyword_set(r0) then begin

    r0file = 'r0.data.full-'+strjoin(datearr, '')

    if ~file_test(logdir+r0file) or keyword_set(overwrite) then begin
      downloadOK = red_geturl('http://www.royac.iac.es/Logfiles/R0/' + r0file + '.xz' $
                              , file = r0file + '.xz' $
                              , dir = logdir $
                              , overwrite = overwrite $
                              , path = pathr0)
      
      if ~downloadOK then begin ; also try in subfolder /{year}/
        downloadOK = red_geturl('http://www.royac.iac.es/Logfiles/R0/' + datearr[0] $
                                + '/' + r0file + '.xz' $
                                , file = r0file + '.xz' $
                                , dir = logdir $
                                , overwrite = overwrite $
                                , path = pathr0) 
      endif
      
      if downloadOK then begin
        spawn, 'cd '+logdir+'; xz -df '+file_basename(pathr0)
        file_delete, pathr0, /allow_nonexistent

        ;; Where did the uncompressed file end up?
        pathr0 = logdir + file_basename(pathr0,'.xz')
      endif

    endif else pathr0 = logdir + r0file

  endif

  
  ;; PIG log file
  if keyword_set(pig) then begin
    pigfile = 'rmslog_guidercams'

    DownloadOK = red_geturl('http://www.royac.iac.es/Logfiles/PIG/' + isodate $
                            + '/' + pigfile $
                            , file = pigfile+'_'+isodate $
                            , dir = logdir $
                            , overwrite = overwrite $
                            , path = pathpig)

    if ~DownloadOK then begin
      dotdate = strjoin(datearr, '.')
      DownloadOK = red_geturl('http://www.royac.iac.es/Logfiles/PIG/' + dotdate $
                              + '/' + pigfile $
                              , file = pigfile+'_'+isodate $
                              , dir = logdir $
                              , overwrite = overwrite $
                              , path = pathpig)     
    endif

    if DownloadOK then begin
      ;; We actually want the logfile converted to time and x/y
      ;; coordinates (in arcseconds).
      pathpig += '_converted'
      if ~file_test(pathpig) then begin
        print, 'red_download : Converting PIG log file...'
        rdx_convertlog, logdir+pigfile+'_'+isodate, logdir+pigfile+'_'+isodate+'_converted' $
                        , dx=31.92, dy=14.81 $
                        , rotation=84.87, scale=4.935, average=16
      endif else begin
        print, 'red_download : Converted PIG log file already exists.'
      endelse
    endif else begin
      ;; We tried to download but failed. So any existing files may
      ;; be corrupt or not correspond to the current state.
      if pathpig ne '' then begin
        file_delete, pathpig, /allow_nonexistent
        file_delete, pathpig + '_' + isodate + '_converted', /allow_nonexistent
        pathpig = ''
      endif
    endelse
  endif

  
  ;; Turret log file

  if keyword_set(turret) then begin

    ;; Turret log data for a particular day can actually be in the
    ;; turret log file of an earlier day. So we need to search days
    ;; backwards until we find one. Then we should concatenate the
    ;; two files and filter the result to get rid of data for another
    ;; days and header info that are interspersed with the data.

    ;; The name of the concatenated and filtered file
    pathturret = logdir+'positionLog_'+red_strreplace(isodate, '-', '.', n = 2)+'_final'
    
    if ~file_test(pathturret) or keyword_set(overwrite) then begin
      
      ;; First try the particular date:
      
      pathturret1 = 'positionLog_'+red_strreplace(isodate, '-', '.', n = 2)
      OK1 = red_geturl('http://www.royac.iac.es/Logfiles/turret/' $
                       + datearr[0]+'/'+pathturret1 $
                       , contents = contents1 $
;                         , dir = logdir $
                       , /overwrite $
                      ) 

      ;; Try previous days until one is found
      predatearr = datearr
      repeat begin
        
        ;; Make isodate for one day earlier
        caldat,julday(predatearr[1],predatearr[2],predatearr[0])-1,month,day,year
        predatearr = [string(year,format='(i4)'),$
                      string(month,format='(i02)'),$
                      string(day,format='(i02)')] 
;            predatearr[2] = string(predatearr[2]-1, format = '(i02)')
;            if predatearr[2] eq 0 then begin
;               predatearr[2] = '31'
;               predatearr[1] = string(predatearr[1]-1, format = '(i02)')
;               if predatearr[1] eq 0 then begin
;                  predatearr[2] = '12'
;                  predatearr[0] += string(predatearr[0]-1, format = '(i04)')
;               endif
;            endif
        preisodate = strjoin(predatearr, '-')
        
        ;; Try to download
        pathturret2 = 'positionLog_'+red_strreplace(preisodate, '-', '.', n = 2)
        print, 'Try '+pathturret2
        OK2 = red_geturl('http://www.royac.iac.es/Logfiles/turret/' $
                         + predatearr[0]+'/'+pathturret2 $
                         , contents = contents2 $
;                            , dir = logdir $
                         , /overwrite $
                        ) 
        
      endrep until OK2

      ;; Concatenate the downloaded files (if needed)
      if OK1 then contents = [contents2, contents1] else contents = contents2
      ;; Filter on date
      turretdate = strjoin(datearr, '/')
      idx = where(strmatch(contents, turretdate+'*'),cnt)
      if cnt ne 0 then begin
        contents = contents[idx]
        ;; Filter on line type, don't want the RA/Decl lines
        contents = contents(where(~strmatch(contents, '*h*m*')))

        ;; Write to disk
        openw, wlun, /get_lun, pathturret
        printf, wlun, contents, format = '(a0)'
        free_lun, wlun
      endif else begin
        message,'No log data for the requested date.',/info
        pathturret = ''         ; don't return a path, we didn't write anything
      endelse
;     grepdate = strjoin(datearr, '/')
;     cmd = 'cd downloads/sstlogs ;'       ; Go to download directory
;     cmd += ' cat positionLog_????.??.??' ; Concatenate the turret files
;     cmd += ' | grep '+grepdate           ; Only lines from the relevant day
;     cmd += ' | grep -v h '               ; Remove RA/Decl lines
;     cmd += ' > positionLog'              ; Output to file to be used
;     print, cmd
      
    endif                       ; non-existent or overwrite
    
  endif                         ; keyword_set(turret)
  
  
  ;; Active regions map
  arfile = strjoin(datearr, '')+'.1632_armap.png'
  if keyword_set(armap) then begin
    tmp = red_geturl('http://kopiko.ifa.hawaii.edu/ARMaps/Archive/' + datearr[0] $
                     + '/' + arfile $
                     , file = arfile, dir = dir, overwrite = overwrite) 
  endif

  ;; HMI images and movies
  if keyword_set(hmi) then begin
    hmidir = '/data/hmi/images/'+strjoin(datearr, '/')+'/'
    hmisite = 'http://jsoc.stanford.edu'
    hmitimestamps = ['08', '10', '12', '14', '16', '18']+'0000'
    hmitypes = ['Ic_flat_4k', 'M_color_4k']
    for i = 0, n_elements(hmitimestamps)-1 do begin
      for j = 0, n_elements(hmitypes)-1 do begin
        hmifile = strjoin(datearr, '')+'_'+hmitimestamps[i]+'_'+hmitypes[j]+'.jpg'
        tmp = red_geturl(hmisite+hmidir+hmifile $ 
                         , file = hmifile, dir = dir+'/HMI/', overwrite = overwrite) 
      endfor
    endfor                      ; i
    hmimovies = ['Ic_flat_2d', 'M_2d', 'M_color_2d']+'.mpg'
    for i = 0, n_elements(hmimovies)-1 do begin
      hmifile = hmimovies[i]
      tmp = red_geturl(hmisite+hmidir+hmifile $
                       , file = hmifile, dir = dir+'/HMI/', overwrite = overwrite) 
    endfor                      ; i
  endif
  
end
