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
;      Set this to download without checking if the file already exists.
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
;    2014-01-08 : MGL. Turn it into a method in order to access
;                 members, like the date and directory names. Remove
;                 the date keyword.
;
;    2014-01-10 : MGL. Do not try to remove files if the file name is
;                 the empty string.
;
;    2014-01-22 : MGL. Adapt to string functions moved to the str_
;                 namespace.
;
;    2015-09-25 : MGL. Don't download the PIG logfile unless
;                 asked for.
;
;    2016-02-15 : MGL. Added downloading of backscatter data files for
;                 8542 and 7772 prefilters.
;
;    2016-02-17 : MGL. Use file name provided by red::loadbackscatter
;                 rather than constructiong one.
;
;    2016-02-19 : MGL. Transform backscatter psf the same way as the
;                 backscatter gain.
;
;    2016-02-24 : MGL. Make backscatter renaming work also for years
;                 prior to 2012.
;
;    2017-04-12 : MGL. Pig log can have isodate or dotdate in the
;                 path. Use new convertlog as DLM and not separate
;                 executable.
;
;    2017-08-18 : THI. Also download the FPI calibration values
;                 (linedef.py). Include FPI calibration in implied
;                 downloads.
;
;    2021-08-27 : MGL. Use red_download_log for r0 and PIG log files.
;                 Include the PIG log files in the logs category.
;
;-
pro red::download, overwrite = overwrite $
                   , all = all $
                   , logs = logs $
                   , pig = pig $
                   , pathpig  = pathpig  $
                   , linedefs = linedefs $
                   , r0 = r0 $
                   , pathr0  = pathr0  $
                   , pathturret = pathturret  $
                   , turret = turret $
                   , armap = armap $
                   , hmi = hmi $
                   , backscatter = backscatter
  
  any = keyword_set(pig) $
        or keyword_set(turret)  $
        or keyword_set(armap)  $
        or keyword_set(hmi)  $
        or keyword_set(linedefs)  $
        or keyword_set(r0)  $
        or keyword_set(all)  $
        or keyword_set(logs) $
        or n_elements(backscatter) gt 0

  if ~any then begin 
    logs = 1
    if self.isodate lt red_dates(tag = 'CHROMIS tuning metadata') then linedefs = 1
  endif
  
  if keyword_set(all) then begin
    pig = 1
    if self.isodate lt red_dates(tag = 'CHROMIS tuning metadata') then linedefs = 1
    r0 = 1
    turret = 1
    armap = 1
    hmi = 1
    backscatter = ['8542', '7772']
  endif

  if keyword_set(logs) then begin
    r0 = 1
    turret = 1
    pig = 1
  endif
  
;  if self.isodate ge red_dates(tag = 'CHROMIS tuning metadata') then linedefs = 0


  dir = 'downloads/'            ; Make this part of the crispred class structure?
  
  if any then file_mkdir, self.log_dir

;  if n_elements(date) gt 0 then begin
;     isodate = red_strreplace(date, '.', '-', n = 2)
;  endif else begin
;     date = stregex(getenv('PWD'),'[12][0-9][0-9][0-9][-.][0-1][0-9][-.][0-3][0-9]',/extr)
;     if date eq '' then begin
;        print, 'red_download : No date given and PWD does not contain a date.'
;        return
;     endif
;     isodate = red_strreplace(date, '.', '-', n = 2)
;  endelse

  datearr = strsplit(self.isodate, '-', /extract)


  ;; Backscatter gain and psf
  if n_elements(backscatter) gt 0 then begin

    file_mkdir, self.descatter_dir

    for iback = 0, n_elements(backscatter)-1 do begin
      
      downloadOK = red_geturl('http://www.isf.astro.su.se/data1/backscatter_' $
                              + backscatter[iback] + '/sst_backscatter_' $
                              + backscatter[iback] + '.tgz' $
                              , dir = self.descatter_dir $
                              , overwrite = overwrite $
                              , path = backscatter_tarfile)

      if downloadOK then begin
        spawn, 'cd '+self.descatter_dir+'; tar xzf '+file_basename(backscatter_tarfile)
        file_delete, backscatter_tarfile, /allow_nonexistent
        gfiles = file_search(self.descatter_dir+'cam*.backgain.'+ backscatter[iback] + '_2012.f0' $
                             , count = Nfiles)
        ;; Due to changes in the way the Sarnoff cameras are read
        ;; out in different years, we have to make versions of the
        ;; backscatter gain for the particular year. The files
        ;; downloaded are for 2012 (and earlier), the
        ;; backscatter_orientations matrix defined below has the
        ;; value of the parameter needed to make the rotate()
        ;; command do the needed transformation for 2013 and later.
        if datearr[0] ne '2012' then begin
          backscatter_cameras = 'cam'+['XVIII', 'XIX', 'XX', 'XXV']
          backscatter_years = '20'+['08', '09', '10', '11', '12', '13', '14', '15' $
                                    , '16', '17', '18', '19', '20', '21', '22']
          backscatter_orientations = bytarr(n_elements(backscatter_years) $
                                            , n_elements(backscatter_cameras))
          ;; Change orientations here if needed. Let's hope the
          ;; orientation never changed *during* an observation
          ;; season...
          backscatter_orientations[ 0, *] = [0, 0, 0, 0]  ; 2008 - same as 2012?
          backscatter_orientations[ 1, *] = [0, 0, 0, 0]  ; 2009 - same as 2012?
          backscatter_orientations[ 2, *] = [0, 0, 0, 0]  ; 2010 - same as 2012?
          backscatter_orientations[ 3, *] = [0, 0, 0, 0]  ; 2011 - same as 2012?
          backscatter_orientations[ 4, *] = [0, 0, 0, 0]  ; 2012
          backscatter_orientations[ 5, *] = [0, 0, 7, 0]  ; 2013
          backscatter_orientations[ 6, *] = [0, 0, 7, 0]  ; 2014 - same as 2013
          backscatter_orientations[ 7, *] = [0, 7, 7, 0]  ; 2015
          backscatter_orientations[ 8, *] = [0, 7, 7, 0]  ; 2016 - same as 2015?
          backscatter_orientations[ 9, *] = [0, 7, 7, 0]  ; 2017 - same as 2015?
          backscatter_orientations[10, *] = [0, 7, 7, 0]  ; 2018 - same as 2015?
          backscatter_orientations[11, *] = [0, 7, 7, 0]  ; 2019 - same as 2015
          backscatter_orientations[12, *] = [0, 7, 7, 0]  ; 2020 - same as 2015
          backscatter_orientations[13, *] = [0, 7, 7, 0]  ; 2021 - same as 2015              
          backscatter_orientations[14, *] = [0, 7, 7, 0]  ; 2022 - same as 2015              
          for ifile = 0, Nfiles-1 do begin
            icam = where(backscatter_cameras $
                         eq (strsplit(file_basename(gfiles[ifile]),'.', /extract))[0], Ncam)
            iyear = where(backscatter_years eq datearr[0], Nyear)
            if Ncam eq 0 or Nyear eq 0 then begin
              print, 'red::download : Backgain orientations unknown for ' + backscatter_cameras[icam] $
                     + ' in year'+datearr[0]+'.'
              print, '               Please do "git pull" in your crispred directory and try again.'
              print, '               Contact crispred maintainers if this does not help.'
              stop
            endif else begin
              ;; Read the downloaded 2012 files
              fzread, bgain, gfiles[ifile], gheader
              pfile = red_strreplace(gfiles[ifile], 'backgain', 'psf')
              fzread, bpsf, pfile, pheader
              ;; Write them with the transformation given by backscatter_orientations.
;                    self -> loadbackscatter, backscatter_cameras[icam], backscatter[iback] $
;                                             , rotate(bgain, backscatter_orientations[iyear, icam]) $
;                                             , rotate(bpsf,  backscatter_orientations[iyear, icam]) $
;                                             , /write
;
;
              rogain = rotate(bgain, backscatter_orientations[iyear, icam])
              ropsf = rotate(bpsf, backscatter_orientations[iyear, icam])
              self -> loadbackscatter, backscatter_cameras[icam], backscatter[iback] $
                                       , rogain $
                                       , ropsf $
                                       , /write
            endelse             ; Known year and camera?
          endfor                ; ifile
        endif                   ; 2012?
      endif else begin
        print, "red::download : Couldn't download the backscatter data for " + backscatter[iback]
      endelse        
    endfor                      ; iback
  endif
    
  ;; linedefs
  if keyword_set(linedefs) && (min(where(*self.cameras eq 'Chromis-N')) ne -1) then begin
    dotdate = strjoin(datearr, '.')
    ldpath = 'http://www.sst.iac.es/Logfiles/CHROMIS/linedef/'
    downloadOK = red_geturl(ldpath , contents=ldfiles )
    if downloadOK then begin
      todays_linedefs = ldfiles[where(strmatch(ldfiles, 'linedef.py-'+dotdate+'*', /FOLD_CASE) EQ 1)]
      todays_linedef_times = todays_linedefs
      for i=0,n_elements(todays_linedef_times)-1 do begin
        todays_linedef_times[i] = red_time_conv( strmid(todays_linedef_times[i],strlen(todays_linedef_times[i])-8,8) )
      endfor
      ntlt = n_elements( todays_linedef_times )
      if ntlt gt 0 then begin 
        if ptr_valid(self.dark_dir)  then red_append,all_dirs,*self.dark_dir ; list all data directories of the day
        if ptr_valid(self.flat_dir)  then red_append,all_dirs,*self.flat_dir 
        if ptr_valid(self.pinh_dirs) then red_append,all_dirs,*self.pinh_dirs
        if ptr_valid(self.data_dirs) then red_append,all_dirs,*self.data_dirs
        idx = where(strmatch(all_dirs, '*CHROMIS-calib*', /FOLD_CASE) EQ 1, compl=cidx)
        all_dirs = all_dirs[cidx] ; skip calib-directories (often saved/copied by accident)
        for i=0,n_elements(all_dirs)-1 do begin
          red_append, all_times, red_Time_conv( strmid(all_dirs[i],strlen(all_dirs[i])-8,8) )
        endfor
        data_times = all_times[uniq(all_times, sort(all_times))]
        for i=0,n_elements(data_times)-1 do begin
          lastidx = max(where( todays_linedef_times lt data_times[i] ))
          if lastidx lt 0 then begin ; TODO: if no calib exists, get one from previous day?
            print, 'red::download : There is not calibration done before the data: ' + $
                   all_dirs[ where(all_times eq data_times[i]) ]
            print, 'red::download : It is safe to ignore this warning if this dataset will not be used for science ' + $
                   ' (it might be an early-morning test or calibration).'
            continue
          endif else begin
            red_append, used_linedefs, todays_linedef_times[lastidx]
          endelse
        endfor
        used_linedefs = used_linedefs[uniq(used_linedefs, sort(used_linedefs))]
        todays_linedefs = todays_linedefs[ where(todays_linedef_times eq used_linedefs) ]
        for i=0,n_elements(todays_linedefs)-1 do begin
          downloadOK = red_geturl( ldpath+todays_linedefs[i], dir=dir, overwrite = overwrite )
        endfor
      endif
    endif
  endif

  ;; R0 log file
  if keyword_set(r0) then begin
    red_download_log, 'r0', self.isodate, self.log_dir $
                      , localpath = pathr0 $                      
                      , status = status
  endif

  ;; PIG log file
  if keyword_set(pig) then begin
    red_download_log, 'pig', self.isodate, self.log_dir $
                      , localpath = pathpig $                      
                      , status = status
  endif
  
  ;; Turret log file

  if keyword_set(turret) then begin

    ;; Turret log data for a particular day can actually be in the
    ;; turret log file of an earlier day. So we need to search days
    ;; backwards until we find one. Then we should concatenate the two
    ;; files and filter the result to get rid of data for another days
    ;; and header info that are interspersed with the data.

    ;; The name of the concatenated and filtered file
;     turretfile = self.log_dir+'positionLog_'+red_strreplace(self.isodate, '-', '.', n = 2)+'_final'
    if ~file_test(self.telog) or keyword_set(overwrite) then begin
      
      ;; First try the particular date:
      
      turretfile1 = 'positionLog_'+red_strreplace(self.isodate, '-', '.', n = 2)
      OK1 = red_geturl('http://www.sst.iac.es/Logfiles/turret/' $
                       + datearr[0]+'/'+turretfile1 $
                       , contents = contents1 $
;                         , dir = self.log_dir $
                       , /overwrite $
                      ) 

      ;; Try previous days until one is found
      predatearr = datearr
      repeat begin
        
        ;; Make isodate for one day earlier
        predatearr[2] = string(predatearr[2]-1, format = '(i02)')
        if predatearr[2] eq 0 then begin
          predatearr[2] = '31'
          predatearr[1] = string(predatearr[1]-1, format = '(i02)')
          if predatearr[1] eq 0 then begin
            predatearr[2] = '12'
            predatearr[0] += string(predatearr[0]-1, format = '(i04)')
          endif
        endif
        preisodate = strjoin(predatearr, '-')
        
        ;; Try to download
        turretfile2 = 'positionLog_'+red_strreplace(preisodate, '-', '.', n = 2)
        print, 'Try '+turretfile2
        OK2 = red_geturl('http://www.sst.iac.es/Logfiles/turret/' $
                         + datearr[0]+'/'+turretfile2 $
                         , contents = contents2 $
;                            , dir = self.log_dir $
                         , /overwrite $
                        ) 
        
      endrep until OK2

      ;; Concatenate the downloaded files (if needed)
      if OK1 then contents = [contents2, contents1] else contents = contents2
      ;; Filter on date
      turretdate = strjoin(datearr, '/')
      contents = contents(where(strmatch(contents, turretdate+'*')))
      ;; Filter on line type, don't want the RA/Decl lines
      ;; from when the telescope is parked
      contents = contents(where(~strmatch(contents, '*h*m*')))

      ;; Write to disk
      openw, wlun, /get_lun, self.telog
      printf, wlun, contents, format = '(a0)'
      free_lun, wlun
      
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
    tmp = red_geturl('http://kopiko.ifa.hawaii.edu/ARMaps/Archive/' + datearr[0]+'/'+arfile $
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
    endfor
    hmimovies = ['Ic_flat_2d', 'M_2d', 'M_color_2d']+'.mpg'
    for i = 0, n_elements(hmimovies)-1 do begin
      hmifile = hmimovies[i]
      tmp = red_geturl(hmisite+hmidir+hmifile $
                       , file = hmifile, dir = dir+'/HMI/', overwrite = overwrite) 
    endfor
  endif
  
end
