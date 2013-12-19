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
;    pig : in, optional, type=boolean
;
;      Try to download the SST/PIG log file.
;
;    r0 : in, optional, type=boolean
;
;      Try to download the SST/AO r0 log file.
;
;    turret : in, optional, type=boolean
;
;      Try to download the SST/Turret log file. 
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
; :History:
; 
; 
; 
;
;
;
;-
pro red_download, date = date, overwrite = overwrite, all = all, pig = pig, r0 = r0, turret = turret, armap = armap, hmi = hmi

  if keyword_set(all) then begin
     pig = 1
     r0 = 1
     turret = 1
     armap = 1
     hmi = 1
  endif

  any = keyword_set(pig) or keyword_set(turret) or keyword_set(armap) or keyword_set(hmi) or keyword_set(r0)

  dir = 'downloads/'            ; Make this part of the crispred class structure?

  if any then file_mkdir, dir

  if n_elements(date) gt 0 then begin
     isodate = strreplace(date, '.', '-', n = 2)
  endif else begin
     date = stregex(getenv('PWD'),'[12][0-9][0-9][0-9][-.][0-1][0-9][-.][0-3][0-9]',/extr)
     if date eq '' then begin
        print, 'red_download : No date given and PWD does not contain a date.'
        return
     endif
     isodate = strreplace(date, '.', '-', n = 2)
  endelse

  datearr = strsplit(isodate, '-', /extr)

  ;; R0 log file
  if keyword_set(r0) then begin
     r0file = 'r0.data.full-'+strjoin(datearr, '')
     if file_test(dir+r0file) and ~keyword_set(overwrite) then begin
        print, 'red_download : R0 log file already downloaded.'
     endif else begin
        print, 'red_download : Trying to download R0 log file.'
        url = 'http://www.royac.iac.es/Logfiles/R0/' + r0file
        print, url
        if red_geturl(url, file = r0file, dir = dir) then begin
           print, 'Downloaded OK.'
        endif else begin
           print, 'Download failed.'
        endelse
     endelse
  endif

  ;; PIG log file
  pigfile = 'rmslog_guidercams'
  if keyword_set(pig) then begin
     if file_test(dir+pigfile) and ~keyword_set(overwrite) then begin
        print, 'red_download : PIG log file already downloaded.'
     endif else begin
        print, 'red_download : Trying to download PIG log file.'
        url = 'http://www.royac.iac.es/Logfiles/PIG/' + isodate + '/' + pigfile
        print, url
        if red_geturl(url, file = pigfile, dir = dir) then begin
           print, 'Downloaded OK.'
        endif else begin
           print, 'Download failed.'
        endelse
        ;; Convert the logfile to time and x/y coordinates (in
        ;; arcseconds).
        pig_N = 16              ; # of positions to average when converting. Originally ~16 pos/s.
        convertcmd = 'cd '+dir+'; convertlog --dx 31.92 --dy 14.81' $
                     + ' --rotation 84.87 --scale 4.935 '
        if pig_N gt 1 then convertcmd += '-a ' + strtrim(pig_N, 2) + ' '
        print, 'Converting PIG log file...'
        spawn, convertcmd+' rmslog_guidercams > '+pigfile+'_'+isodate+'_converted'
     endelse
  endif

  ;; Turret log file
  turretfile = 'positionLog_'+strreplace(isodate, '-', '.', n = 2)
  if keyword_set(turret) then begin
     if file_test(dir+turretfile) and ~keyword_set(overwrite) then begin
        print, 'red_download : Turret log file already downloaded.'
     endif else begin
        print, 'red_download : Trying to download Turret log file.'
        url = 'http://www.royac.iac.es/Logfiles/turret/' + datearr[0]+'/'+turretfile
        print, url
        if red_geturl(url, file = turretfile, dir = dir) then begin
           print, 'Downloaded OK.'
        endif else begin
           print, 'Download failed.'
        endelse
     endelse
  endif
  
  ;; Active regions map
  arfile = strjoin(datearr, '')+'.1632_armap.png'
  if keyword_set(armap) then begin
     if file_test(dir+arfile) and ~keyword_set(overwrite) then begin
        print, 'red_download : AR map already downloaded.'
     endif else begin
        print, 'red_download : Trying to download AR map.'
        url = 'http://kopiko.ifa.hawaii.edu/ARMaps/Archive/' + datearr[0]+'/'+arfile
        print, url
        if red_geturl(url, file = arfile, dir = dir) then begin
           print, 'Downloaded OK.'
        endif else begin
           print, 'Download failed.'
        endelse
     endelse
  endif

  ;; HMI images and movies
;  months = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov', 'dec']
;  hmidate = datearr[2]+months[long(datearr[1])-1]+strmid(datearr[0], 2, 2)
  if keyword_set(hmi) then begin
     hmidir = '/data/hmi/images/'+strjoin(datearr, '/')+'/'
     hmisite = 'http://jsoc.stanford.edu'
     hmitimestamps = ['08', '10', '12', '14', '16', '18']+'0000'
     hmitypes = ['Ic_flat_4k', 'M_color_4k']
     for i = 0, n_elements(hmitimestamps)-1 do begin
        for j = 0, n_elements(hmitypes)-1 do begin
           hmifile = strjoin(datearr, '')+'_'+hmitimestamps[i]+'_'+hmitypes[j]+'.jpg'
          if file_test(dir+'/HMI/'+hmifile) and ~keyword_set(overwrite) then begin
              print, 'red_download : SDO/HMI image '+hmifile+' already downloaded.'
           endif else begin
              print, 'red_download : Trying to download SDO/HMI image '+hmifile+'.'
              url = hmisite+hmidir+hmifile
              print, url
              if red_geturl(url, file = hmifile, dir = dir+'/HMI/') then begin
                 print, 'Downloaded OK.'
              endif else begin
                 print, 'Download failed.'
              endelse
           endelse
        endfor
     endfor
     hmimovies = ['Ic_flat_2d', 'M_2d', 'M_color_2d']+'.mpg'
     for i = 0, n_elements(hmimovies)-1 do begin
        hmifile = hmimovies[i]
        if file_test(dir+'/HMI/'+hmifile) and ~keyword_set(overwrite) then begin
           print, 'red_download : SDO/HMI movie '+hmifile+' already downloaded.'
        endif else begin
           print, 'red_download : Trying to download SDO/HMI movie '+hmifile+'.'
           url = hmisite+hmidir+hmifile
           print, url
           if red_geturl(url, file = hmifile, dir = dir+'/HMI/') then begin
              print, 'Downloaded OK.'
           endif else begin
              print, 'Download failed.'
           endelse
        endelse
     endfor
  endif

end
