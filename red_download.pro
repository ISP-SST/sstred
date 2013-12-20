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
;
;
;
;-
pro red_download, date = date $
                  , overwrite = overwrite $
                  , all = all $
                  , pig = pig $
                  , pathpig  = pathpig  $
                  , r0 = r0 $
                  , pathr0  = pathr0  $
                  , pathturret = pathturret  $
                  , turret = turret $
                  , armap = armap $
                  , hmi = hmi

  if keyword_set(all) then begin
     pig = 1
     r0 = 1
     turret = 1
     armap = 1
     hmi = 1
  endif

  any = keyword_set(pig) or keyword_set(turret) or keyword_set(armap) or keyword_set(hmi) or keyword_set(r0)

  dir = 'downloads/'            ; Make this part of the crispred class structure?
  logdir = dir+'sstlogs/'

  if any then file_mkdir, logdir

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

  datearr = strsplit(isodate, '-', /extract)

  ;; R0 log file
  if keyword_set(r0) then begin
     r0file = 'r0.data.full-'+strjoin(datearr, '')
     downloadOK = red_geturl('http://www.royac.iac.es/Logfiles/R0/' + r0file $
                             , file = r0file $
                             , dir = logdir $
                             , overwrite = overwrite $
                             , path = pathr0) 
  endif

  ;; PIG log file
  if keyword_set(pig) then begin
     pigfile = 'rmslog_guidercams'
     DownloadOK = red_geturl('http://www.royac.iac.es/Logfiles/PIG/' + isodate + '/' + pigfile $
                             , file = pigfile $
                             , dir = logdir $
                             , overwrite = overwrite $
                             , path = pathpig)

     if DownloadOK then begin
        ;; We actually want the logfile converted to time and x/y
        ;; coordinates (in arcseconds).
        pathpig += '_'+isodate+'_converted'
        if ~file_test(pathpig) then begin
           pig_N = 16           ; # of positions to average when converting. Originally ~16 pos/s.
           convertcmd = 'cd '+logdir+'; convertlog --dx 31.92 --dy 14.81' $
                        + ' --rotation 84.87 --scale 4.935 '
           if pig_N gt 1 then convertcmd += '-a ' + strtrim(pig_N, 2) + ' '
           print, 'red_download : Converting PIG log file...'
           spawn, convertcmd+' '+pigfile+' > '+pigfile+'_'+isodate+'_converted'
;        file_link, logdir+pigfile+'_'+isodate+'_converted', 'log_pig'
;        print, 'red_download : Linked to ' + link
        endif else begin
           print, 'red_download : Converted PIG log file already exists.'
        endelse
     endif else begin
        ;; We tried to download but failed. So any existing files may
        ;; be corrupt or not correspond to the current state.
        file_delete, pathpig, /allow_nonexistent
        file_delete, pathpig + '_' + isodate + '_converted', /allow_nonexistent
        pathpig = ''
     endelse
  endif

  ;; Turret log file

  ;; The turret log data for a particular day can actually be in the
  ;; turret log file of an earlier day. So if todays file does not
  ;; exist we should search days backwards until we find one. Then we
  ;; should filter it to get rid of data for aother days and header
  ;; info that are interspersed with the data. Not implemented yet.
  if keyword_set(turret) then begin
     turretfile = 'positionLog_'+strreplace(isodate, '-', '.', n = 2)
     downloadOK = red_geturl('http://www.royac.iac.es/Logfiles/turret/' + datearr[0]+'/'+turretfile $
                             , file = turretfile $
                             , dir = logdir $
                             , overwrite = overwrite $
                             , path = pathturret) 
  endif
  
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
