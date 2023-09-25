; docformat = 'rst'

;+
; Download log files for PIG and r0, possibly with uncompression and
; other processing.
; 
; Compressed versions of the log files are preferred.
; 
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Params:
; 
;    instrument : in, type=string
; 
;      The instrument to download the log file for, 'PIG' or 'R0'.
; 
;    date : in, type=string
; 
;      The date in ISO format.
;
;    outdir : in, type=string 
; 
;      The directory in which to write the log file.
; 
; 
; :Keywords:
; 
;    localpath : out, optional, type=string
; 
;      The path to the downloaded and processed log file.
; 
;    overwrite : in, optional, type=booldean
; 
;      Set this to overwrite a previously downloaded file. If a
;      compressed version has not been downloaded before, an existing
;      local file will be overwritten without /overwrite.
; 
;    status : out, optional, type=integer
; 
;      Zero for success, non-zero when problems.
; 
; :History:
; 
;    2021-08-27 : MGL. First version.
; 
;    2023-08-14 : MGL. Change parameters for rdx_convertlog.
; 
;-
pro red_download_log, instrument, date, outdir $
                      , localpath = localpath $
                      , overwrite = overwrite $
                      , status = status

  if n_elements(instrument) eq 0 then return  
  if n_elements(date) eq 0 then return  
  if n_elements(outdir) eq 0 then return  
  
  datearr = strsplit(date, '-.', /extract)
  isodate = strjoin(datearr, '-')  
  dotdate = strjoin(datearr, '.')  

  file_mkdir, outdir
  
  case strupcase(instrument) of

    'PIG' : begin
      remotedirs = 'http://www.sst.iac.es/Logfiles/PIG/' + [dotdate, isodate] + '/' 
      remotefile = 'rmslog_guidercams'      
      localfile  = 'rmslog_guidercams' + '_' + isodate + '_converted'
    end

    'R0' : begin
      remotedirs = 'http://www.sst.iac.es/Logfiles/R0/' + ['', datearr[0]+'/']
      remotefile = 'r0.data.full-'+strjoin(datearr, '')      
      localfile  = 'r0.data.full-'+strjoin(datearr, '')      
    end

    'WFWFS' : begin
      remotedirs = 'http://www.sst.iac.es/Logfiles/wfwfs/' 
      remotefile = 'wfwfs_r0.log-'+strjoin(datearr, '')
      localfile  = remotefile
    end

    else : begin
      print, 'Instrument '+instrument+' not implemented.'
      stop
    end
    
  endcase

  downloadpath = outdir + '/' + remotefile          
  localpath    = outdir + '/' + localfile          
  
  if ~keyword_set(overwrite) then begin
    if file_test(downloadpath+'.xz') then begin
      print, 'Compressed file exists'

      if file_test(downloadpath) then begin
        print, 'Uncompressed file also exists'
      endif else begin
        spawn, 'cd '+file_dirname(localpath) $
               + '; xzcat '+file_basename(downloadpath)+ '.xz > ' + file_basename(downloadpath)
      endelse
      
      status = 0
      return
    endif
  endif
    
  ;; Try first for an xz compressed version

  for ipath = 0, n_elements(remotedirs)-1 do begin
    print, ipath
    downloadOK = red_geturl(remotedirs[ipath] + remotefile+'.xz' $
                            , file = remotefile+'.xz' $
                            , dir = outdir $
                            , overwrite = overwrite $
                            , path = actual_downloadpath)
    if downloadOK then break
  endfor                        ; ipath

  if downloadOK then begin
    print, 'uncompress'
    spawn, 'cd '+file_dirname(localpath) $
           + '; xzcat '+file_basename(downloadpath)+ '.xz > ' + file_basename(downloadpath)
    
  endif else begin

    ;; Look for a file without compression
    
    for ipath = 0, n_elements(remotedirs)-1 do begin
      print, ipath
      downloadOK = red_geturl(remotedirs[ipath] + remotefile $
                              , file = remotefile $
                              , dir = outdir $
                              , overwrite = overwrite $
                              , path = actual_downloadpath)
      if downloadOK then break
    endfor                      ; ipath

  endelse

  ;; Any further processing needed?

  case strupcase(instrument) of

    'PIG' : begin
      if ~file_test(localfile) then begin
        
        ;; PIG logfile needs to be converted. For that we need the
        ;; limb calibration results.
        print, 'red_download_log : Downloading 4-limb calibration file...'
        downloadOK = red_geturl('http://www.sst.iac.es/Logfiles/PIG/PIG_center' $
                                , file = 'PIG_center' $
                                , dir = outdir $
                                , overwrite = overwrite $
                                , path = pigcfile)
        pigctemplate = { VERSION: 1.0, $
                         DATASTART: 0L, $
                         DELIMITER: 32B, $
                         MISSINGVALUE: 0.0, $
                         COMMENTSYMBOL: "#", $
                         FIELDCOUNT: 4, $
                         FIELDTYPES: [7L, 5L, 5L, 5L], $
                         FIELDNAMES: [ "date", "dx", "dy", "scale"], $
                         FIELDLOCATIONS: [0L, 12L, 31L, 50L], $
                         FIELDGROUPS: [0L, 1L, 2L, 3L] $
                       }
;        pigcstruct = {date:'', dx:0.d, dy:0.d, scale:0.d}
        pig_center = read_ascii(pigcfile, TEMPLATE=pigctemplate, count = count)

        ;; In theory we'd use the dx,dy,scale parameters from the
        ;; calibration that is nearest in time. But there are outliers
        ;; so we better remove them by use of median filtering. Also,
        ;; there are differences from year to year (sometimes large)
        ;; so we first select just the year.
        indx = where(strmid(pig_center.date, 0, 4) eq datearr[0], Nwhere)
        if Nwhere eq 0 then begin
          ;; Average first half of 2023
          dx=20.24
          dy=18.35
          scale=4.31830
        endif else begin
          pig_center_date = pig_center.date[indx]
          pig_center_dx = median(pig_center.dx[indx], 5)
          pig_center_dy = median(pig_center.dy[indx], 5)
          pig_center_scale = median(pig_center.scale[indx], 5)
          juldate = date_conv(isodate, 'J')
          pig_center_juldates = dblarr(Nwhere)
          for i = 0, Nwhere-1 do pig_center_juldates[i] = date_conv(pig_center_date[i], 'J') ; Convert to Julian
          dist = min(abs(pig_center_juldates - juldate), ii)                                 ; Nearest date ii
          dx = pig_center_dx[ii]
          dy = pig_center_dy[ii]
          scale = pig_center_scale[ii]
        endelse
        
        ;; Now do the converting
        print, 'red_download_log : Converting PIG log file...'
        if file_test(downloadpath) then begin
          rdx_convertlog, downloadpath, localpath, average=16 $
;                          , dx=31.92, dy=14.81, rotation=84.87, scale=4.935
                          , dx=dx, dy=dy, scale=scale $
                          , rotation = 84.599 
          ;; Switched to dx,dy,scale parameters averaged from this
          ;; year's limb calibrations by Pit on 2 Aug 2023. Also
          ;; switched to rotation parameter used in the PIG/turret
          ;; code.
        endif else begin
            print, 'Download failed, no PIG log file to convert!'
            status = 1
            return
        endelse
      endif
    end

    else :
    
  endcase

  status = ~downloadOK
  
end

if 1 then begin
  red_download_log, 'pig', '2021-04-30', 'tmp/', /over $
                    , localpath = localpath $                      
                    , status = status
endif else begin
  red_download_log, 'r0', '2016-09-11', 'tmp/', /over $
                    , localpath = localpath $                      
                    , status = status
endelse


if status eq 0 then begin
  help, localpath
endif


end

