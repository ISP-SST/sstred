; docformat = 'rst'

;+
; Return SST log data.
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
; :Returns:
; 
; 
; :Params:
; 
;    date : in, type=string
; 
;      The date for which the log data should be read.
; 
; 
; :Keywords:
; 
;    r0 : out, optional, type=structarr
;   
;       Shack-Hartmann r0 data returned in this variable.
;   
;    turret : out, optional, type=structarr
;   
;       Turret positionlog data returned in this variable.   
;   
;    pig : out, optional, type=structarr 
;   
;       PIG data returned in this variable, undefined if unsuccessful
;       read.   
;   
;    shabar : out, optional, type=structarr
;   
;       SHABAR r0 data returned in this variable.

;    temp : out, optional, type=structarr
;   
;       Data from temperature sensors returned in this variable.   
;   
;    weather : out, optional, type=structarr
;   
;       Weather data returned in this variable.     
;   
;    wfwfs : out, optional, type=structarr
;   
;       Wide-field wavefront sensor r0 data returned in this variable.
; 
; 
; :History:
; 
;    2017-04-21 : MGL. First version.
; 
;    2017-04-22 : MGL. Take the different types of pointing info in
;                 the turret log into account. Implement reading log
;                 files for the SHABAR, for the weather station, and
;                 for the temperature sensors.
; 
;    2017-10-02 : MGL. Check pig data for successful read.
; 
;    2023-09-25 : MGL. New keyword wfwfs.
; 
;-
pro red_getlog, date $
                , pig = pig $
                , r0 = r0 $
                , shabar = shabar $
                , turret = turret $
                , temp = temp $
                , weather = weather $
                , wfwfs = wfwfs

  ;; Make sure we have the date in ISO format
  isodate = (strsplit(date_conv(red_strreplace(date, '.', '-', n = 2), 'F'), 'T', /extract))[0]
  
  ;; Times in several log files are given in Unix system time = # of
  ;; seconds since 1 Jan 1970 UTC. Need the following to convert them
  ;; to seconds after midnight.
  yr = long((strsplit(isodate, '-', /extract))[0])
  mo = long((strsplit(isodate, '-', /extract))[1])
  dy = long((strsplit(isodate, '-', /extract))[2])
  CDF_EPOCH, epoch, 1970, 1, 1, 0, /COMPUTE_EPOCH    ; In milliseconds
  CDF_EPOCH, midnight, yr, mo, dy, 0, /COMPUTE_EPOCH ; In milliseconds
  midnight = (midnight-epoch)/1000.d                 ; Unix system time for midnight today.
  
  ;; WFWFS - Wide-field wavefront sensor r0 data
  if arg_present(wfwfs) then begin

    print, 'red_getlog : Get WFWFS data'
    ;; Columns: timestamp "triangle" (5s average) "triangle" (30s
    ;; average) "triangle" (180s average) "triangle_sub" (5s average)
    ;; "triangle_sub" (30s average) "triangle_sub" (180s average)
    red_download, date = isodate, /wfwfs, pathwfwfs = wfwfsfile
    if wfwfsfile then begin
      wfwfstemplate = { VERSION:1.0 $
                        , DATASTART: 7 $
                        , DELIMITER: 32B $
                        , MISSINGVALUE: !VALUES.F_NAN $
                        , COMMENTSYMBOL: "#" $
                        , FIELDCOUNT: 7 $
                        , FIELDTYPES: [5, 4, 4, 4, 4, 4, 4] $
                        , FIELDNAMES: [ "time", "r0_5s", "r0_30s", "r0_180s" $
                                        , "r0_sub_5s", "r0_sub_30s", "r0_sub_180s"] $
                        , FIELDLOCATIONS: [0L, 18L, 27L, 36L, 45L, 54L, 63L] $
                        , FIELDGROUPS: [0L, 1L, 2L, 3L, 4L, 5L, 6L] $
                      }
      wfwfsstruct = {time:0d, r0_5s:0., r0_30s:0., r0_180s:0., r0_sub_5s:0., r0_sub_30s:0., r0_sub_180s:0.}
      wfwfsdata = read_ascii(wfwfsfile, TEMPLATE=wfwfstemplate)
      indx = where(wfwfsdata.time ne 0, Nwfwfs)
      if Nwfwfs ne 0 then begin
        wfwfs = replicate(wfwfsstruct, Nwfwfs)

        wfwfs.time        = wfwfsdata.time[indx] - midnight ; In seconds since midnight
        wfwfs.r0_5s       = wfwfsdata.r0_5s[indx]           ; Ground layer r0, 5s average  
        wfwfs.r0_30s      = wfwfsdata.r0_30s[indx]          ; Ground layer r0, 30s average 
        wfwfs.r0_180s     = wfwfsdata.r0_180s[indx]         ; Ground layer r0, 180s average
        wfwfs.r0_sub_5s   = wfwfsdata.r0_sub_5s[indx]       ; AO WFS FOV r0, 5s average  
        wfwfs.r0_sub_30s  = wfwfsdata.r0_sub_30s[indx]      ; AO WFS FOV r0, 30s average 
        wfwfs.r0_sub_180s = wfwfsdata.r0_sub_180s[indx]     ; AO WFS FOV r0, 180s average
      endif
    endif
  endif
  
  
  ;; SHABAR
  if arg_present(shabar) then begin
    
    print, 'red_getlog : Get SHABAR data'
    
    red_download, date = isodate, /shabar, pathshabar = shabarfile
    if shabarfile then begin
      
      shabartemplate = { VERSION:1.0 $
                         , DATASTART:0L $
                         , DELIMITER:32B $
                         , MISSINGVALUE:!VALUES.F_NAN $
                         , COMMENTSYMBOL:'' $
                         , FIELDCOUNT:24L $
                         , FIELDTYPES:[5L, replicate(4L, 23)] $
                         , FIELDNAMES:['time', 'what', 'why', 'trump' $
                                       , replicate('val', 19), 'intensity'] $
                         , FIELDLOCATIONS:[0L, 18L, 31L, 44L, 14L*lindgen(20)+58L] $
                         , FIELDGROUPS:[0L, 1L, 2L, 3L, replicate(4L, 19), 5L] $
                       }
      shabarstruct = {time:0d, what:0., why:0., trump:0., val:fltarr(19), intensity:0.}
      shabardata = read_ascii(shabarfile, TEMPLATE=shabartemplate)

      ;; Remove glitches
      indx = where(shabardata.time ne 0, Nshabar)
      
      if Nshabar ne 0 then begin

        ;; Change struct of arrays to array of structs

        shabar = replicate(shabarstruct, Nshabar)

        shabar.time      = shabardata.time[indx] - midnight ; In seconds since midnight
        shabar.what      = shabardata.what[indx]
        shabar.why       = shabardata.why[indx]
        shabar.trump     = shabardata.trump[indx]
        shabar.val       = shabardata.val[*, indx]
        shabar.intensity = shabardata.intensity[indx]
        
      endif                     ; Nshabar

    endif
  endif
  
  ;; Temperature sensors in observing room
  if arg_present(temp) then begin
    
    print, 'red_getlog : Get temperature sensor data.'

    red_download, date = isodate, /temp, pathtemp = tempfile
    if tempfile then begin

      temptemplate = { VERSION:1.0 $
                       , DATASTART:0L $
                       , DELIMITER:32B $
                       , MISSINGVALUE:!VALUES.F_NAN $
                       , COMMENTSYMBOL:'' $
                       , FIELDCOUNT:7L $
                       , FIELDTYPES:[replicate(7L, 6), 4L] $
                       , FIELDNAMES:['month', 'day', 'time', 'ignore', 'sensor', 'unit', 'temperature'] $
                       , FIELDLOCATIONS:[0L, 4L, 7L, 16L, 23L, 25L, 28L] $
                       , FIELDGROUPS:lindgen(7) $
                     }
      tempstruct = {time:0d, temperature:0., sensor:'', unit:''}
      tempdata = read_ascii(tempfile, TEMPLATE=temptemplate)

      ;; Change struct of arrays to array of structs
      Ntemp = n_elements(tempdata.time)

      if Ntemp gt 0 then begin

        temp = replicate(tempstruct, Ntemp)
        
        temp.time = red_time2double(tempdata.time) ; In seconds since midnight
        temp.temperature = tempdata.temperature
        temp.sensor = strtrim(tempdata.sensor, 2)
        temp.unit = red_strreplace(tempdata.unit, ':', '')

      endif                     ; Ntemp
      
    endif                       ; tempfile
  endif                         ; temperature
  
  ;; Weather
  if arg_present(weather) then begin

    print, 'red_getlog : Get weather data.'
      
    red_download, date = isodate, /weather, pathweather = weatherfile
    if weatherfile then begin

; For the 1-min average file:      
;      weathertemplate = { VERSION:1.0 $
;                          , DATASTART:0L $
;                          , DELIMITER:32B $
;                          , MISSINGVALUE:!VALUES.F_NAN $
;                          , COMMENTSYMBOL:'#' $
;                          , FIELDCOUNT:8L $
;                          , FIELDTYPES:replicate(5L, 8) $
;                          , FIELDNAMES:['time', 'winddir', 'windspeed', 'temperature' $
;                                        , 'humidity', 'pressure', 'rain', 'hail'] $
;                          , FIELDLOCATIONS:[0L,18L,29L,38L,47L,57L,68L,77L] $
;                          , FIELDGROUPS:lindgen(8) $
;                        }
      weathertemplate = { VERSION:1.0 $
                          , DATASTART:0L $
                          , DELIMITER:32B $
                          , MISSINGVALUE:!VALUES.F_NAN $
                          , COMMENTSYMBOL:'#' $
                          , FIELDCOUNT:8L $
                          , FIELDTYPES:[5, 3, replicate(4L, 6)] $
                          , FIELDNAMES:['time', 'winddir', 'windspeed', 'temperature' $
                                        , 'humidity', 'pressure', 'rain', 'hail'] $
                          , FIELDLOCATIONS:[0L,18L,22L,26,30L,35L,41L,45L] $
                          , FIELDGROUPS:lindgen(8) $
                        }
      weatherstruct = {time:0d, winddir:0d, windspeed:0d, temperature:0d $
                       , humidity:0d, pressure:0d, rain:0d, hail:0d}
      weatherdata = read_ascii(weatherfile, TEMPLATE=weathertemplate)

      ;; Remove glitches
      indx = where(weatherdata.time ne 0, Nweather)
      
      if Nweather ne 0 then begin

        ;; Change struct of arrays to array of structs

        weather = replicate(weatherstruct, Nweather)

        weather.time        = weatherdata.time[indx] - midnight ; In seconds since midnight
        weather.winddir     = weatherdata.winddir[indx]         ; [deg]
        weather.windspeed   = weatherdata.windspeed[indx]       ; [m/s]
        weather.temperature = weatherdata.temperature[indx]     ; [C]
        weather.humidity    = weatherdata.humidity[indx]        ; [%]
        weather.pressure    = weatherdata.pressure[indx]        ; [hPa]
        weather.rain        = weatherdata.rain[indx]            ; [mm/h]
        weather.hail        = weatherdata.hail[indx]            ; [mm/h]

      endif                     ; Nweather
      
    endif                       ; weatherfile
  endif                         ; weather
  
  
  ;; Fried parameter r0
  
  if arg_present(r0) then begin

    print, 'red_getlog : Get r0 data.'
    
    red_download, date = isodate, /r0, pathr0 = r0file
    if r0file then begin

      ;; Read r0 file
      if isodate ge red_dates(tag = 'AO 8x8') then begin
        ;; 8x8 measurements exist from this date.
        have8x8 = 1
        r0template = { VERSION:1.0 $
                       , DATASTART:0L $
                       , DELIMITER:32B $
                       , MISSINGVALUE:!VALUES.F_NAN $
                       , COMMENTSYMBOL:'' $
                       , FIELDCOUNT:13L $
                       , FIELDTYPES:[5L, replicate(4L, 12)] $
                       , FIELDNAMES:['time', 'r0_24x24', 'closedloop', 'subimages' $
                                     , replicate('r0_8x8', 9)] $
                       , FIELDLOCATIONS:[0L, 18L, 27L, 36L, 9L*lindgen(9)+46L] $
                       , FIELDGROUPS:[0L, 1L, 2L, 3L, replicate(4L, 9)] $
                     }
        r0struct = {time:0d, r0_8x8:fltarr(9), r0_24x24:0. $
                    , closedloop:0., subimages:0.}
        endif else begin
        have8x8 = 0
        r0template = { VERSION:1.0 $
                       , DATASTART:0L $
                       , DELIMITER:32B $
                       , MISSINGVALUE:!VALUES.F_NAN $
                       , COMMENTSYMBOL:'' $
                       , FIELDCOUNT:4L $
                       , FIELDTYPES:[5L, 4L, 4L, 4L] $
                       , FIELDNAMES:['time', 'r0_24x24', 'closedloop', 'subimages'] $
                       , FIELDLOCATIONS:[0L, 18L, 27L, 36L] $
                       , FIELDGROUPS:[0L, 1L, 2L, 3L] $
                     }
        r0struct = {time:0d, r0_24x24:0. $
                    , closedloop:0., subimages:0.}
      endelse

      r0data = read_ascii(r0file, TEMPLATE=r0template)

      ;; Remove glitches
      indx = where(r0data.time ne 0, Nr0)
      if Nr0 ne 0 then begin

        ;; Change struct of arrays to array of structs

        r0 = replicate(r0struct, Nr0)

        r0.time = r0data.time[indx] - midnight ; In seconds since midnight
        r0.r0_24x24 = r0data.r0_24x24[indx]
        r0.closedloop = r0data.closedloop[indx]
        r0.subimages = r0data.subimages[indx]

        if have8x8 then r0.r0_8x8 = r0data.r0_8x8[*, indx]
                        
      endif                     ; Nr0

    endif                       ; r0file
  endif                         ; r0

  
  ;; PIG log
  if arg_present(pig) then begin

    red_download, date = isodate, /pig, pathpig = pigfile

    if pigfile then begin

      pigtemplate = { VERSION:1.0 $
                      , DATASTART:0L $
                      , DELIMITER:32B $
                      , MISSINGVALUE:!VALUES.F_NAN $
                      , COMMENTSYMBOL:'#' $
                      , FIELDCOUNT:3L $
                      , FIELDTYPES:[5L, 4L, 4L] $
                      , FIELDNAMES:['time', 'x', 'y'] $
                      , FIELDLOCATIONS:[0L, 18L, 27L] $
                      , FIELDGROUPS:[0L, 1L, 2L] $
                    }

      pigstruct = {time:0d, x:0., y:0.} 

      pigdata = read_ascii(pigfile, TEMPLATE=pigtemplate, count = count)

      if count gt 0 then begin
        
        ;;  Remove glitches
        indx = where(pigdata.time ne 0, Npig) 
        
        if Npig ne 0 then begin

          ;; Change struct of arrays to array of structs
          pig = replicate(pigstruct, Npig)

          pig.time = pigdata.time[indx] - midnight ; [s]
          pig.x = pigdata.x[indx]                  ; [arcsec]
          pig.y = pigdata.y[indx]                  ; [arcsec]
                                
        endif                   ; Npig
      endif                     ; Successful read      
    endif                       ; pigfile
  endif                         ; PIG
 
  
  ;; Turret log
  if arg_present(turret) then begin
    
    red_download, date = isodate, /turret, pathturret = turretfile

    if turretfile then begin

      turrettemplate = { VERSION:1.0 $
                         , DATASTART:0L $
                         , DELIMITER:32B $
                         , MISSINGVALUE:!VALUES.F_NAN $
                         , COMMENTSYMBOL:'' $
                         , FIELDCOUNT:12L $
                         , FIELDTYPES:[7L, 7L, 4L, 4L, 7L, 4L, 7L, 4L, 4L, 4L, 4L, 4L] $
                         , FIELDNAMES:['date', 'time', 'az', 'el' $
;                                       , 'stony_ns_direction', 'stony_ns' $
;                                       , 'stony_ew_direction', 'stony_ew' $
                                       , 'pointtag1', 'pointing1' $
                                       , 'pointtag2', 'pointing2' $
                                       , 'parala_angle', 'solar_local_tilt' $
                                       , 'az_th', 'el_th'] $
                         , FIELDLOCATIONS:[0L, 11L, 22L, 32L, 41L, 44L, 50L, 53L, 59L, 69L, 80L, 89L] $
                         , FIELDGROUPS:lindgen(12) $
                       }

      turretstruct = {time:0d $
                      , az:0., el:0. $
                      , pointingtype1:'', pointingtype2:'' $
                      , pointing1:0., pointing2:0. $
                      , parala_angle:0. $
                      , solar_local_tilt:0. $
                      , az_th:0., el_th:0.}
      
      turretdata = read_ascii(turretfile, template = turrettemplate)

      ;; Turretdata is a structure with the following arrays: [date
      ;; YYYY/MM/DD, time HH:MM:SS, az, el, N or S, Stonyhurst N/S,
      ;; E or W, Stonyhurst E/W, parala. angle, solar local tilt,
      ;; Az(th), El(th)]. Angles are in degrees.
      
      ;; Az(th), El(th) are the coordinates where you want to point? 
      
      ;; In the Stonyhurst system the zero point is set at the
      ;; intersection of the Sun's equator and central meridian as
      ;; seen from the Earth. Longitude increases towards the Sun's
      ;; western limb. A solar feature will have a fixed latitude as
      ;; it rotates across the solar disk, but its longitude will
      ;; increase. This is in contrast to the Carrington
      ;; heliographic coordinate system, where the longitude remains
      ;; approximately fixed in time.
      ;; (http://en.wikipedia.org/wiki/Stonyhurst_Observatory#Stonyhurst_heliographic_coordinates) 
      
      Nturret = n_elements(turretdata.time)
      if Nturret ne 0 then begin

        ;; Change struct of arrays to array of structs
        turret = replicate(turretstruct, Nturret)
        
        turret.time = red_time2double(turretdata.time) ; Seconds since midnight
        
        turret.az = turretdata.az
        turret.el = turretdata.el
        turret.parala_angle = turretdata.parala_angle
        turret.solar_local_tilt = turretdata.solar_local_tilt
        turret.az_th = turretdata.az_th
        turret.el_th = turretdata.el_th
        turret.pointing1 = turretdata.pointing1
        turret.pointing2 = turretdata.pointing2

        ;; The pointing info comes in varying forms, as defined in the
        ;; pointtag<i> fields. We have for {pointtag1, pointing1,
        ;; pointtag2, pointing2}:
        
        ;; {'[NS]', value, '[EW]', value} : Stonyhurst
        ;; {'A',    value, 'E',    value} : Az/El tracking [deg]       
        ;; {'a',    value, 'e',    value} : Az/El wanted   [deg]       
        ;; {'X',    value, 'Y',    value} : Disk position tracking ["] 
        ;; {'x',    value, 'y',    value} : Disk position wanted   ["]  
        ;; {'f',    value, 'f',    value} : Flat field mode?       ["]

        ;; Stonyhurst
        indx = where( (turretdata.pointtag1 eq 'N' or turretdata.pointtag1 eq 'S') $
                      and (turretdata.pointtag2 eq 'E' or turretdata.pointtag2 eq 'W'), count )
        if count gt 0 then begin
          turret[indx].pointingtype1 = 'Stonyhurst N [deg]'
          turret[indx].pointingtype2 = 'Stonyhurst W [deg]'

          ;; Negative sign can be indicated by the opposite wind from
          ;; N or W:
          negindx = where(turretdata.pointtag1[indx] eq 'S', count)
          if count gt 0 then turret[indx[negindx]].pointing1 *= -1
          negindx = where(turretdata.pointtag2[indx] eq 'E', count)
          if count gt 0 then turret[indx[negindx]].pointing2 *= -1
        endif
        
        ;; Az/El tracking [deg]       
        indx = where( (turretdata.pointtag1 eq 'A') and (turretdata.pointtag2 eq 'E'), count )
        if count gt 0 then begin
          turret[indx].pointingtype1 = 'Azimuth [deg]'
          turret[indx].pointingtype2 = 'Elevation [deg]'
        endif
        
        ;; Az/El wanted [deg]           
        indx = where( (turretdata.pointtag1 eq 'a') and (turretdata.pointtag2 eq 'e'), count )
        if count gt 0 then begin
          turret[indx].pointingtype1 = 'Wanted Azimuth [deg]'
          turret[indx].pointingtype2 = 'Wanted Elevation [deg]'
        endif
        
        ;; Cartesian disk position tracking ["]
        indx = where( (turretdata.pointtag1 eq 'X') and (turretdata.pointtag2 eq 'Y'), count )
        if count gt 0 then begin
          turret[indx].pointingtype1 = 'Disk X ["]'
          turret[indx].pointingtype2 = 'Disk Y ["]'
        endif
        
        ;; Cartesian disk position wanted ["]
        indx = where( (turretdata.pointtag1 eq 'x') and (turretdata.pointtag2 eq 'y'), count )
        if count gt 0 then begin
          turret[indx].pointingtype1 = 'Wanted Disk X ["]'
          turret[indx].pointingtype2 = 'Wanted Disk Y ["]'
        endif
        
        ;; Flat field mode? ["]
        indx = where( (turretdata.pointtag1 eq 'f' or turretdata.pointtag1 eq 'F') $
                      and (turretdata.pointtag2 eq 'f' or turretdata.pointtag2 eq 'F'), count )
        if count gt 0 then begin
          turret[indx].pointingtype1 = 'Flat 1 ["]'
          turret[indx].pointingtype2 = 'Flat 2 ["]'
        endif
        
      endif                     ; Nturret
      
    endif                       ; turretfile
  endif                         ; Turret

end

cd, '/scratch/mats/2016.09.19/CRISP-aftersummer/'

date = '2016-09-19'
red_getlog, date, shabar = shabar
cgwindow
red_timeplot,/add, shabar.time,shabar.intensity, xtitle = 'time [UT]', ytitle = 'SHABAR intensity', title = date
;cgcontrol, out = 'shabar_intensity_'+date+'.pdf'


end

