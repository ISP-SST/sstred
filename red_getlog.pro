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
; 
; 
; 
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;    2017-04-21 : MGL. First version.
; 
; 
; 
; 
;-
pro red_getlog, date $
                , r0 = r0 $
                , turret = turret $
                , pig = pig $
                , temp = temp $
                , weather = weather

  ;; Make sure we have the date in ISO format
  isodate = (strsplit(date_conv(red_strreplace(date, '.', '-', n = 2), 'F'), 'T', /extract))[0]
  
  ;; Times in the log files are given in Unix system time = # of
  ;; seconds since 1 Jan 1970 UTC. Need the following to convert them
  ;; to seconds after midnight.
  yr = long((strsplit(isodate, '-', /extract))[0])
  mo = long((strsplit(isodate, '-', /extract))[1])
  dy = long((strsplit(isodate, '-', /extract))[2])
  CDF_EPOCH, epoch, 1970, 1, 1, 0, /COMPUTE_EPOCH    ; In milliseconds
  CDF_EPOCH, midnight, yr, mo, dy, 0, /COMPUTE_EPOCH ; In milliseconds
  midnight = (midnight-epoch)/1000.d                 ; Unix system time for midnight today.

  ;; Temperature sensors in observing room
  if arg_present(temp) then begin
    print, 'red_getlog : Temperature not implemented yet'
  endif
  
  ;; Weather
  if arg_present(weather) then begin
    print, 'red_getlog : Weather not implemented yet'
  endif
  
  
  ;; Fried parameter r0
  
  if arg_present(r0) then begin

    print, 'red_getlog : Get r0 data.'
    
    red_download, date = isodate, /r0, pathr0 = r0file
    if r0file then begin

      ;; Read r0 file
      if isodate ge '2013-10-28' then begin
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

      pigdata = read_ascii(pigfile, TEMPLATE=pigtemplate)

      ;;  Remove glitches
      indx = where(pigdata.time ne 0, Npig) 
      
      if Npig ne 0 then begin

        ;; Change struct of arrays to array of structs
        pig = replicate(pigstruct, Npig)

        pig.time = pigdata.time[indx] - midnight
        pig.x = pigdata.x[indx]
        pig.y = pigdata.y[indx]
        
      endif                     ; Npig
      
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
                                       , 'stony_ns_direction', 'stony_ns' $
                                       , 'stony_ew_direction', 'stony_ew' $
                                       , 'parala_angle', 'solar_local_tilt' $
                                       , 'az_th', 'el_th'] $
                         , FIELDLOCATIONS:[0L, 11L, 22L, 32L, 41L, 44L, 50L, 53L, 59L, 69L, 80L, 89L] $
                         , FIELDGROUPS:lindgen(12) $
                       }

      turretstruct = {time:0d $
                      , az:0., el:0. $
                      , stony_n:0., stony_w:0. $
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

        turret.stony_n = turretdata.stony_ns
        negindx = where(turretdata.stony_ns_direction eq 'S', count)
        if count gt 0 then turret[indx].stony_n *= -1
        
        turret.stony_w = turretdata.stony_ew
        negindx = where(turretdata.stony_ew_direction eq 'E', count)
        if count gt 0 then turret[indx].stony_w *= -1

        
      endif                     ; Nturret
      
    endif                       ; turretfile
  endif                         ; Turret

end
