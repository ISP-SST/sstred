; docformat = 'rst'

;+
; Returns various quantities related to SST observations for specific
; moments in time. This includes data that various instruments write
; in their log files (AO, PIG). Also calculates solar distance and
; radius in arcsec.
; 
; :Categories:
;
;    SST observations
; 
; 
; :Author:
; 
;    Mats Löfdahl, 2013-11-11
; 
; 
; :Params:
; 
;    date : in, out, optional, type=string
; 
;      Should be a string representing the day for which the logged
;      data is wanted. If an image file header is given with the
;      header keyword, date will be calculated from this and returned
;      with a proper date and time string.
; 
;    time : in, out, optional, type="string or double or dblarr(N)"
; 
;      Should be a string (hh:mm:ss) or float (seconds since midnight)
;      representing the time of day for which the logged data is
;      wanted. If an image file header is given with the header
;      keyword, time will be calculated from this and returned. Arrays
;      are also permitted, then arrays with the return values will be
;      generated. If not given, time coordinates of r0 log file will
;      be used if r0 keyword is given. Else time coordinates of PIG
;      log file.
; 
; :Keywords:
;
;   azel : out, optional, type=array
;
;      Azimuth and elevation.
; 
;   header : in, optional, type=string
;
;     An SST image header (or an array of N headers). If time is not
;     given, it will be extracted from this header. Same for date.
;
;   r0 : out, optional, type="fltarr(2) or fltarr(2,N)"
; 
;     The r0 values in meters from the AO (24x24 and 8x8 pixels) will
;     be returned for 2013-10-28 and later. Before that day no 8x8
;     data exist and the dimensions of the returned array is changed
;     accordingly. 
;
;   turret : out, optional, type="fltarr(2) or fltarr(2,N)"
;
;     Turret pointing info in xxx coordinates.
;
;   pig : out, optional, type="fltarr(2) or fltarr(2,N)"
;
;     PIG pointing info in xxx coordinates.
;
;   diskpos : out, optional, type="fltarr(2) or fltarr(2,N)"
;
;     SST pointing info in xxx coordinates. Based on PIG log data but
;     will use turret data if PIG data is not available or is deemed
;     unreliable. The latter can happen when the PIG software has
;     problems fitting a circle to the solar disk, like near the
;     horizon and possibly with clouds. If no pig_T or pig_N keywords
;     are given, pig_T=1 (second) will be assumed.
; 
; 
;   mu : out, optional, type="float or fltarr(N)"
; 
;       mu=cos(theta) based on PIG coordinates.
; 
;   Rsun : out, optional, type=double
;
;       Solar radius in arcsec.
;
;
;   Sundist : out, optional, type=double
;
;       Solar distance in meters.
;
;   Stonyhurst : out, optional, type=array
;
;       Stonyhurst coordinates.
;
;   zenithangle : out, optional, type=float
;
;       Pointing angular distance from zenith (= 1-elevation) in
;       degrees. 
;
;
; :History:
; 
;    2013-11-14 : MGL. Two new keywords, pig_T and pig_N. Work
;                 on the PIG part of the code. Made it work with
;                 arrays and not just a single point in time.
;
;    2013-11-15 : MGL. Added code from R. Sterner's and S. Freeland's
;                 get_sun.pro to calculate solar distance in meters
;                 and solar radius in arcsec for a given date.
;
;    2013-12-10 : MGL. Renamed from sst_logdata for inclusion in the
;                 crispred pipeline. Added turret keyword and as
;                 code to download and read the turret log file. Added
;                 pointing keyword and code to evaluate the PIG log
;                 values and use turret log values as fallback.
;
;    2013-12-20 : MGL. Remove pig_N and pig_T keywords. Let
;                 red_download do the downloading and the PIG logfile
;                 post processing. Get date from PWD. Add code for the
;                 Turret log. Start implementing the "pointing"
;                 keyword.
;
;    2014-01-22 : MGL. Adapt to string functions moved to the str_
;                 namespace.
;
;    2014-04-30 : MGL. Implemented the turret keyword, added and
;                 implemented the zenithangle keyword.
;
;    2014-05-28 : MGL. Give the date when calling red_download.
;
;    2014-10-10 : MGL. May have to look for r0 log files in
;                 YYYY-subdirectories on the web site. 8x8 r0 data
;                 does not exist before 2013-10-28.
;
;     2014-10-10 : MGL. Use read_ascii rather than mgl_rd_tfile for
;                  reading the r0, PIG, and turret log files.
;
;     2015-08-18 : MGL. Corrected a bug that popped up when you want
;                  interpolated r0 values at specified times.
;
;     2017-04-22 : MGL. Use red_logdata to read the log files.
;
;     2017-04-23 : MGL. New keywords use_r0_time, use_turret_time,
;                  use_pig_time, azel, stonyhurst, diskpos. Removed
;                  keyword pointing.
;
;
;-
pro red_logdata, date, time $
                 , azel = azel $
                 , diskpos = diskpos $
                 , header = header $
                 , mu = mu $
                 , pig = pig $
                 , r0 = r0 $
                 , Rsun = Rsun $
                 , stonyhurst = stonyhurst $
                 , Sundist = Sundist $
                 , turret = turret $
                 , zenithangle = zenithangle $
                 , use_pig_time = use_pig_time $
                 , use_turret_time = use_turret_time $
                 , use_r0_time = use_r0_time

  ;; Date in ISO format
  if n_elements(date) gt 0 then begin
    isodate = (strsplit(date_conv(red_strreplace(date, '.', '-', n = 2), 'F') $
                        , 'T', /extract))[0]
  endif else begin
    if n_elements(header) ne 0 then begin
      Ts = strmid(header[0], strpos(header[0], 'Ts=')+3, 26)
      isodate = red_strreplace((strsplit(Ts, ' ', /extract))[0], '.', '-', n = 2)
      date = isodate
    endif else begin
      print, 'sst_logdata : No date given.'
      date = stregex(getenv('PWD'),'[12][0-9][0-9][0-9][-.][0-1][0-9][-.][0-3][0-9]' $
                     , /extract)
      if date eq '' then begin
        print, 'red_logdata : No date given and PWD does not contain a date.'
        retall
      endif
      isodate = red_strreplace(date, '.', '-', n = 2)
    endelse
  endelse

  ;; Times in the log files are given in Unix system time = # of
  ;; seconds since 1 Jan 1970 UTC. Need the following to convert them
  ;; to seconds after midnight.
  yr = long((strsplit(isodate, '-', /extract))[0])
  mo = long((strsplit(isodate, '-', /extract))[1])
  dy = long((strsplit(isodate, '-', /extract))[2])
;  CDF_EPOCH, epoch, 1970, 1, 1, 0, /COMPUTE_EPOCH    ; In milliseconds
;  CDF_EPOCH, midnight, yr, mo, dy, 0, /COMPUTE_EPOCH ; In milliseconds
;  midnight = (midnight-epoch)/1000.d                 ; Unix system time for midnight today.

  
  ;; Get the times (in seconds after midnight) for which the logdata are wanted.
  if n_elements(header) ne 0 then begin
    ;; From header(s) if given.
    Ntimes = n_elements(header)
    T = dblarr(Ntimes)
    for i = 0, Ntimes-1 do begin
      Ts = strmid(header[i], strpos(header[i], 'Ts=')+3, 26)
      Te = strmid(header[i], strpos(header[i], 'Te=')+3, 26)
      ;; Time in seconds after midnight
      Ts = (strsplit(Ts, ' ', /extract))[1] ; Just time part
      Te = (strsplit(Te, ' ', /extract))[1]
      T[i] = (red_time2double(Ts) + red_time2double(Te)) / 2.d
    endfor
    ;; Return value
    time = T
  endif else if n_elements(time) gt 0 then begin
    ;; From time keyword.
    Ntimes = n_elements(time)
    if size(time, /type) eq 7 then begin
      ;; String
      T = dblarr(Ntimes)
      for i = 0, Ntimes-1 do begin
        T[i] = red_time2double(time[i])
      endfor
    endif else begin
      ;; Numerical
      T = double(time)
    endelse
  endif else begin
    ;; No explicit time info given.
    Ntimes = 0
    print, 'sst_logdata : No time info given.'
    print, '              Will return data for time coordinates in r0 log file'
    print, '              if keyword_set(r0), else time coordinates in PIG file.'
  endelse
  
  if Ntimes eq 0 then begin
    print, 'Get SST log info for ' + isodate + '.'
  endif else if Ntimes eq 1 then begin
    print, 'Get SST log info for ' + isodate + ' at ' $
           + red_time2double(T, /dir) $
           + ' (' + strtrim(t, 2) + ' seconds after midnight).'
  endif else begin
    print, 'Get SST log info for ' + isodate + ' at ' $
           + strtrim(Ntimes, 2) + ' different points in time.'
  endelse

  
  ;; Calculation of solar radius and distance.
  
  AU = 149597870700.d0          ; [m] 
  ;;Rsun_m = 6.955e8                        ; [m]
  ;;Rsun = atan(Rsun_m/AU)*180/!pi*3600. ; [arcsec]
  ;; = 958.931"

  ;; BEGIN get_sun  ---------------------------------------------------;
  ;; Here follows calculation of solar radius, code copied from        ;
  ;; get_sun.pro by  R. Sterner and S. Freeland.                       ;
  ;; http://hesperia.gsfc.nasa.gov/ssw/gen/idl/solar/get_sun.pro       ;
  ;;                                                                   ;
  ;; Copyright (C) 1991, Johns Hopkins University/Applied Physics      ;
  ;; Laboratory This software may be used, copied, or redistributed as ;
  ;; long as it is not sold and this copyright notice is reproduced on ;
  ;; each copy made.                                                   ;
  ;;                                                                   ;
  jd = JULDAY(Mo, Dy, Yr)                                              ;
                                                                       ;
  ;; Julian Centuries from 1900.0:                                     ;
  jc = (jd - 2415020d)/36525d                                          ;
                                                                       ;
  ;; Carrington Rotation Number:                                       ;
  carr = (1.d0/27.2753d0)*(jd-2398167.d0) + 1.d0                       ;
                                                                       ;
  ;; Geometric Mean Longitude (deg):                                   ;
  mnl = 279.69668d0 + 36000.76892d0*jc + 0.0003025*jc^2                ;
  mnl = mnl mod 360d0                                                  ;
                                                                       ;
  ;; Mean anomaly (deg):                                               ;
  mna = 358.47583d0 + 35999.04975d0*jc $                               ;
        - 0.000150d0*jc^2 - 0.0000033d0*jc^3                           ;
  mna = mna mod 360d0                                                  ;
                                                                       ;
  ;; Eccentricity of orbit:                                            ;
  e = 0.01675104d0 - 0.0000418d0*jc - 0.000000126d0*jc^2               ;
                                                                       ;
  ;; Sun's equation of center (deg):                                   ;
  c = (1.919460d0 - 0.004789d0*jc - 0.000014d0*jc^2)*sin(mna/!radeg) $ ;
      + (0.020094d0 - 0.000100d0*jc)*sin(2*mna/!radeg) $               ;
      + 0.000293d0*sin(3*mna/!radeg)                                   ;
                                                                       ;
  ;; Sun's true geometric longitude (deg)                              ;
  true_long = (mnl + c) mod 360d0                                      ;
                                                                       ;
  ;; Sun's true anomaly (deg):                                         ;
  ta = (mna + c) mod 360d0                                             ;
                                                                       ;
  ;; Sun's radius vector (AU).                                         ;
  dist = 1.0000002d0*(1.d0 - e^2)/(1.d0 + e*cos(ta/!radeg))            ;
                                                                       ;
  ;; Semidiameter (arc sec):                                           ;
  Rsun = 959.63/dist                                                   ;
                                                                       ;
  ;; Distance to the sun in meters                                     ;
  Sundist = AU*dist                                                    ;
  ;; END get_sun ------------------------------------------------------;

  
  ;; Decide what time coordinates to use
  if ~keyword_set(use_r0_time) $
     && ~keyword_set(use_turret_time) $
     && ~keyword_set(use_pig_time) then begin
    ;; If we didn't ask for a particular time, then decide here:
    case 1 of
      n_elements(T) gt 0 : use_input_time = 1           ; Use input
      arg_present(r0) ne 0 : use_r0_time = 1         ; Every s
      arg_present(pig) ne 0 : use_pig_time = 1       ; Every s?
      arg_present(turret) ne 0 : use_turret_time = 1 ; Every 30 s
      else : begin
        print, 'red_logdata : No time coordinates!'
        return
      end
    endcase
  endif
  
  ;; Download log file data

  get_r0_file = keyword_set(use_r0_time) $
                || arg_present(r0)

  get_pig_file = keyword_set(use_pig_time) $
                 || arg_present(pig) $
                 || arg_present(mu) $
                 || arg_present(diskpos) $
                 || arg_present(pointing)

  get_turret_file = keyword_set(use_turret_time) $
                    || arg_present(turret) $
                    || arg_present(mu) $
                    || arg_present(azel) $
                    || arg_present(stonyhurst) $
                    || arg_present(diskpos) $
                    || arg_present(pointing) $
                    || arg_present(zenithangle)
  
  if get_r0_file then begin
    red_getlog, isodate, r0 = r0data
    if n_elements(r0data) eq 0 then $
       print, 'red_logdata : No r0 log file.'
  endif

  if get_pig_file then begin
    red_getlog, isodate, pig = pigdata
    if n_elements(pigdata) eq 0 then $
       print, 'red_logdata : No pig log file.'
  endif
  
  if get_turret_file then begin
    red_getlog, isodate, turret = turretdata
    if n_elements(turretdata) eq 0 then $
       print, 'red_logdata : No turret log file.'
  endif

  
  ;; Use the selected time
  case 1 of
    keyword_set(use_input_time) :
    keyword_set(use_r0_time) : begin
      T = r0data.time
      time = T
    end
    keyword_set(use_pig_time) : begin
      T = pigdata.time
      time = T
    end
    keyword_set(use_turret_time) : begin
      T = turretdata.time
      time = T
    end
    else : begin
      print, 'red_logdata : This should not happen!'
      stop
    end
  endcase
  Ntimes = n_elements(T)
  
  if n_elements(r0data) ne 0 then begin
   
    if keyword_set(use_r0_time) then begin
      ;; Return all values
      if n_elements(r0data.r0_8x8) ne 0 then begin
        r0 = fltarr(2, Ntimes)
        r0[1, *] = median(r0data.r0_8x8, dim = 1)
      endif else begin
        r0 = fltarr(1, Ntimes)
      endelse
      r0[0, *] = r0data.r0_24x24
    endif else begin
      ;; Get interpolated values
      if n_elements(r0data.r0_8x8) ne 0 then begin
        r0 = fltarr(2, Ntimes)
        r0[1, *] = interpol(median(r0data.r0_8x8, dim = 1), r0data.time, T)
      endif else begin
        r0 = fltarr(1, Ntimes)
      endelse
      r0[0, *] = interpol(r0data.r0_24x24, r0data.time, T)
    endelse

  endif

  
  if n_elements(pigdata) gt 0 then begin
    pig = fltarr(2, Ntimes)
    if keyword_set(use_pig_time) then begin
      ;; Return all values
      pig[0, *] = pigdata.x
      pig[1, *] = pigdata.y
    endif else begin
      ;; Get interpolated values
      pig[0, *] = interpol(pigdata.x, pigdata.time, T)
      pig[1, *] = interpol(pigdata.y, pigdata.time, T)
    endelse 

  endif

  
  if n_elements(turret) ne 0 then begin
    
    if keyword_set(use_turret_time) then begin
      ;; Return all values
      turret = transpose([[turretdata.az],[turretdata.el]]) ; az/el on sky
    endif else begin
      ;; Get interpolated values
      turret = fltarr(2, Ntimes) ; az/el on sky
      turret[0, *] = interpol(turretdata.az, turret_time, T)
      turret[1, *] = interpol(turretdata.el, turret_time, T)
    endelse 
    
  endif

  
  ;; Fill keywords

  ;; Azimuth and elevation
  if (arg_present(azel) || arg_present(zenithangle)) $
     && n_elements(turret) then azel = turret
  
  ;; Disk coordinates (Helioprojective-Cartesian coordinates)
  if arg_present(diskpos) || arg_present(mu) then begin
    case 1 of
      n_elements(pig) gt 0 : begin
        ;; Use PIG coordinates if available.
        diskpos = pig
      end
      n_elements(turret) gt 0 : begin
        ;; Use Disk X/Y from the turret log file
        diskpos = red_turret_select_pointing(turretdata, 'Disk', time = T)
        ;; If we did not get data for all T points, maybe try to fill
        ;; in data by calculating from Stonyhurst?
      end
      else : stop
    endcase
  endif
      
  ;; Stonyhurst cordinates
  if arg_present(stonyhurst) then begin
    case 1 of
      ;; Prefer if available to calculate Stonuhurst coordinates from
      ;; pig disk coordinates?
      n_elements(turret.stony) : stonyhurst = red_turret_select_pointing(turretdata, 'Stony', time = T)
    endcase
  endif

  ;; Zenith angle
  if arg_present(zenithangle) && n_elements(azel) then $
     zenithangle = 90.0 - reform(azel[1, *])
  
  ;; µ
  if arg_present(mu) && n_elements(diskpos) then begin
    diskr = sqrt(total(diskpos^2,1))/Rsun ; Normalized radial coordinate
    mu = sqrt(1 - diskr^2)
;    theta = asin(diskr)                   ; 
;    mu = cos(theta)
  endif


  ;; Get the "best" pointing info from the PIG log file and the Turret
  ;; log file. The idea is to use PIG pointing info whenever possible
  ;; and use Turret pointing info when a) there is no PIG log file
  ;; (old data) or b) the PIG lost tracking for an extended period of
  ;; time. 
  if arg_present(pointing) then begin

    print, 'red_logdata : keyword "pointing" is not implemented yet.'
    stop
    
  endif

end
