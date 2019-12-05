; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
; 
; 
; 
; :Returns:
; 
;    The image rotation angle of the SST, in radians.
; 
; :Params:
; 
;    time : in, type=strarr
;   
;       The time of day as strings, hh:mm:ss.ddd.
;   
;    date :  in, type="string or strarr"
;   
;       The date(s) in iso format, yyyy-mm-dd. If it is an array, it
;       should have (at least) the same number of elements as time.
;   
; 
; :Keywords:
; 
;    from_log : in, optional, type=boolean
;
;       Base calculations on data from the turret log file.
;
;    offset_angle : in, optional, type=scalar
;
;       Angle offset in degrees to make the output equal to "the angle
;       between N/S of the sun and the horizontal, on the optical
;       table, pointing right looking at the sun." Aka "table
;       constant".(Used only together with /from_log.)
;
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2014-05-07 : MGL. Allow scalar date. Add to the documentation. 
; 
;   2019-12-05 : MGL. New keywords from_log and offset_angle.
;-
function red_lp_angles, time, date, from_log = from_log, offset_angle = TC

  if keyword_set(from_log) then begin

    ;; Get az, el, and tilt angle (in radians) from the turret log 
    red_logdata, date, time, tilt = TL, azel = azel
     
    if n_elements(TC) eq 0 then TC = 48.d0 ; table constant in degrees
  
    ;; From get_sun_angle.pro by van Noort, de la Cruz Rodriguez, and
    ;; Schnerr, 2010:
    ang = (azel[0, *] - azel[1, *] - TL - TC) * !dtor

    return, ang
    
  endif else begin

    ave = dblarr(n_elements(time))

    for i = 0L, n_elements(time) - 1 do begin

      tt = double(strsplit(time[i], ':', /extract))
      
      if n_elements(date) eq 1 then begin
        dat = double(strsplit(date, '-.', /extract))
      endif else begin
        dat = double(strsplit(date[i], '-.', /extract))
      endelse

      ;; Get the position of the Sun on the sky.
      red_get_sun, dat[0], dat[1], dat[2] $
                   , red_reform_frac_time(tt[0], tt[1], tt[2]) $
                   , ha, dec

      ;; Convert equatorial coordinates to azimut and elevation.
      red_get_azel, ha, dec, az, el

      ;; Compute rotation angle
      ave[i] = red_get_rot(az, el, dec)

    endfor

    return, ave

  endelse
  
end
