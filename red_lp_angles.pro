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
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2014-05-07 : MGL. Allow scalar date. Add to the documentation. 
;
;-
function red_lp_angles, time, date

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

end
