; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    time : 
;   
;   
;   
;    date : 
;   
;   
;   
; 
; :Keywords:
; 
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2014-05-07 : MGL. Allow scalar date. 
;
;-
function red_lp_angles, time, date

  ave = dblarr(n_elements(time))

  for i = 0L, n_elements(time) - 1 do begin
                                ;head=fzhead(names[i])
                                ;head=strsplit(head,' ,=,-,: ',/extract)
                                ;get_time,head,yy,mm,dd,hh,mi,ss

     tt = double(strsplit(time[i], ':', /extract))
     
     if n_elements(date) eq 1 then begin
        dat = double(strsplit(date, '-.', /extract))
     endif else begin
        dat = double(strsplit(date[i], '-.', /extract))
     endelse
                                ;
     red_get_sun, dat[0], dat[1], dat[2] $
                  , red_reform_frac_time(tt[0], tt[1], tt[2]) $
                  , ha, dec

     red_get_azel, ha, dec, az, el
     ave[i] = red_get_rot(az, el, dec)

  endfor

  return, ave

end
