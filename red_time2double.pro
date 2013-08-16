; docformat = 'rst'

;+
; Convert a time string to the number of seconds after midnight (or
; the reverse operation).
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
;    The number of seconds after midnight or (if the dir keyword is
;    set) the time as a string.
; 
; :Params:
; 
;    t : in
;   
;       Either a time as a string or (if the dir keyword is set) the
;       number of seconds after midnight.
;   
; 
; :Keywords:
; 
;    dir  : in, optional, type=boolean
;   
;       Set this to do the reverse operation, make a string.
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-08-16 : MGL. Added documentation.
; 
; 
;-
function red_time2double, t, dir = dir
;
  if(~keyword_set(dir)) then begin
     tmp = strsplit(t, ':.', /extract)
     res = double(tmp[0]) * 3600d0 + double(tmp[1]) * 60d0 + $
           double(tmp[2]) + double('0.'+tmp[3])
  endif else begin
     it = long(t)
     frac = t - it
     hours = it / 3600L
     min = (it  - hours * 3600L) / 60L
     secs = t - hours * 3600L - min * 60L
     ni = '(I02)'
     res = red_stri(hours, ni = ni) + ':' +red_stri(min, ni = ni) + $
           ':' + string(secs+frac, format = '(F06.3)')
  endelse
;
  return, res
end
