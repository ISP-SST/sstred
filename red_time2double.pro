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
;    t : in, type="string, float, or double or array of any of those"
;   
;       Scalar or array representing time, either as a string or (if
;       the dir keyword is set) the number of seconds after midnight.
;       In the string case, the following format is assumed:
;       hh[:mm[:ss[.ddd]]].
;   
; 
; :Keywords:
; 
;    dir : in, optional, type=boolean
;   
;       Set this to do the reverse operation, make a string.
;   
; 
;    inverse : in, optional, type=boolean
;   
;       Set this to do the reverse operation, make a string.
;   
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-08-16 : MGL. Added documentation.
; 
;   2013-11-11 : MGL. In the string making part, fraction of seconds
;                was added to number of seconds already including the
;                fraction. Fixed this.
; 
;   2014-04-30 : MGL. Input data t can now be an array. Also, be more
;                relaxed about t's format.
; 
;   2014-05-31 : MGL. New inverse keyword, does the same as dir
;                keyword. 
; 
;-
function red_time2double, t, dir = dir, inverse = inverse

  isscalar = size(t, /n_dim) eq 0
  n = n_elements(t)

  if keyword_set(dir) then inverse = 1
  
  if(~keyword_set(inverse)) then begin
     res = dblarr(n)
     for i = 0, n-1 do begin
        tmp = strsplit(t[i], ':.', /extract, count = Nfields)
        res[i] = double(tmp[0]) * 3600d0
        if Nfields ge 2 then res[i] += double(tmp[1]) * 60d0
        if Nfields ge 3 then res[i] += double(tmp[2])
        if Nfields ge 4 then res[i] += double('0.'+tmp[3])
     endfor
  endif else begin
     res = strarr(n)
     for i = 0, n-1 do begin
        it = long(t[i])
        hours = it / 3600L
        min = (it  - hours * 3600L) / 60L
        secs = t[i] - hours * 3600L - min * 60L
        ni = '(I02)'
        res[i] = red_stri(hours, ni = ni) + ':' +red_stri(min, ni = ni) + $
                 ':' + string(secs, format = '(F06.3)')
     endfor
  endelse
  
  if isscalar then return, res[0] else return, res
  
end
