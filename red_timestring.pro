; docformat = 'rst'

;+
; Make a string representation of time in seconds since midnight.
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl
; 
; 
; :Returns:
; 
;    A string on the form hh:mm:ss.[ddd...], where 
; 
; :Params:
; 
;     t : in, type=double
; 
; 
; :Keywords:
; 
;     Nsecdecimals : in, optional, type=integer, default=5
;   
;        Use this number of decimals for the seconds.
;
;     interval : in, optional, type=boolean
;
;        Set this keyword to interpret t as a time interval, producing
;        a different type of string. 
; 
; 
; :History:
; 
;     2016-09-22 : MGL. First version, based on red_time2double.
;                  Specify number of decimals by using keyword. Allow
;                  to wrap to day(s) before (negative numbers) or
;                  day(s) after (large numbers).
; 
;     2016-11-16 : MGL. New keyword short.
; 
;     2016-11-28 : MGL. New keyword interval.
; 
; 
;-
function red_timestring, t, Nsecdecimals = Nsecdecimals, short = short, interval = interval

  isscalar = size(t, /n_dim) eq 0
  n = n_elements(t)

  if n_elements(Nsecdecimals) eq 0 then Nsecdecimals = 5

  ;; Make a format string for hours and minutes (and maybe seconds)
  format = '(I02)'

  ;; Make a format string for the seconds
  if Nsecdecimals gt 0 then begin
     secspace = Nsecdecimals+3
     secformat = '(F0'+strtrim(secspace, 2)+'.'+strtrim(Nsecdecimals, 2)+')'
  endif else begin
     secformat = format
  endelse

  res = strarr(n)
  for i = 0, n-1 do begin
     secperday = 24d * 3600d
     tt = t[i] - floor(t[i]/secperday) * secperday ; Like mod secperday but works for negative numbers
     it = long(tt)
     hours = it / 3600L
     min = (it  - hours * 3600L) / 60L
     secs = tt - hours * 3600L - min * 60L

     if keyword_set(short) then begin
       if hours gt 0 then res[i] += string(hours, format = format) + ':'
       if hours gt 0 or min gt 0 then res[i] += string(min, format = format) + ':' 
       res[i] += string(secs, format = secformat)
     endif else if keyword_set(interval) then begin
       if hours gt 0 then res[i] += string(hours, format = format) + 'h'
       if hours gt 0 or min gt 0 then res[i] += string(min, format = format) + 'm' 
       if min gt 0 then begin
         res[i] += string(secs, format = format) + 's'
       endif else begin
         res[i] += string(secs, format = secformat) + 's'
       endelse
     endif else begin
       res[i] += red_stri(hours, ni = format) + ':'
       res[i] += red_stri(min, ni = format) + ':' 
       res[i] += string(secs, format = secformat)
     endelse


;     if hours gt 0 or ~keyword_set(short) then res[i] += red_stri(hours, ni = format) + ':'
;     if hours gt 0 or min gt 0 or ~keyword_set(short) then res[i] += red_stri(min, ni = format) + ':' 
;     res[i] += string(secs, format = secformat)
;     res[i] = red_stri(hours, ni = format) + ':' +red_stri(min, ni = format) + $
;              ':' + string(secs, format = secformat)
  endfor
  

  if isscalar then return, res[0] else return, res
  
end
