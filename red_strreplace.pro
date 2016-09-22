; ANA utility function

; Substitute pat2 for pat1 in st. If n is set, do it n times,
; otherwise do it once.  

; 2014-01-22 : MGL. Renamed for inclusion in the red_ namespace.
;
; 2014-03-05 : THI. Handle both arrays of strings and single strings.
;
; 2016-09-20 : JLF. Fix bug in handling arrays of strings where some of the
;		strings don't contain the replaced character.

FUNCTION red_strreplace,st,pat1,pat2,n=n

  if n_elements(n) eq 0 then n=1

  m = red_strcount(st,pat1)

  if max(m) le 0 then return,st

  idx = where(m gt n, count)
  if count gt 0 then m[idx] = n

  st0 = st
  
  while max(m) gt 0 do begin
      idx = where(m gt 0)
      pos = strpos(st0[idx],pat1)
      st0[idx] = strmid(st0[idx],0,transpose(pos))+pat2+red_strskp(st0[idx],pat1)
      m -= 1
  endwhile

  return,st0

END
