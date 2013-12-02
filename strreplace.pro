; ANA utility function

; Substitute pat2 for pat1 in st. If n is set, do it n times,
; otherwise do it once.  

FUNCTION strreplace,st,pat1,pat2,n=n

  if n_elements(n) eq 0 then n=1

  m = strcount(st,pat1)

  if m eq 0 then return,st

  ;;print,m,n

  st0 = st

  for i=0,min([m,n])-1 do begin
     pos = strpos(st0,pat1)
     st0 = strmid(st0,0,pos)+pat2+strskp(st0,pat1)
     ;;print,i,' ',st0
  end

  return,st0

END
