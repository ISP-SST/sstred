; Limited form of ANA's strcount: Counts the number of
; occurences of a character in a string.

; Could make it more general using strskp.

function strcount,st,pat1

  if strlen(pat1) eq 0 then return,0
  if strlen(pat1) eq 1 then return,long(total(byte(st) eq (byte(pat1))[0]))

  ;; pat1 is not a single character

  cnt=0
  st0 = st
  st1 = strskp(st0,pat1)

  while st0 ne st1 do begin
     
     cnt += 1

     if strlen(st1) eq 0 then return,cnt

     st0 = st1
     st1 = strskp(st0,pat1)
 
     ;print,cnt,' "',st0,'"   "',st1,'"'

  end

  return,cnt


end
