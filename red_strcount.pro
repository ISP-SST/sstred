; Limited form of ANA's strcount: Counts the number of
; occurences of a character in a string.

; Could make it more general using strskp.

; 2014-01-22 : MGL. Renamed for inclusion in the red_ namespace.
;
; 2014-03-05 : THI. Handle both arrays of strings and single strings.
;
; 2014-04-02 : MGL. Return scalar if input was a single string. 
;
function red_strcount,st,pat1

  if strlen(pat1) eq 0 then return,0

  ;; pat1 is not a single character

  cnt = intarr(size([st],/dim))
  st0 = st
  st1 = red_strskp(st0,pat1)

  while not array_equal(st0,st1) do begin
     
     cnt += (st0 ne st1)

     st0 = st1
     st1 = red_strskp(st0,pat1)
 
     ;print,cnt,' "',st0,'"   "',st1,'"'

  end

  if size(st, /n_dim) eq 0 then cnt = cnt[0]

  return,cnt


end
