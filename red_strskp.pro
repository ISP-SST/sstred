; ANA's strskp function

; 2014-01-22 : MGL. Renamed for inclusion in the red_ namespace.

function red_strskp,st,del

  pos = strpos(st,del)

  if pos ge 0 then begin
     return,strmid(st,pos+strlen(del))
  end else begin
     return,st
  end

end
