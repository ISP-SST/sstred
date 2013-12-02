; ANA's strskp function

function strskp,st,del

  pos = strpos(st,del)

  if pos ge 0 then begin
     return,strmid(st,pos+strlen(del))
  end else begin
     return,st
  end

end
