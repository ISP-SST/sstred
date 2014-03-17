; ANA's strskp function

; 2014-01-22 : MGL. Renamed for inclusion in the red_ namespace.
;
; 2014-03-05 : THI. Handle both arrays of strings and single strings.

function red_strskp,st,del

    pos = [strpos(st,del)]

    if size(pos,/dim) gt 0 then begin
        return,strmid(st,transpose(pos+(pos ge 0)*strlen(del)))
    end else begin
       return,st
    end

end
