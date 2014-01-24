; docformat = 'rst'

;+
; Get a free portnumber
;
; :Author:
;
;   Mats Löfdahl, Institute for Solar Physics, mats@astro.su.se.
;
; :Keywords:
;
;   port : in, type=integer
;       the first free port >= port will be returned
;
; :History:
;
;   2014-01-24 : TH: Made into a free function
;
;-
function red_freeport, port=port
  
    if keyword_set(port) then begin 
        freeport = long(port)
    end else begin 
        freeport = long(32765)          ; some arbitrary default value > 1024
    end
    
    spawn,'netstat -atn',netstat
    Nn = dimen(netstat, 0)
    used_ports = lonarr(Nn)
    for i = 2, Nn-1 do used_ports[i] = long(last(strsplit(token(netstat[i],4),':',/extract)))
    used_ports = used_ports[uniq(used_ports,sort(used_ports))] ; Uniquify
    
    while inset(freeport, used_ports) do freeport += 1 

    return, freeport

end
