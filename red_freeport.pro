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
;   2014-04-26 : MGL. Remove dependence on token() and inset()
;                functions. 
;
;   2014-05-15 : MGL. Don't use non-standard function "dimen". 
;   
;-
function red_freeport, port=port
  
    if keyword_set(port) then begin 
        freeport = long(port)
    end else begin 
       freeport = 32765L        ; some arbitrary default value > 1024
    end

    ;; Find used ports
    spawn,'netstat -atn',netstat
    Nn = n_elements(netstat)
    used_ports = lonarr(Nn)
    for i = 2, Nn-1 do $
       used_ports[i] = long((strsplit((strsplit(netstat[i],/extract))[3] $
                                      , ':', /extract $
                                      , count = Nsplit))[Nsplit-1])

    ;; Uniquify
    used_ports = used_ports[uniq(used_ports,sort(used_ports))] 

    ;; Increment freeport until it is no longer in used_ports
    while min(abs(freeport-used_ports)) eq 0 do freeport += 1
    
    return, freeport

end
