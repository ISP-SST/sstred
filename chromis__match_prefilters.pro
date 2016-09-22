function chromis::match_prefilters, pf1, pf2
;                         NB      WB
    prefilter_table = [['3999', '3950'], $
                       ['3969', '3950'], $
                       ['3978', '3950'], $
                       ['3934', '3950'], $
                       ['3925', '3950'], $
                       ['4862', '4846']]

    Npref = n_elements(pf1)
    if( n_elements(pf2) ne Npref ) then begin
        print, 'chromis::match_prefilters: prefilter lists must be of the same length.'
        return, 0
    endif
    
    ret = bytarr(Npref)
    for ip=0, Npref-1 do begin
        ret[ip] = max(where( (prefilter_table[0,*] eq pf1[ip] and prefilter_table[1,*] eq pf2[ip]) $
                          or (prefilter_table[1,*] eq pf1[ip] and prefilter_table[0,*] eq pf2[ip]) $
                          or (prefilter_table[1,*] eq pf1[ip] and prefilter_table[1,*] eq pf2[ip]) ) $
                  ) ge 0
    endfor

    if n_elements(ret) eq 1 then return, ret[0]
    
    return, ret

end
