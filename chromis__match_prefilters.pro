function chromis::match_prefilters, pf1, pf2
;                          NB     WB
    prefilter_table = [['3999', '3950'], $
                       ['3969', '3950'], $
                       ['3978', '3950'], $
                       ['3934', '3950'], $
                       ['3925', '3950'], $
                       ['4862', '4846']]

    idx = where( (prefilter_table[0,*] eq pf1 and prefilter_table[1,*] eq pf2) $
              or (prefilter_table[1,*] eq pf1 and prefilter_table[0,*] eq pf2) $
              or (prefilter_table[1,*] eq pf1 and prefilter_table[1,*] eq pf2) )    ; this row will match the 2 wideband channels

    if max(idx) gt 0 then return, 1
    
    return, 0

end
