function chromis::match_prefilters, pf1, pf2
;                          NB           WB
    prefilter_table = [['CaH-cont', 'CaHK-cont'], $
                       ['CaH-core', 'CaHK-cont'], $
                       ['CaH-red',  'CaHK-cont'], $
                       ['CaK-core', 'CaHK-cont'], $
                       ['CaK-blue', 'CaHK-cont'], $
                       ['Hb-core',  'Hb-cont']]

    idx = where( (prefilter_table[0,*] eq pf1 and prefilter_table[1,*] eq pf2) $
              or (prefilter_table[1,*] eq pf1 and prefilter_table[0,*] eq pf2) $
              or (prefilter_table[1,*] eq pf1 and prefilter_table[1,*] eq pf2) )    ; this row will match the 2 wideband channels

    if max(idx) gt 0 then return, 1
    
    return, 0

end
