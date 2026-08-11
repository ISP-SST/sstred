;       2013-07-12 : MGL. Renamed to red_m_sst for inclusion in
;                    crispred pipeline.
;
FUNCTION red_m_sst, az, el, par, DEGREE=deg

e = el + (keyword_set(deg) ? 90 : !pi/2)

L = red_m_lens(par(0:4))
Re = red_m_rot(e, DEG=deg)
Mae = red_m_free_mirror(par(5), par(6), DEG=deg)
Ra = red_m_rot(az, DEG=deg)
Rfp = red_m_rot(par(11), DEG=deg)
Mf = red_m_free_mirror(par(7), par(8), DEG=deg)
Ms = red_m_free_mirror(par(9), par(10), DEG=deg)
Rfm = red_m_rot(-par(11), DEG=deg)
Rlp = red_m_rot(par(13), DEG=deg)
Roa = red_m_rot(par(14), DEG=deg)

;; Thats like implemented by Michiel - wrong order of L and Rlp?
;;MM = Roa##Rfp##Ms##Mf##Rfm##Ra##Mae##Re##Mae##Rlp##L
MM = Roa##Rfp##Ms##Mf##Rfm##Ra##Mae##Re##Mae##L##Rlp

return, MM

END
