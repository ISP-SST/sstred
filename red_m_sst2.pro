FUNCTION red_m_sst2, az, el, par, DEGREE=deg

 ;;; note: The sign of the (az/el) angles is physically wrong, they should go
 ;;; in negative.  But polcal is wrong, too, and to stay consistent......

 ;;; 20120222:  Polcal is correct now, so also change this here
;       2013-07-12 : MGL. Renamed to red_m_sst2 for inclusion in
;                    crispred pipeline.
;

e = -el + (keyword_set(deg) ? 90 : !pi/2)

L = red_m_lens(par(0:4))
Re = red_m_rot(e, DEG=deg)
Mae = red_m_free_mirror(par(5), par(6), DEG=deg)
Ra = red_m_rot(-az, DEG=deg)
Rfp = red_m_rot(par(11), DEG=deg)
Mf = red_m_free_mirror(par(7), par(8), DEG=deg)
Ms = red_m_free_mirror(par(9), par(10), DEG=deg)
Rfm = red_m_rot(-par(11), DEG=deg)
Roa = red_m_rot(par(14), DEG=deg)

MM = Roa##Rfp##Ms##Mf##Rfm##Ra##Mae##Re##Mae##L
 ;;; Hmm, the Schuppmann should be after rotating back?
;MM = Roa##Ms##Rfp##Mf##Rfm##Ra##Mae##Re##Mae##L

return, MM

END
