FUNCTION red_polcal_sim, qwp, lp, par, EXTINCTION=elp, STOKES_IN=sp0


;;;  par = [ mm(*), d_qwp, delta_qwp ]

IF NOT keyword_set(elp) THEN elp = 0.
IF NOT keyword_set(sp0) THEN sp0 = [1., 0., 0., 0.]

mm = reform(par[0:15], 4, 4)

sp = red_m_retard(par[16], qwp+par[17], /deg) ## red_m_rot(-lp, /deg) ## $
     red_m_linpol(elp) ## red_m_rot(lp, /deg) ## sp0

; The first rotation does nothing, therefore simpler:

;sp = m_retard(par(16), qwp+par(17), /deg) ## m_rot(-lp, /deg) ## m_linpol(elp) ## sp0

return, mm ## sp

END
