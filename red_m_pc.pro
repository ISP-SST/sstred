;       2013-07-12 : MGL. Renamed to red_m_pc for inclusion in
;                    crispred pipeline.
;
FUNCTION red_m_pc, mr, sa, ca, sd, cd, sc, cc
q = cc
u = sc
a = ca*ca+sa*sa*cd
b = sa*sa+ca*ca*cd
c = ca*sa*(1.-cd)
return, 0.5*(mr(0)+mr(1)*(a*q+c*u)+mr(2)*(c*q+b*u)+mr(3)*(sa*q-ca*u)*sd)
END
