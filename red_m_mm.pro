;       2013-07-12 : MGL. Renamed to red_m_mm for inclusion in
;                    crispred pipeline.
;
FUNCTION red_m_mm, sc, pbs=pbs
; return modulation matrix for a LC scheme
; scheme should be array [2,n] 
IF NOT keyword_set(pbs) THEN pbs = 0
s = size(sc, /dim)
ns = s(1)
res = fltarr(4, ns)
FOR i=0, ns-1 DO res(*, i) = $
  (2*m_pbs(pbs)##red_m_retard(sc(1, i),45.,/deg)##red_m_retard(sc(0, i),0.,/deg))(*,0)
return, res
END
