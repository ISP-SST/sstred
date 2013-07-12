;       2013-07-12 : MGL. Renamed to red_m_pbs for inclusion in
;                    crispred pipeline.
;
FUNCTION red_m_pbs, d
;;; polarizing beamsplitter.  d determines which beam
IF d EQ 0 THEN a = 1 ELSE a = -1
return, 0.5*[[1., a, 0, 0], [a, 1., 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]
END

  
