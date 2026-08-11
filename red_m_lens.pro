;       2013-07-12 : MGL. Renamed to red_m_lens for inclusion in
;                    crispred pipeline.
;
FUNCTION red_m_lens, par
L = [[1.0, 0.0,     0.0,     0.0], $
     [0.0, par(0),  par(1), -par(2)], $
     [0.0, par(1),  par(3),  par(4)], $
     [0.0, par(2), -par(4),  par(0)+par(3)-1.0]]
return, L
END
