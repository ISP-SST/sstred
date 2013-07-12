;       2013-07-12 : MGL. Renamed to red_m_rot for inclusion in
;                    crispred pipeline.
;
FUNCTION red_m_rot, aa, DEGREE=deg

IF keyword_set(deg) THEN a = 2*aa*!dtor ELSE a = 2*aa

Re = [[1.0,  0.0,    0.0,    0.0], $
      [0.0,  cos(a), sin(a), 0.0], $
      [0.0, -sin(a), cos(a), 0.0], $
      [0.0,  0.0,    0.0,    1.0]]
return, Re
END
