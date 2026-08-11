;       2013-07-12 : MGL. Renamed to red_m_free_mirror for inclusion
;                    in crispred pipeline.
;
FUNCTION red_m_free_mirror, R, dd, DEGREE=deg
IF keyword_set(deg) THEN d = dd*!dtor ELSE d = dd
r1 = 2*sqrt(R)
Mae = 0.5*[[1+R, 1-R,  0.0,        0.0], $
           [1-R, 1+R,  0.0,        0.0], $
           [0.0, 0.0, -r1*cos(d), -r1*sin(d)], $
           [0.0, 0.0,  r1*sin(d), -r1*cos(d)]]
return, Mae
END
