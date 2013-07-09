FUNCTION red_m_retard, dd, aa, DEGREE=deg

IF n_params() EQ 1 THEN aa = 0
IF keyword_set(deg) THEN a = 2*aa*!dtor ELSE a = 2*aa
IF keyword_set(deg) THEN d = dd*!dtor ELSE d = dd

sa = sin(a)
sa2 = sa*sa
ca = cos(a)
ca2 = ca*ca
sd = sin(d)
cd = cos(d)

Re = [[1.0, 0.0,          0.0,           0.0], $
      [0.0, ca2+sa2*cd,   sa*ca*(1-cd), -sa*sd], $
      [0.0, sa*ca*(1-cd), sa2+ca2*cd,    ca*sd], $
      [0.0, sa*sd,       -ca*sd,         cd]]
return, Re
END
