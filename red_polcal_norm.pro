FUNCTION red_polcal_norm, data_r, par

sd = size(data_r, /dim)

dmm = invert(reform(par(0:15), 4, 4))
d1 = data_r
FOR j=0,sd(2)-1 DO BEGIN
    FOR i=0,sd(1)-1 DO BEGIN
        c = 2*total(dmm(*, 0)*reform(d1(*, i, j)))
        d1(*, i, j) /= c 
    ENDFOR
ENDFOR

return, d1

END
