PRO red_clean_warp, map, kx, ky, ord, FACTOR=fak, NITER=nit, DOUBLE=double, $
                    VERBOSE = verb, REFSTAT = rs
;; wrapper around polywarp that checks for outliers that don't properly
;; match the model, and excludes them from the fit

IF n_params() LT 4 THEN ord = 2
IF NOT keyword_set(nit) THEN nit = 3
IF NOT keyword_set(fak) THEN fak = 3.
IF n_elements(double) EQ 0 THEN double = 1

ni = n_elements(map)/4
dmap = reform(map, ni, 4)
xin = dmap(*, 0)
yin = dmap(*, 1)
xou = dmap(*, 2)
you = dmap(*, 3)

IF keyword_set(rs) THEN BEGIN
    mask = (rs(*, *, 0)/rs(*, *, 1))(*)
    ix = where(mask GT median(mask))
ENDIF ELSE $
  ix = indgen(ni)
    
FOR it=1, nit DO BEGIN
    polywarp, xin(ix), yin(ix), xou(ix), you(ix), ord, kx, ky, DOUBLE=double
    ;;; apply model - and use double due to the large numbers!
    xm = 0 & ym = 0
    FOR i=0, ord DO BEGIN
        yy = (double(you)^i)
        FOR j=0, ord DO BEGIN
            xx = (double(xou)^j)*yy
            xm += kx(i, j)*xx
            ym += ky(i, j)*xx
        ENDFOR
    ENDFOR
;    FOR i=0, ord DO FOR j=0, ord DO xm += kx(i, j)*(xou^j)*(you^i)
;    FOR i=0, ord DO FOR j=0, ord DO ym += ky(i, j)*(xou^j)*(you^i)
    dr = sqrt((xm-xin)^2. + (ym-yin)^2.)
    sr = stdev(dr)
    ix1 = where(dr LE fak*sr, ni1)
    IF ni1 GE ni THEN BREAK  ;;; no further improvement, exit
    ix = ix1
    ni = ni1
ENDFOR
IF keyword_set(verb) THEN $
  message, 'Used '+strtrim(ni1, 2)+' of '+strtrim(n_elements(map)/4, 2)+ $
           ' data points after '+strtrim(it, 2)+' iterations', /info
END
