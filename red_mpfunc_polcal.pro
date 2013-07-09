; docformat = 'rst'

FUNCTION red_mpfunc_polcal, par, X=x, Y=y, NORM=norm,  EXTINCTION=elp, STOKES_IN=sp0

  sq = size(X.qq, /dim)
  sl = size(X.lp, /dim)

  sim = make_array(4, sq(0), sl(0), type=size(par(0)*y(0), /type))

  FOR q=0, sq(0)-1 DO BEGIN
     FOR l=0, sl(0)-1 DO BEGIN
        sp = red_polcal_sim(X.qq(q), X.lp(l), par, EXTINCTION=elp, STOKES_IN=sp0)
        sim(*, q, l) = sp
     ENDFOR
  ENDFOR
   ;;; with the variable extinction and stokes input I also need to
   ;;; normalize the simulations, as they contain intensity variations
   ;;; that are corrected for in the raw data...
  IF keyword_set(elp) THEN sim = red_polcal_norm(sim, par)

  IF keyword_set(norm) THEN BEGIN
     y1 = y
     CASE norm OF 
        1:  nrm = max(par(0:15))
        2:  nrm = total(par(indgen(4)*4))/4.
        3:  nrm = rebin(par(indgen(4)*4), 16, /sam)
        4:  nrm = 1./total((invert(reform(par(0:15), 4, 4)))(*, 0))
     ENDCASE
     dmm = invert(reform(par[0:15]/nrm, 4, 4))
     FOR j=0,sl(0)-1 DO BEGIN
        FOR i=0,sq(0)-1 DO BEGIN 
           c = 2*total(dmm(*, 0)*reform(y1(*, i, j)))
           y1(*, i, j) /= c
        ENDFOR
     ENDFOR
     return, (y1-sim)(*)
  ENDIF ELSE $
     return, (y-sim)(*)

END
