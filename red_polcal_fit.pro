; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    data_r : 
;   
;   
;   
;    qq : 
;   
;   
;   
;    lp : 
;   
;   
;   
; 
; :Keywords:
; 
;    NORM : 
;   
;   
;   
;    init : 
;   
;   
;   
;    INI_NORM : 
;   
;   
;   
;    CHISQR : 
;   
;   
;   
;    FIX_QL : 
;   
;   
;   
;    PAR1 : 
;   
;   
;   
;    EXTINCTION : 
;   
;   
;   
;    STOKES_IN : 
;   
;   
;   
;    DOUBLE : 
;   
;   
;   
; 
; 
; :history:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
; 
;-
FUNCTION red_polcal_fit, data_r, qq, lp, NORM=norm, init=par0, INI_NORM=inrm, $
                         CHISQR=chi, FIX_QL=fql, PAR1=par1, EXTINCTION=elp, $
                         STOKES_IN=sp0, DOUBLE=dbl

  nq = n_elements(qq)
  nl = n_elements(lp)

  IF NOT keyword_set(norm) THEN norm = 1
  IF NOT keyword_set(elp) THEN elp = 0.
  IF NOT keyword_set(sp0) THEN sp0 = [1., 0., 0., 0.]

 ;;; determine data type - mpfit is very picky!
  IF keyword_set(dbl) THEN dtyp = 5 ELSE $
     dtyp = max([size(data_r, /type), size(qq, /type), $
                 size(lp, /type), size(par0, /type)])

  one = (make_array(1, type=dtyp, val=1))(0)
  null = 0*one

 ;;; first simple max normalization
;d1 = float(data_r)
  d1 = one*data_r
  FOR i=0, nl-1 DO d1(*, *, i) /= max(d1(*, *, i))

 ;;; set up fitting
  parinfo = replicate({fixed:0b, limited: [0b, 0b], limits: [null, null]}, 18)
  parinfo(0:15).limits(0) = (norm GT 1 ? -1.2*one : -one)
  parinfo(0:15).limits(1) = (norm GT 1 ?  1.2*one :  one)
  parinfo(0:15).limited = replicate(1b, 2, 16)

  IF keyword_set(par0) THEN BEGIN
     par = par0
     IF keyword_set(inrm) THEN BEGIN
        dmm = invert(reform(par(0:15), 4, 4))
        FOR j=0, nl-1 DO BEGIN
           FOR i=0, nq-1 DO BEGIN 
              c = 2*total(dmm(*, 0)*reform(d1(*, i, j)))
              d1(*, i, j) /= c 
           ENDFOR
        ENDFOR
     ENDIF
  ENDIF ELSE $
     par = [replicate(null, 16), 90., 0]

  IF keyword_set(fql) THEN BEGIN
      ;;; /fixed.  use to fix parameter in INIT
     parinfo(16:17).fixed = 1
      ;;; 2-element vector:  values for ret and off
     IF n_elements(fql) GE 2 THEN par(16:17) = fql(0:1)
      ;;; 4-element:  also specify which to fix, [ret,off,f_ret,f_off]
     IF n_elements(fql) EQ 4 THEN parinfo(16:17).fixed = fql(2:3)
  ENDIF

  xx = {qq:qq*one, lp:lp*one}
  f_args={X:xx, y:d1}

  ;;; initial fit
  par1 = mpfit('red_mpfunc_polcal', par, functargs=f_args, parinfo=parinfo, /quiet)

 ;;; now there is a reasonable estimate of the MM/DMM, and we can switch to
 ;;; using the normalisation within the evaluation function

  f_args={X:xx, y:one*data_r, norm: norm, EXTINCTION: elp, STOKES_IN: sp0}

  par2 = mpfit('red_mpfunc_polcal', par1, functargs=f_args, $
               parinfo=parinfo, best=chi, /quiet)

 ;;; final normalization

  CASE norm OF 
     1:  par2(0:15) /= max(par2(0:15))
     2:  par2(0:15) /= total(par2(indgen(4)*4))/4.
     3:  par2(0:15) /= rebin(par2(indgen(4)*4), 16, /sam)
     4:  par2(0:15) *= total((invert(reform(par2(0:15),4,4)))(*,0))
  ENDCASE

  return, par2

END
