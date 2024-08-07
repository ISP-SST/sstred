; docformat = 'rst'

;+
; Determine the 2D polarimetric calibration map of CRISP
; 
; Performs a Levenberg-Marquardt fit for the optimum MM of the CRISP
; instrument. The actual computation is performed in an external C
; routine that uses the C version of MPFIT by Craig Markward. The
; external C code by Jaime de la Cruz Rodriguez (IFA-UU).
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
;    2012-01-20 : P. Suetterlin, ISP 
; 
; 
; 
; :returns:
; 
;       ModMat: (float) array [4, 4, Sx, Sy] with the fitted modulation matrix
;               (MM) in each pixel
; 
; :Params:
; 
;    dat : in, type="fltarr(4, N_q, N_lp, Sx, Sy)"
;   
;      polcal data, 4 LC states, N_q angles of the Quarter wave plate
;      (QWP), N_lp angles of the linear polarizer (LP).
;   
;    guess : in, optional, type=fltarr(18)
;   
;       Initial values. If omitted, they are computed from the
;       spatially averaged data
;   
;    lp : in, type=fltarr(N_lp)
;   
;       The actual LP angles
;   
;    qq : in, type=fltarr(N_q)
;   
;      The actual QWP angles
;   
; 
; :Keywords:
; 
;    chisq : out
;   
;      Named variable to pass back the chi-square values of the fit
;
;    mask : in, optional, type="bytarr(Sx,Sy)"
;
;    nthreads : in, optional, type=integer, default="# of CPUs"
;   
;      Number of threads to use.
;   
; :Restrictions:
;
;       The external C library creduc.so has to be available.
; 
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2013-07-11 : MGL. Use red_polcal_fit, not polcal_fit.
; 
;   2022-09-02 : PS. New keyword mask.
; 
;-
FUNCTION red_cpolcal_2D, dat, qq, lp, guess, chisq = chisq, nthreads=nt, mask=mask

  on_error, 2
  
  IF n_params() LT 3 THEN BEGIN
     message, 'Usage: res = red_cpolcal_2d(data, qq, lp [, guess] )', /info
     return, -1
  ENDIF
  
  libfile = red_libfile('creduc.so')
  
  ;;; Need parameter type checking, else things might (will) break
  IF (size(dat, /type) NE 4) OR $
     (size(qq, /type) NE 4) OR $
     (size(lp, /type) NE 4) THEN BEGIN
     message, 'Input data DAT, QQ and LP have to be type FLOAT', /info
     return, -1
  ENDIF
  ;;; # of threads = # of CPUs
  IF NOT keyword_set(nt) THEN nt = !CPU.HW_NCPU ELSE nt = long(nt)
  
  nqq = n_elements(qq)
  nlp = n_elements(lp)
  lc = bindgen(4)
  nlc = 4
  
  dim = size(dat, /dimension)
  nx = dim[3]
  ny = dim[4]
  npix = nx * ny
  res = dblarr(18, nx, ny)
  IF keyword_set(mask) THEN BEGIN
    IF total((size(mask, /dim) EQ [nx, ny])) NE 2 THEN BEGIN
      message, 'Supplied mask must have the same dimensions as the data!', /info
      return, -1
    ENDIF
    IF (size(mask, /type) NE 1) OR (max(mask) GT 1b) THEN BEGIN
      message, 'Supplied mask must be a logical BYTE array', /info
      return, -1
    ENDIF
    tm = total(mask GT 0)
  ENDIF
  
  ;;; If no initial guess is supplied, create one from averaged data
  IF n_params() LT 4 THEN BEGIN
    IF keyword_set(mask) THEN BEGIN
      d1d = fltarr(nlc, nqq, nlp)
      FOR k=0, nlp-1 DO FOR j=0, nqq-1 DO FOR i=0, nlc-1 DO $
        d1d[i, j, k] = total(reform(dat[i, j, k, *, *])*mask)/tm
    ENDIF ELSE BEGIN
      d1d = total(total(dat, 5), 4)
    ENDELSE
    guess = red_polcal_fit(d1d, qq, lp, norm=4)
    ;;print, reform(guess[0:15], 4, 4)
    ;;print, guess[16:*]
  ENDIF
  FOR ii=0, 17 DO res[ii, *, *] = guess[ii]
  chisq = fltarr(nx, ny)
  
  print, 'red_cpolcal_2d : using libfile -> '+libfile
  IF keyword_set(mask) THEN BEGIN
    
    dum = call_external(libfile, 'cpolcal_2d_mask', nlc, nqq, nlp, npix, $
                        dat, qq, lp, res, chisq, nt, mask)

  ENDIF ELSE BEGIN
    
    dum = call_external(libfile, 'cpolcal_2d', nlc, nqq, nlp, npix, $
                        dat, qq, lp, res, chisq, nt)
  ENDELSE
                                ;
  res = float(reform(res[0:15, *, *], 4, 4, nx, ny))
  
  return, res
  
END

