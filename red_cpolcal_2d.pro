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
;    dat : 
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
;    guess : 
;   
;   
;   
; 
; :Keywords:
; 
;    chisq  : 
;   
;   
;   
;    NTHREADS : 
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
FUNCTION red_cpolcal_2D, dat, qq, lp, guess, chisq = chisq, NTHREADS=nt
;+
; NAME:
;       RED_CPOLCAL_2D
; PURPOSE:
;       Determine the 2D polarimetric calibration map of CRISP
; CALLING SEQUENCE:
;       ModMat = red_cpolcal_2d(data, qq, lp [, guess] )
; INPUTS:
;       data:  (float) array [4, N_q, N_lp, Sx, Sy] polcal data:
;              4 LC states, N_q angles of the Quarter wave plate (QWP), 
;              N_lp angles of the linear polarizer (LP)
;       qq:    (float) array [N_q] with the actual QWP angles
;       lp:    (float) array [N_lp] with the actual LP angles
; OPTIONAL PARAMETERS:
;       guess: 18-element array with initial values.  If omitted, they are
;              computed from the spatially averaged data
; KEYWORDS:
;       chisq: named variable to pass back the chi-square values of the fit
;       NTHREADS: number of threads to use [def:# of CPUs]
; OUTPUTS:
;       ModMat: (float) array [4, 4, Sx, Sy] with the fitted modulation matrix
;               (MM) in each pixel
; RESTRICTIONS:
;       The external C library cpolcal_simple.so has to be available.
; PROCEDURE:
;       Performs a Levenberg-Marquardt fit for the optimum MM of the CRISP
;       instrument.  The actual computation is performed in an external C
;       routine that uses the C version of MPFIT by Craig Markward.
;       The external code is based on a routine by Jaime de la Cruz
; MODIFICATION HISTORY:
;       20-Jan-2012  P.Suetterlin, ISP
;-

on_error, 2

IF n_params() LT 3 THEN BEGIN
    message, 'Usage: res = red_cpolcal_2d(data, qq, lp [, guess] )', /info
    return, -1
ENDIF

  ;;; locate the library.  Either in the same dir as the procedure, in the DLM
  ;;; path or the normal search path 
libname = 'cpolcal_simple.so'
dir = file_dirname((routine_info('cpolcal_2d',/func,/source)).path)
libfile = file_search(dir, libname)
libfile = file_search(strsplit(!DLM_PATH, ':', /extr), libname)
IF libfile EQ '' THEN $
  libfile = file_search(strsplit(!PATH, ':', /extr), libname)
IF libfile EQ '' THEN $
  message, 'Could not locate library file cpolcal_simple.so; Exiting'
libfile = libfile(0)

  ;;; Need parameter type checking, else things might (will) break
IF (size(dat, /type) NE 4) OR $
   (size(qq, /type) NE 4) OR $
   (size(lp, /type) NE 4) THEN BEGIN
    message, 'Input data DAT, QQ and LP have to be type FLOAT', /info
    return, -1
ENDIF
  ;;; # of threads = # of CPUs
IF NOT keyword_set(nt) THEN nt = !CPU.HW_NCPU ELSE nt = long(nt)

 ;;; If no initial guess is supplied, create one from averaged data
IF n_params() LT 4 THEN BEGIN
    d1d = total(total(dat, 5), 4)
    guess = polcal_fit(d1d, qq, lp, norm=4)
ENDIF

nqq = n_elements(qq)
nlp = n_elements(lp)
lc = bindgen(4)
nlc = 4

dim = size(dat, /dimension)
nx = dim[3]
ny = dim[4]
npix = nx * ny
res = dblarr(18, nx, ny)
FOR ii=0, 17 DO res[ii, *, *] = guess[ii]
chisq = fltarr(nx, ny)


dum = call_external(libfile, 'cpolcal_2d', nlc, nqq, nlp, npix, $
                    dat, qq, lp, res, chisq, nt)
  ;
res = float(reform(res(0:15, *, *), 4, 4, nx, ny))

return, res

END

