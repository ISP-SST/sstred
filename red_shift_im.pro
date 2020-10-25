; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author: J. de la Cruz Rodriguez, Institute for Solar Physics
; 
; 
; 
; 
; :Returns:
; 
; 
; :Params:
; 
;    var : 
;   
;   
;   
;    dx : 
;   
;   
;   
;    dy : 
;   
;   
;   
; 
; :Keywords:
; 
;    cubic : in, optional, type=float, default=-0.5
;
;      Keyword to determine the kind of cubic interpolation done.
;   
;    missing : in, optional, type=float, default=median
;
;      Set pixels that are shifted in from outside the array to this
;      value.
; 
; :History:
; 
;   2013-06-04 : Split from monolithic version of crispred.pro.
; 
;   2014-01-13 : PS  Use (faster) Poly_2D
; 
;   2018-05-29 : MGL. New keyword missing.
;
;   2020-10-01 : JdlCR. Changed the routine to use internal bilinear
;                or nearest neighbor interpolation.
;
;-
function red_shift_im, var, dx, dy, cubic = cubic, missing = missing, nearest = nearest, nthreads = nthreads

;;  if n_elements(cubic) eq 0 then cubic = -0.5
  if n_elements(missing) eq 0 then missing = median(var)
  
;;  p = [-dx, 0., 1., 0.] & q = [-dy, 1., 0., 0.]
  
;;  return, poly_2d(var, p, q, 2, cubic = cubic, missing = missing)

  dim  = size(var, /dim)
  print, dim[1]
  xx = dindgen(dim[0])
  yy = dindgen(dim[1])

  xgrid = xx # (dblarr(dim[1]) + 1.0d) - dx
  ygrid = (dblarr(dim[0]) + 1.0d) # yy - dy

  ;; interpolate
  
  res = red_interpolate2D(xx, yy, var, xgrid, ygrid, nthreads = nthreads, nearest = nearest)

  ;; Find missing values
  
  xx -= dx
  yy -= dy
  
  xidx = where((xx gt dim[0]-1) or (xx lt 0), count)
  if(count gt 0) then res[xidx,*] = missing
  
  yidx = where((yy gt dim[1]-1) or (yy lt 0), count)
  if(count gt 0) then res[*,yidx] = missing

  return, temporary(res)
end
