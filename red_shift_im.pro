; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :Author:
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
;-
function red_shift_im, var, dx, dy, cubic = cubic, missing = missing

  if n_elements(cubic) eq 0 then cubic = -0.5
  if n_elements(missing) eq 0 then missing = median(var)
  
  p = [-dx, 0., 1., 0.] & q = [-dy, 1., 0., 0.]
  
  return, poly_2d(var, p, q, 2, cubic = cubic, missing = missing)

end
