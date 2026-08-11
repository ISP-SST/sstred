; docformat = 'rst'

;+
; Computes analytical solution to a parabolic fit.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Jaime de la Cruz Rodriguez (ITA-UiO 2011)
; 
; 
; :Returns:
; 
;    Returns the 3 coefficients of the parabola: y = cf[0] + cf[1] * x
;    + cf[2] * x^2
;
; :Params:
; 
;    x : in, type=array
; 
;      X values to fit.
; 
;    y : in, type=array
;
;      Y values to fit.
; 
; :History:
; 
;    2017-11-16 : MGL. Import parab_fit.pro from Jaime.
; 
; 
;-
function red_parab_fit, x, y
  
  cf = fltarr(3)
  
  d = x[0]
  e = x[1]
  f = x[2]

  yd = y[0]
  ye = y[1]
  yf = y[2]
  
  cf[1] = ((yf - yd) - (f^2. - d^2.) * ((ye - yd) / (e^2. - d^2.)))/$
          ((f - d) - (f^2. - d^2.) * ((e - d) / (e^2. - d^2.)))
  
  cf[2] = ((ye - yd) - cf[1] * (e - d)) / (e^2.0 - d^2.0)
  
  cf[0] = yd - cf[1] * d - cf[2] * d^2.0
  
  return, cf

end
