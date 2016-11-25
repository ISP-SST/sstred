; docformat = 'rst'

;+
; Weighted smoothing.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, ISP
; 
; 
; :Returns:
; 
;    A smoothed version of the y, taking the weights into account.
; 
; :Params:
; 
; 
; :Keywords:
; 
;   
;   
;   
; 
; 
; :History:
; 
;    2016-11-25 : MGL. First version.
; 
; 
;-

function red_wsmooth, x, y, window, ndeg, weight = weight, select_fraction = select_fraction

  if n_elements(select_fraction) eq 0 then select_fraction = 0.5
  
  ;; Smoothed value is based on fitting an ndeg degree polynomial to
  ;; the best half of the data points in the window, based on the
  ;; weights. 
  
  ;; This function is inspired by lowess.pro in astrolib.
  if n_elements(weight) eq 0 then stop

  m = round(window/2)           ; Prefer odd window!
  Npoints = n_elements(x)
  z = fltarr(Npoints)

  ;; For use with mpfitexpr:
  case ndeg of
     1: expr = 'P[0] + P[1]*X'
     2: expr = 'P[0] + P[1]*X + P[2]*X*X'
     else: stop
  endcase

  for ipoint = 0, Npoints-1 do begin

     case 1 of
        ipoint lt m/2           : indx = indgen(m)
        ipoint lt m             : indx = indgen(ipoint*2+1)
        ipoint gt Npoints-m/2-1 : indx = indgen(m) - m + Npoints 
        ipoint gt Npoints-m-1   : indx = indgen((Npoints-ipoint)*2) - (Npoints-ipoint)*2 + Npoints
        else                    : indx = indgen(window) + ipoint - m
     endcase
     
     Ninterval = n_elements(indx)

     u = x[indx] 
     v = y[indx]
     w = weight[indx]

     indx = (reverse(sort(w)))[0:round(Ninterval*select_fraction)-1]
     u = u[indx]
     v = v[indx]
     w = w[indx]

     u0 = x[ipoint] - mean(u)
     u = u-mean(u)
     
     p = mpfitexpr(expr, u, v, weights = w)

     case ndeg of
        1: yfit = p[0] + p[1]*u0
        2: yfit = p[0] + p[1]*u0 + p[2]*u0*u0
        else: stop
     endcase
     
     z[ipoint] = yfit
     
  endfor                        ; ipoint

  return, smooth(z, 5)

end
