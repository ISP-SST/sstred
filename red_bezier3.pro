;; Bezier Spline Interpolants: red_Bezier3
;; 
;; Purpose: Interpolate arrays y(x) onto new grid of values xx, using
;;          Bezier-Spline interpolants (quadratic and cubic).
;;          The optional keyword /linear performs linear extrapolation
;;          at the points outside the domain of x, y. Otherwise, the
;;          boundary value is repeated.
;;
;;          yy = red_bezier3(x, y, xx, /linear)
;;
;; Dependencies: cent_deriv (included)
;;
;; References: Auer 2003
;;
;; Author: Jaime de la Cruz Rodriguez (IFA-UU 2012)
;;
;; Modifications:
;;         2012-09-24, JdlCR: Created!
;;
;; -
function cent_deriv, x, y
  n = n_elements(x)
  yp = dblarr(n)

  ;;
  ;; dx, dx1
  ;;
  dx = (x[1:n-2] - x[0:n-3])
  dx1 =  (x[2:*] - x[1:n-2])
  
  ;;
  ;; Derivatives from both sides (displaced 1 elements compared to
  ;; their position in yp)
  ;;
  der = (y[1:n-2] - y[0:n-3]) / dx
  der1 = (y[2:*] - y[1:n-2]) / dx1
  
  
  ;;
  ;; Where the slope does not change sign, use formula. Else the
  ;; derivative is zero
  ;;
  idx = where(der*der1 gt 0.0, count)
  if(count gt 0) then begin
     lambda = (1.d0 + dx1[idx] / (dx1[idx] + dx[idx])) / 3.d0
     yp[idx+1] = (der[idx] * der1[idx] ) / (lambda * der1[idx] + (1.d0 - lambda) * der[idx]) 
  endif
  
  ;;
  ;; Fill borders
  ;; 
  yp[0] = der[0]  
  yp[n-1] = der1[n_elements(der1)-1]

  return, yp
end
function red_bezier3, x, y, xx, linear = linear
  n = n_elements(x)
  res = dblarr(n_elements(xx))
  
  ;;
  ;; derivatives
  ;;
  yp = cent_deriv(x, y)

  ;;
  ;; loop 
  ;;
  dx = x[1:*] - x[0:n-2]

  for k = 0L, n - 2 do begin
     k1 = k + 1
     idx = where((xx GT x[k]) AND (xx LE x[k1]), count)
     if(count eq 0) then continue

     ;;
     ;; Control point
     ;; 
     kdx = dx[k] / 3.0d0
     cntr = y[k] + kdx * yp[k] 
     cntr1 = y[k1] - kdx * yp[k1]

     ;;
     ;; Bezier Interpolation
     ;;
     u = (xx[idx] - x[k]) / dx[k]
     u1 = 1.d0 - u
     res[idx] =  y[k] * u1*u1*u1 + y[k1] * u*u*u + $
                 3.d0 * cntr * u * u1*u1 + 3.0d0 * cntr1 * u*u * u1

  endfor

  ;;
  ;; Values outside range -> repeat last value or linear extrapolation (dangerous!)
  ;;
  idx = where(xx le x[0], count)
  idx1 = where(xx gt x[n-1], count1)

  if(keyword_set(linear)) then begin
     if(count gt 0) then begin
        a = yp[0]
        b = y[0] - a * x[0]
        res[idx] = a * x[0] + b
     endif
     if(count1 gt 0) then begin
        a = yp[n-1]
        b = y[n-1] - a * x[n-1]
        res[idx1] = a * x[n-1] + b
     endif
  endif else begin
     if(count gt 0) then res[idx] = y[0]
     if(count1 gt 0) then res[idx1] = y[n-1]
  endelse

  return, temporary(res)
end
