; docformat = 'rst'

;+
; The 2QI method from Yi & Molowny Horas, page 72-73, see Löfdahl (2010).
; 
; :Categories:
;
;    SST observations
; 
; 
; :author:
; 
;    Mats Löfdahl, 2011
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    c : in, type=array
;
;      The 2D array in which we want to find the peak to subpixel
;      accuracy. If c is a 3x3 array, we assume the center pixel has
;      the largest value. If c is larger, we crop it to such a 3x3
;      array.
;
; :history:
; 
;    2013-09-11 : MGL. Renamed to red_findmax2qi for inclusion in
;                 crispred pipeline.
; 
; 
;-
function red_findmax2qi, c, verbose = verbose

  cdim = size(c, /dim)

  if min(cdim) lt 3 then begin
     print,'findmax2QI : Input array too small in at least one dimension' 
     help, c
     retall
  endif

  ;; Whole pixel max point
  mx = max(c, lc)
  xp = lc mod cdim(0)
  yp = lc / cdim(0) 

  ;; Peak not surrounded by smaller values?
  csdim = size(cs, /dim)
  if xp eq 0 or xp eq csdim(0)-1 or yp eq 0 or yp eq csdim(1)-1 then begin
     ;; We could choose to return the max position anyway, but for now
     ;; we do not.
     print,'findmax2QI : Peak on border'
     retall
  endif

  if max(cdim) gt 3 then begin
     if keyword_set(verbose) then print,'findmax2QI : larger than 3x3, cropping'
     cs = c[xp-1:xp+1, yp-1:yp+1]
  endif else cs = c
  
  ;; 2D parabolic fit, b array is coefficients in this fit:
  ;; b0 + b1 X + b2 Y + b3 XY + b4 XX + b5 YY.
  ;; Fit is done explicitly.
  b = dblarr(6)

  b(1) = (cs( 2, 1) - cs( 0, 1))/2.
  b(2) = (cs( 1, 2) - cs( 1, 0))/2.
  b(3) = (cs( 2, 2) - cs( 2, 0) - cs( 0, 2) + cs( 0, 0))/4.
  b(4) = (cs( 2, 1) - 2*cs( 1, 1) + cs( 0, 1))/2.
  b(5) = (cs( 1, 2) - 2*cs( 1, 1) + cs( 1, 0))/2.

  xmax = (2*b(1)*b(5)-b(2)*b(3))/(b(3)^2-4*b(4)*b(5))
  ymax = (2*b(2)*b(4)-b(1)*b(3))/(b(3)^2-4*b(4)*b(5))

;;  tvimg, cs, xp+[-1.5, 0, 1.5], yp+[-1.5, 0, 1.5], /sample
;;  oplot, [xmax+xp], [ymax+yp], color=fsc_color('red'), psym = symcat(9)
;;  print, [xmax,ymax] + [xp, yp]
;;
;;  plot, xp+[-1, 0, 1], cs[*, 1]
;;  oplot, [xmax+xp], [mean(cs[*, 1])], color=fsc_color('red'), psym = symcat(9)
;;
;;  plot, yp+[-1, 0, 1], cs[1, *]
;;  oplot, [ymax+yp], [mean(cs[1, *])], color=fsc_color('red'), psym = symcat(9)

  return, [xmax,ymax] + [xp, yp]

end; findmax2QI
