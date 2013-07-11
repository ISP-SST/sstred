; docformat = 'rst'

;+
; Vectorized version adapted from Mats Carlsson routine intep.pro
; 
; Reference: Publications of the Dominion astrophysical observatory,
; xvi,6,67 graham hill: intep, an effective interpolation subroutine.
;
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
;    Jaime de la Cruz Rodriguez (ITA-UiO 2011).
; 
; 
; :returns:
;
;        Returns the values of the function interpolated
;        at the xp values.
; 
; :Params:
; 
;    x :
;
;       array to interpolate in
;
;    y :
;
;       array to interpolate in 
;
;    xp :
;
;       array where y values are wanted
; 
; :Keywords:
; 
;    linear : in, type=boolean
;
;        if present, it performs linear extrapolation at the
;        out-of-bound points. By default it repeats the bound values
;        of y.
;
;
; 
; 
; :history:
; 
;        2011-04-13 : Jaime. The loop in ii should be up to np (not
;                     np-1).
;
;        2011-09-17 : Jaime. i-1 forced to be >0. Also replaced: (xp
;                     LT x[i>0]) by (xp LE x[i>0])
;
;        2013-07-11 : MGL. Renamed to red_intepf and added to the
;                     crispred pipeline. Changed documentation to rst
;                     docinfo format.
;
;-
function red_intepf, x, y, xp, linear = linear ;, lp1 = lp1, lp2 = lp2, fp2 = fp2, fp1 = fp1
  ;
  np = n_elements(x)
  np1 = np - 1L
  nxp = n_elements(xp)
  ;
  yp = dblarr(nxp)
  lp1 = dblarr(np1)
  lp2 = dblarr(np1)
  fp1 = dblarr(np1)
  fp2 = dblarr(np1)
  ;
  ; Compute interpolation coefficients (vectorized)
  ; 
  lp1[0:np-2] = 1.0d0 / (x[0:np-2] - x[1:np-1])
  lp2[0:np-2] = 1.0d0 / (x[1:np-1] - x[0:np-2])
  ;
  fp1[1:np-2] = (y[2:np-1] - y[0:np-3]) / (x[2:np-1] - x[0:np-3])
  fp1[0] = (y[1] - y[0]) / (x[1] - x[0])
  ;
  fp2[np-2] = (y[np-1] - y[np-2]) / (x[np-1] - x[np-2])
  fp2[0:np-3] = fp1[1:np-2] ;(y[2:np-1] - y[0:np-3]) / (x[2:np-1] - x[0:np-3])
  ;
  ; Must do at least one loop to detect which value of x must be
  ; subtracted to xp in order to create xpi and xpi1
  ;
  k=0L
  for ii = 1L, np  do begin
     ;
     i = ii - 1L
     i1 = (i - 1L)>0
     pos = where((xp LE x[i>0]) AND (xp GE x[i1]), ct)
     ;
     if ct gt 0 then begin
        ;
        xpi = xp[pos] - x[i1]
        xpi1 = xp[pos] - x[i]

        ; Square it right away (when computing yp we 
        ; don't need to do l1*l1, l2*l2)
        l1 = (xpi1 * lp1[i1])^2.
        l2 = (xpi * lp2[i1])^2.

        ; Interpolated values
        yp[pos] = y[i1] * (1.0d0 - 2.0d0 * lp1[i1] * xpi) * l1 + $
                  y[i] * (1.0d0 - 2.0d0 * lp2[i1] * xpi1) * l2  + $
                  fp2[i1] * xpi1 * l2 + fp1[i1] * xpi * l1
        ;
     endif
     ;
  endfor
  ;
  ; Values out of bounds 
  ;
  if keyword_set(linear) then begin
     ; Linear Extrapolation
     pos = where(xp gt max(x), ct)
     if ct gt 0 then begin
        a = (y[np-1] - y[np-2]) / (x[np-1] - x[np-2])
        b = y[np-1] - a * x[np-1]
        ;
        yp[pos] = a * xp[pos] + b
     endif
     ;
     pos = where(xp lt min(x), ct)
     if ct gt 0 then begin
        a = (y[1] - y[0]) / (x[1] - x[0])
        b = y[0] - a * x[0]
        ;
        yp[pos] = a * xp[pos] + b
     endif
  endif else begin
     ; Repeat bound values (DEFAULT)
     pos = where(xp gt max(x), ct)
     if ct gt 0 then yp[pos] = y[np1]
     pos = where(xp lt min(x), ct)
     if ct gt 0 then yp[pos] = y[0]
     ;
  endelse
  ;
  return, yp
end

