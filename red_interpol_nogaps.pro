; docformat = 'rst'

;+
; Do interpolation with interpol(), but only where the input data do
; not have gaps. No extrapolation into the gaps.
;
; Parameters and keywords as for interpol() using the three parameters
; form.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;    Mats Löfdahl, Institute for Solar Physics
; 
; 
; :Returns:
; 
;    The interpolated Y values, NaN for intervals with no input data.
;
; :Keywords:
;
;    tol : in, optional, type=float, default=0.1
;
;       "Gaps" are defined as distances between x points that are
;       larger than dx*(1+tol), where dx is the median distance
;       between points.
; 
; :History:
; 
;    2017-04-24 : MGL. First version.
; 
;    2017-05-31 : MGL. New keyword tol.
; 
;    2017-06-08 : MGL. Improve the gap finding.
; 
; 
; 
; 
;-
function red_interpol_nogaps, y, x, xp, tol = tol, test = test, _ref_extra = extra

  if n_elements(tol) eq 0 then tol = 0.1
  
  ;; Output variable, filled with NaNs for now.
  yp = fltarr(n_elements(xp))
  yp[*] = !Values.F_NaN

  ;; Distance between data points.
  dx = red_differential(x)
  mdx = median(dx)

  ;; Find endpoints of OK intervals. This is based on the assumption
  ;; that the data have more or less constant dx, except for gaps.

  indx_largedist = where(dx gt mdx*(1.+tol), Nlarge)
  if Nlarge eq 0 then begin
    ;; No large distances detected, just do normal interpolation.
    return, interpol(y, x, xp, _strict_extra = extra)
  endif

  
  ;; Some large distances detected. We want to use only the
  ;; stretches of data where the distances are small enough.
  mask_smalldist = dx lt median(dx)*(1.+tol)
  onoff = red_differential(float(dx lt median(dx)*(1.+tol)))
;  on = [0, where(onoff gt 0)]
;  off = where(onoff lt 0)
  ;; Looks like the on and off points were off by one.
  on = [0, where(onoff gt 0)-1]
  off = where(onoff lt 0)-1
  if n_elements(off) lt n_elements(on) then off = [off, n_elements(x)-1]
  
  Nintervals = n_elements(on)
;  
;  indx_smalldist = where(dx lt median(dx)*(1.+tol), Nsmall)
;  mask_gaps =  deriv(indx_smalldist) gt 1
;  indx_gaps = [0, where(deriv(indx_smalldist) gt 1, Ngaps)]
;    
;  , n_elements(x)-1]
;  if count eq 0 then return, yp
;  dindx = where(deriv(indx) ne 1., Nd)
;  intervals = indx[dindx]
;
;  if red_odd(n_elements(intervals)) then stop

  for iinterval = 0, Nintervals-1 do begin
;    print, iinterval
    Npoints = off[iinterval] - on[iinterval]
    if Npoints gt 0 then begin
      ;;  indx = where((xp ge x[on[iinterval]]) and (xp le
      ;;  x[off[iinterval]-1]), count)
      ;; Include a margin of half the median dx in the data intervals.
      indx = where((xp ge x[on[iinterval]]-mdx/2.) and (xp le x[off[iinterval]]+mdx/2.), count)
      if count gt 0 then begin
        yp[indx] = interpol(y[on[iinterval]:off[iinterval]-1] $
                            , x[on[iinterval]:off[iinterval]-1] $
                            , xp[indx], _strict_extra = extra)
        if keyword_set(test) then begin
          print, 'Interpolation'
          print, 'X coordinates : ', x[on[iinterval]:off[iinterval]-1]
          print, 'Y values : ', y[on[iinterval]:off[iinterval]-1]
          print, 'Xp coordinates : ', xp[indx]
          print, 'Yp values : ', yp[indx] 
          stop
        endif
      endif
    endif 

  endfor                        ; iinterval


;  cgplot, x/3600, y, psym = 9, xrange = [10.75, 11.50]
;  cgplot, /over, xp/3600, yp, psym = 16
;  for i = 0, n_elements(xp) -1 do cgplot, /over, xp[i]/3600*[1, 1], !y.crange
;  indx = where(~finite(yp))
;  cgplot, /over, xp[indx]/3600,  0*indx, psym = 16, color = 'red'
;  
;  print, xp[indx]/3600
;  print, yp[indx]
;  stop
  
  return, yp
  
end


cd, '/scratch/alex/2017.09.06/CHROMIS'
a = chromisred(/dev)
a->fit_wb_diskcenter,pref='3950',tmax='12:00'


stop

x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 14, 16, 17, 18, 19]/20.*!pi
y = sin(x)

xp = findgen(200)/200*!pi

yp = interpol(y, x, xp, /quad)
yp1 = red_interpol_nogaps(y, x, xp, /quad)
yp2 = red_interpol_nogaps(y, x, xp, /quad, tol = .6)

window, 1
cgplot, x, y, psym = 16
cgplot, /over, xp, yp, psym = 1
cgplot, /over, xp, yp1, color = 'red', psym = 9

window, 2
cgplot, x, y, psym = 16
cgplot, /over, xp, yp, psym = 1
cgplot, /over, xp, yp2, color = 'red', psym = 9


end



