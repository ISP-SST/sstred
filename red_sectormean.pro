;	   File: 	sectormean_f.ana
;	   Created:	<2002-08-15 16:50:26 mats>
;	   Author: 	mats@astro.su.se
;	   Time-stamp: <2014-05-15 09:59:21 mats>

; Returns angular averages of |u| from zero frequency to limfreq in na
; sectors of the half plane from -90 to 90 degrees.

; The returned result is ordered from -90 to +90 degrees.

; If <angles> scalar or not specified, then compute default based on
; number of angles and return in keyword.

; If <width> zero or not specified, then compute default based on
; number of angles and return in keyword.

; Number of angles is given by the length of the first dimension in
; <angles> (or the smaller of that and <na> if <na> given).

FUNCTION red_sectormean, u, limfreq, na, mask=mask, xOffset=xOffset, yOffset=yOffset, angles = angles, width = width
  
  v = u

  IF n_elements(mask) eq 0 THEN maskr = 0*v+1 else maskr = mask
  if n_elements(limfreq) eq 0 then limfreq = (size(u, /dim))[0]/2

  if (size(angles, /dim))[0] gt 0 then na = (size(angles, /dim))[0]
  if n_elements(na) eq 0 then na = 4
    
  if n_elements(angles) eq 0 then angles = (indgen(na)/float(na)*!pi)-!pi/2
  if (size(angles, /dim))[0] eq 0 then angles = (indgen(na)/float(na)*!pi)-!pi/2

  if n_elements(width) eq 0 then width = !pi/na
  if width eq 0 then width = !pi/na

;  print, angles,  angles*180/!pi

  c = (size(v, /dim))[0]/2 
  
  rAbs = fltarr(limfreq+1,na)
  
  r = round(red_radial_coordinate(c, 1, xOffset=xOffset, yOffset=yOffset))
  n = lonarr(limfreq+1, na)
  
  ang = red_angular_coordinate((size(v, /dim))[0]/2, 1, xOffset=xOffset, yOffset=yOffset) 
  
 
  for i=0,limfreq do begin
     for j=0,na-1 do begin
        ;;help, r eq i, ang gt a(j), ang le a(j+1), maskr ne 0
        ;indx = where((r eq i) and (ang gt angles(j)) and (ang le angles(j+1)) and (maskr ne 0))
        indx = where((r eq i) and (maskr ne 0) and (abs(ang-angles(j)) lt width/2.0))
        if (size(indx))[0] gt 0 then begin ; Number of dimensions
           n(i,j) = n_elements(indx)
           rAbs(i,j) = total(v(indx))
        endif
     endfor
  endfor 

  n = float(n + (n EQ 0))       ; prevent division by zero
  
  return, rabs/n
  
end                             ; SectorMean
