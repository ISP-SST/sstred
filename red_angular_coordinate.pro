;	   File: 	angle_f.ana
;	   Created:	Wed, Nov 30 1994, 13:18:28
;	   Author: 	mats@astro.su.se

; Returns a 2r x 2r matrix with the angular coordinates of each pixel 
; as measured from the origin in pixel (r,r). The r_c parameter is a 
; dummy.

; Ported from ANA.

FUNCTION red_angular_coordinate,r,r_c,xOffset=xOffset,yOffset=yOffset,rotangle=rotangle
						   
  if n_elements(xOffset) eq 0 then xOffset = 0
  if n_elements(yOffset) eq 0 then yOffset = 0
  if n_elements(rotangle) eq 0 then rotangle = 0.0
  
  x = red_horizontal_coordinate(r,r_c,xOffset=xOffset)             ; x coordinates
  y = transpose(red_horizontal_coordinate(r,r_c,xOffset=yOffset))  ; y coordinates
  
  IF rotangle THEN BEGIN
    xold = x
    yold = y
    x = xold*cos(rotangle) + yold*sin(rotangle)
    y = yold*cos(rotangle) - xold*sin(rotangle)
  END
  
  mask = abs(x) LT 1e-40
  x = x - mask*(x - 1e-40)
  ang = atan(y,x)
  
  x = 0					   ; Save space!
  y = 0
  
  Return, ang

END; Angle

