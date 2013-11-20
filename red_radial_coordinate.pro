;	   File: 	radius_f.ana
;	   Created:	Wed, Nov 30 1994, 13:14:52
;	   Author: 	mats@astro.su.se

; Returns a 2r x 2r matrix with the radial coordinates of each pixel 
; as measured from the origin in pixel (r,r). The value is unity at 
; the distance r_c from the origin.

; Ported from ANA

function red_radial_coordinate,r,r_c,xOffset=xOffset,yOffset=yOffset 

  if n_elements(xOffset) eq 0 then xOffset = 0
  if n_elements(yOffset) eq 0 then yOffset = 0
  
   x = red_horizontal_coordinate(r,r_c,xOffset=xOffset)
   y = transpose(red_horizontal_coordinate(r,r_c,xOffset=yOffset))

   rr = sqrt(x*x+y*y)

   x = 0					   ; Save space!
   y = 0

   return, rr

end
