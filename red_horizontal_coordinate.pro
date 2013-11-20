;	   File: 	X_Coord_F.ANA
;	   Created:	Tue, Sep 3 1991, 17:43:09
;	   Author: 	LOFDAHL@royacs.astro.su.se

; Return 2D matrix with x coordinates relative to origin at (r,r).          
; Coordinate is unity at r_c

; Rewrite this to use IndGen of array with second argument to give 
; x- or y- (0 or 1, resp.) coordinates.


FUNCTION red_horizontal_coordinate,r,r_c,xOffset=xOffset

  x = dindgen(2*r,2*r) mod (2*r)  ; 2D index in x-direction.

  x = x - r
  if n_elements(xOffset) ne 0 then x = x - xOffset
  x = x / float(r_c)

  Return, x 

END; X_Coord
