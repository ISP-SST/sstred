;	   File: 	round_matrix_f.ana
;	   Created:	Wed, Nov 30 1994, 13:27:33
;	   Author: 	mats@astro.su.se

; Ported from ANA

; Returns a 2r x 2r matrix with the values of <arr> rotated around the
; origin in pixel (r,r).



FUNCTION red_RoundMatrix,arr,r 
                           
  r_c = (size(arr,/dim))[0]          
  arr2 = fltarr(r+1)
  
  ;d,arr,arr2
  
  arr2(0) = arr(0:r<(r_c-1))
  
  xxx = ((findgen(2*r)-r)) # replicate(1.,2*r)
  rad = sqrt(xxx^2+transpose(xxx)^2)

  ;rad = Radius(r,r_c) * r_c			   ; The coordinates of arr
  sz = (size(rad,/dim))[0] 

  rad = reform(rad,sz*sz,/over)                    ; are taken to range from
  						   ; zero to r_c.

  Rnd_Mtrx = interpol(arr2,IndGen((size(arr2,/dim))[0]),rad)
  
  Rnd_Mtrx = reform(Rnd_Mtrx,sz,sz,/over)*(rad lt r_c);aperture(r,r_c)
  
  Return, Rnd_Mtrx

END; RoundMatrix
