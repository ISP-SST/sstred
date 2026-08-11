; docformat = 'rst'

;+
; Shift an image with subpixel accuracy.
; 
; :Categories:
;
;    SST pipeline
; 
; 
; :Author:
; 
;     J. Chae
; 
; 
; :Returns:
; 
;    A shifted version of the input image.
; 
; 
; :Params:
; 
;    image : in
;
;       The image to be shifted.
;
;    x0 : in, type=float
;
;       The first-dimension shift in pixels.
;
;    y0 : in, type=float
;
;       The second-dimension shift in pixels.
; 
; 
; 
; :History:
; 
;     2012-05-04 : MGL. Use bicubic interpolation.
;
;     2016-11-29 : MGL. Renamed into the red_ namespace.
;   
; 
;-
function red_shift_sub, image, x0, y0

  if fix(x0)-x0 eq 0. and fix(y0)-y0 eq 0. then return, shift(image, x0, y0)
  
  s = size(image)
  x = findgen(s(1))#replicate(1., s(2))
  y = replicate(1., s(1))#findgen(s(2))

  x1 = (x-x0)>0<(s(1)-1.)  
  y1 = (y-y0)>0<(s(2)-1.)  

  return, interpolate(image, x1, y1, cubic = -0.5)

end  
