; docformat = 'rst'

;+
; Rotate and shift an image using linear or cubic interpolation.
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
;    The rotated and shifted image.
; 
; :Params:
; 
;    arr : in, type=array
;
;       The image array to be rotated and shifted.
;
;    angle : in, type=float
;
;       The rotation angle in radians.
;
;    dx : in, type=float
;
;       The X shift in pixels.
;
;    dy : in, type=float
;
;       The Y shift in pixels.       
; 
; 
; :Keywords:
; 
;    background : in, optional, type=float, default="median(arr)"
;   
;       Padding for non-data pixels.
;   
;    cubic : in, optional, type=float
;   
;       Cubic interpolation parameter, best value is -0.5.
;   
;    linear : in, optional, type=boolean
;   
;       Use linear interpolation, takes precedence over using the
;       cubic keyword.
;   
; 
; 
; :History:
; 
;    2020-06-22 : MGL. First version.
; 
;-
function red_rotshift, arr, angle, dx, dy $
                       , background = background $
                       , cubic = cubic $
                       , linear = linear

  if n_elements(background) eq 0 then background = median(arr)

  dim = size(arr, /dim)
  xsi = dim[0] * 0.5
  ysi = dim[1] * 0.5

  xgrid = findgen(dim[0]) # (fltarr(dim[1]) + 1.0)
  ygrid = (fltarr(dim[0]) + 1.0) # findgen(dim[1])

  xgrid *= float(dim[0]) / max(xgrid)
  ygrid *= float(dim[1]) / max(ygrid)

  xgrid1 = cos(angle[0]) * (xgrid - xsi - dx) - sin(angle[0]) * (ygrid - ysi - dy) + xsi
  ygrid1 = sin(angle[0]) * (xgrid - xsi - dx) + cos(angle[0]) * (ygrid - ysi - dy) + ysi

  ;; Interpolation onto new grid
  if keyword_set(linear) then begin
    return, interpolate(arr, xgrid1, ygrid1, missing = background)
  endif else begin
    return, interpolate(arr, xgrid1, ygrid1, missing = background, cubic = cubic)
  endelse

end
