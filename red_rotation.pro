; docformat = 'rst'

;+
; 
; 
; :Categories:
;
;    CRISP pipeline
; 
; 
; :author:
; 
; 
; 
; 
; :returns:
; 
; 
; :Params:
; 
;    img : 
;   
;   
;   
;    angle : 
;   
;   
;   
;    sdx : 
;   
;   
;   
;    sdy : 
;   
;   
;   
; 
; :Keywords:
; 
;    linear : in, optional, type=boolean
;
;      Set this to use bilinear interpolation. Otherwise bicubic
;      interpolation with cubic= -0.5 is used.
; 
; :history:
; 
;    2013-06-04 : Split from monolithic version of crispred.pro.
; 
;    2013-06-10 : Made using bicubic interpolation with -0.5 the
;    default and introduced a new flag for linear interpolation. MGL
; 
; 
;-
function red_rotation, img, angle, sdx, sdy, linear = linear

  if n_elements(sdx) eq 0 then sdx = 0.0
  if n_elements(sdy) eq 0 then sdy = 0.0

  dim = size(img, /dim)

  ;; get the index of each matrix element

  xgrid = findgen(dim[0]) # (fltarr(dim[1]) + 1.0)
  ygrid = (fltarr(dim[0]) + 1.0) # findgen(dim[1])
  
  ;; get rotation of indexes

  xsi = dim[0] * 0.5
  ysi = dim[1] * 0.5
  dx = cos(angle) * (xgrid - xsi) - sin(angle) * (ygrid - ysi) + xsi - sdx
  dy = sin(angle) * (xgrid - xsi) + cos(angle) * (ygrid - ysi) + ysi - sdy
  
  ;; Interpolation onto new grid
  if keyword_set(linear) then begin
     return, interpolate(img, dx, dy, missing = median(img))
  endif else begin
     return, interpolate(img, dx, dy, missing = median(img), cubic = -0.5)
  endelse

end
